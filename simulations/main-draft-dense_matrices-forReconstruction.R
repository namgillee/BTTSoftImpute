## main-draft-dense_matrices-forReconstruction.R
# Simulation for BTTSoftImpute with low-rank dense matrices
# 
# Goal: Compare reconstruction errors of the two methods: 
#   - BTTSoftImpute: The proposed method with tensor data 
#                   (by Namgil Lee & Jong-Min Kim, 2018)
#   - SoftImpute: The standard method with matrix data 
#                   (by Mazumder & Hastie, 2010, JMLR)
#
# Data: Generate dense matrices,
#         Y = USV + E
#     by dense U, S, V (V in block TT format), their sizes are
#         U : I x R
#         S : R x R
#         V : J*J*...*J x R
#     For BTT matrix V, 
#     - use low-order (small N), small-to-large sizes (size J_1,..., J_N)
#
#### NOTE: This simulation code can take a long time to run through. 
####    So please adjust experimental parameters in the biginning of the file. 
#
# Last modified at 2018.08.29. by Namgil Lee (Kangwon National University)
#
# < Regularized Block Tensor Train Decomposition >
# Copyright (C) 2018  Namgil Lee


## Read source code
rm(list=ls())
set.seed(11111)

file.sources = list.files("../softImpute_BlockTT/", pattern="*.R", full.names = TRUE)
sapply(file.sources,source,.GlobalEnv)

outdir = 'dense_matrices'
if (!file.exists(outdir)) {
  dir.create(outdir)
}


## Set Data Sizes
#
# 100 x J^N
#
ORDERS <- seq(2,3,by = 1)     #Tensor orders
sizeI  <- 100                 #Row size of data matrix
sizeJ  <- 12    # J=4:  Considering real image data quantization. 
rankRY <- 4     # RY=5: Number of singular values to estimate is not so small.
ttrankY <- 2    # R_n=4: TT-rank of V_Y is not large.

rankRX <- rankRY
myttrankmax <- 10     # R_n(Y) <= R_n(X) : TT-rank of V_X
mytol_df <- -1        # (unused)
mytol_dx <- 5e-4      # Small tol_dx  ---- sufficient accuracy;; 1e-3 / sqrt(5-1) = 5e-4
mytol_f  <- 0         # lambda=0  ==>  obj -> 0 
maxiter <- 100000     # Sufficiently large, because we expect tol_df works...

rate_obs = 0.4

noisestd <- 0   #noise standard deviation
lambda = 0

num_rep = 30  ##



time_btt = array(-1,c(num_rep,length(ORDERS))) #2: {training, test}
#conv_btt = NULL
time_mat = array(-1,c(num_rep,length(ORDERS))) #2: {training, test}
#conv_mat = NULL
error_btt = array(-1,c(num_rep,2,length(ORDERS))) #2: {training, test}
error_mat = array(-1,c(num_rep,2,length(ORDERS))) #2: {training, test}

## Main iteration
for (ido in 1:length(ORDERS)) 
{
  orderN <- ORDERS[ido]
  print(paste('----- Order:', orderN, '-----'))
  
  for (idrep in 1:num_rep) {    
      ## Create a random [Ix(J*J*...*J)] matrix, Y = UDV';  U~normal(0,1), D~uniform[0,1], V~ttRand(...,'normal')
      U0 = matrix(rnorm(sizeI * rankRY), nrow = sizeI, ncol = rankRY)
      U0 = svd(U0)$u
      d0 = sort(runif(n = rankRY, min = 0.5, max = 1), decreasing = TRUE) #* sqrt(sizeI*sizeJ^orderN/rankRY)
      my_block_tt_ranks = pmin( rankRY * sizeJ^(0:orderN), 
                                sizeJ^(orderN:0),
                                c(rankRY, rep(ttrankY,orderN-1), 1)  )   # block-TT-rank for V IS BOUNDED BY RTT_MAX, 
                                                 # And LR&RL orthogonalization is considered.
      V0 = bttRand(J = rep(sizeJ, orderN), N = orderN, R = my_block_tt_ranks, direction = -1)
      Y0 = list(u = U0, d = d0, v = V0)
      
      

      ## Select observed values randomly from the original matrix Y0. 
      # The number of observed values: (1-missP)*I*prod(J)
      # Y0_omega is in SparseMatrix format
      Y0_omega = select_observed_from_LR(x = Y0, rate = rate_obs)
      

      ## Temporary variable: 'Y0 orthogonal'
      Y0orth = Y0$u %*% diag(Y0$d) %*% t(bttFull(Y0$v)) 
      for ( idnz in 1:length(Y0_omega$v) ) { # We know nnz>0
        m = Y0_omega$i[idnz]   #row index
        j = Y0_omega$j[idnz, ] #column indices
        CumSize = cumprod(c(1,Y0_omega$ncol[1:(orderN-1)]))
        Y0orth[m,1+sum((j-1)*CumSize)] = 0
      }# (Y0 - Y0_omega)
      
            
      
      ## Add small Gaussian noise
      Y_omega = Y0_omega
      vdata = Y_omega$v
      vdata = vdata + noisestd * rnorm(length(vdata), 0, 1)
      Y_omega$v = vdata
      
    
      ## Conduct tensor completion by using BTT(by Lee&Kim) and standard thresholded SVD(by Mazumder&Hastie)
      ## 1) Using block-TT-ALS
      start_time = Sys.time()  
      Z_btt_als = bttdSoftImpute(Y = Y_omega, lambda = 0, 
                                 rank_Y = rankRX, ttrank_max = myttrankmax, 
                                 tol_dx = mytol_dx, tol_df = mytol_df, tol_f = mytol_f,
                                 verbose = TRUE, method = 'als', maxit = maxiter)
      end_time = Sys.time()  
      time_btt[idrep,ido] = difftime(end_time, start_time, units = 'sec')
      #conv_btt[[ido]] = attr(Z_btt_als,'obj')
      
      ## Compute Training Error
      Yres = sp_minus_lr_btt(Y_omega, Z_btt_als)
      error_btt[idrep,1,ido] = sqrt( sum(Yres$v^2) ) / sqrt( sum(Y_omega$v^2) ) 
      
      ## Compute Test Error
      ## (Y0 - R_omega) - X; Here, R_omega is the Yres defined above.
      Zbres = Y0$u %*% diag(Y0$d) %*% t(bttFull(Y0$v))  - 
        Z_btt_als$u %*% diag(Z_btt_als$d) %*% t(bttFull(Z_btt_als$v)) 
      for ( idnz in 1:length(Yres$v) ) { # We know nnz>0
       m = Yres$i[idnz]   #row index
       j = Yres$j[idnz, ] #column indices
       CumSize = cumprod(c(1,Yres$ncol[1:(orderN-1)]))
       Zbres[m,1+sum((j-1)*CumSize)] = Zbres[m,1+sum((j-1)*CumSize)] - Yres$v[idnz]
      }# (Y0 - R_omega)
      error_btt[idrep,2,ido] = sqrt(sum(Zbres^2)) / sqrt(sum(Y0orth^2))
      
      
      ## 2) Using standard thresholded SVD
      Ymat_omega = reduce_to_sparse_matrix(Y_omega)
      
      start_time = Sys.time()  
      Z_mat_als = bttdSoftImpute(Y = Ymat_omega, lambda = 0, 
                                 rank_Y = rankRX, ttrank_max = myttrankmax, 
                                 tol_dx = mytol_dx, tol_df = mytol_df, tol_f = mytol_f,
                                 verbose = TRUE, method = 'als', maxit = maxiter)
      end_time = Sys.time()  
      time_mat[idrep,ido] = difftime(end_time, start_time, units = 'sec')
      # conv_mat[[ido]] = attr(Z_mat_als,'obj')
      
      ## Compute Training Error
      Yres = sp_minus_lr_btt(Ymat_omega, Z_mat_als)
      error_mat[idrep,1,ido] = sqrt( sum(Yres$v^2) ) / sqrt( sum(Y_omega$v^2) ) 
      
      ## Compute Test Error
      ## (Y0 - R_omega) - X; Here, R_omega is the Yres defined above.
      Zbres = Y0$u %*% diag(Y0$d) %*% t(bttFull(Y0$v))  - 
        Z_mat_als$u %*% diag(Z_mat_als$d) %*% t(bttFull(Z_mat_als$v)) 
      for ( idnz in 1:length(Yres$v) ) { # We know nnz>0
        m = Yres$i[idnz]   #row index
        j = Yres$j[idnz]#j = Yres$j[idnz,] #column indices
        #CumSize = cumprod(c(1,Yres$ncol[1:(orderN-1)]))
        Zbres[m,j] = Zbres[m,j] - Yres$v[idnz]
      }# (Y0 - R_omega)
      error_mat[idrep,2,ido] = sqrt(sum(Zbres^2)) / sqrt(sum(Y0orth^2))
      
      
      print(paste('---- btt_trn, btt_tst, mat_trn, mat_tst: ', error_btt[idrep,1,ido], error_btt[idrep,2,ido], 
                  error_mat[idrep,1,ido], error_mat[idrep,2,ido],'----'))
      
      
  }  #end for num_rep
}#end for ORDERS



## Save the results
save(file = paste('resu_main-draft-dense_matrices-forReconstruction.Rdata' , sep=''), 
     time_btt, time_mat,
     error_btt, error_mat, 
     ORDERS
)


## A summary table
dif_error <- error_btt - error_mat
resu_table = matrix(0,6,4)
rownames(resu_table) <- c('N.2..BTT', 
                          'N.2..Matrix',
                          'N.2..D(BTT.Matrix)',
                          'N.3..BTT', 
                          'N.3..Matrix',
                          'N.3..D(BTT.Matrix)')
colnames(resu_table) <- c('Training.Error.Mean','Training.Error.SD','Test.Error.Mean','Test.Error.SD')

resu_table[1,] <- c( mean(error_btt[,1,1]) , sd(error_btt[,1,1]) , mean(error_btt[,2,1]) , sd(error_btt[,2,1]) )
resu_table[2,] <- c( mean(error_mat[,1,1]) , sd(error_mat[,1,1]) , mean(error_mat[,2,1]) , sd(error_mat[,2,1]) )
resu_table[3,] <- c( mean(dif_error[,1,1]) , sd(dif_error[,1,1]) , mean(dif_error[,2,1]) , sd(dif_error[,2,1]) )
resu_table[4,] <- c( mean(error_btt[,1,2]) , sd(error_btt[,1,2]) , mean(error_btt[,2,2]) , sd(error_btt[,2,2]) )
resu_table[5,] <- c( mean(error_mat[,1,2]) , sd(error_mat[,1,2]) , mean(error_mat[,2,2]) , sd(error_mat[,2,2]) )
resu_table[6,] <- c( mean(dif_error[,1,2]) , sd(dif_error[,1,2]) , mean(dif_error[,2,2]) , sd(dif_error[,2,2]) )

print(round(resu_table,5))
