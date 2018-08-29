## main-draft-dense_matrices-forConvergence.R
# Simulation for BTTSoftImpute with low-rank dense matrices
# 
# Goal: Compare convergence of the two methods: 
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
num_rep = 1



time_btt = NULL
conv_btt = NULL
time_mat = NULL
conv_mat = NULL
resu_all = NULL

## Main iteration
for (ido in 1:length(ORDERS)) 
{
  orderN <- ORDERS[ido]
  print(paste('----- Order:', orderN, '-----'))
      
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
  
  
  ## Calculate full noiseless matrix from Y0
  ## Ony when the numel of Y0 is sufficiently small
  ## And check if SVD of Y0 produces back the (U0, d0, V0).
  # if (sizeI * sizeJ^orderN > 1e7) { 
  # } else {
  #   Y0full = Y0$u %*% diag(Y0$d) %*% t(bttFull(Y0$v))
  #   svdY0full = svd(Y0full)
  # }
  
  ## Select observed values randomly from the original matrix Y0. 
  # The number of observed values: (1-missP)*I*prod(J)
  # Y0_omega is in SparseMatrix format
  Y0_omega = select_observed_from_LR(x = Y0, rate = rate_obs)
  
  
  
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
  time_btt[ido] = difftime(end_time, start_time, units = 'sec')
  conv_btt[[ido]] = attr(Z_btt_als,'obj')

  
  ## 2) Using standard thresholded SVD
  Ymat_omega = reduce_to_sparse_matrix(Y_omega)
  
  start_time = Sys.time()  
  Z_mat_als = bttdSoftImpute(Y = Ymat_omega, lambda = 0, 
                             rank_Y = rankRX, ttrank_max = myttrankmax, 
                             tol_dx = mytol_dx, tol_df = mytol_df, tol_f = mytol_f,
                             verbose = TRUE, method = 'als', maxit = maxiter)
  end_time = Sys.time()  
  time_mat[ido] = difftime(end_time, start_time, units = 'sec')
  conv_mat[[ido]] = attr(Z_mat_als,'obj')
  

}#end for ORDERS



## Save the results
save(file = paste('resu_main-draft-dense_matrices-forConvergence.Rdata' , sep=''), 
     time_btt, conv_btt, time_mat, conv_mat,
     ORDERS
)


## Plot convergence
for (idx in 1:length(ORDERS)) { 
  
  ylim_max = 0.7 #Fix the ylim_max;  #max(conv_btt[[idx]])
  ylim_min = 1e-8
  xlim_max = 60  #Fix the xlim_max; 
  
  pdf(paste(outdir, '/convergence_dense_sizeJ',sizeJ,'_N',ORDERS[idx],'.pdf', sep=''))
  plot(c(0,xlim_max), c(ylim_min, ylim_max), type='n', log='y',
       xlab = 'Iteration', ylab = 'Cost function', main = '', 
       cex.lab=1.5, cex.main = 1.5)
  lines(seq(0,length(conv_btt[[idx]])-1), 
        conv_btt[[idx]], type = 'l', col=1)
  lines(seq(0,length(conv_mat[[idx]])-1), 
        conv_mat[[idx]], type = 'l', col=2)
  #----Plot for over-layered points----#
  lines(seq(0,length(conv_btt[[idx]])-1,by=10), 
        conv_btt[[idx]][seq(1,length(conv_btt[[idx]]),by=10)], 
          type = 'p', col=1, pch=1, cex=2)
  lines(seq(0,length(conv_mat[[idx]])-1,by=10), 
        conv_mat[[idx]][seq(1,length(conv_mat[[idx]]),by=10)], 
        type = 'p', col=2, pch=4, cex=2)
  legend('topright', legend = c('BTT','Matrix'), 
         lty=1, pch=c(1,4), col=1:2, cex=1.5)
  dev.off()
  
  png(paste(outdir, '/convergence_dense_sizeJ',sizeJ,'_N',ORDERS[idx],'.png', sep=''))
  plot(c(0,xlim_max), c(ylim_min, ylim_max), type='n', log='y',
       xlab = 'Iteration', ylab = 'Cost function', main = '', 
       cex.lab=1.5, cex.main = 1.5)
  lines(seq(0,length(conv_btt[[idx]])-1), 
        conv_btt[[idx]], type = 'l', col=1)
  lines(seq(0,length(conv_mat[[idx]])-1), 
        conv_mat[[idx]], type = 'l', col=2)
  #----Plot for over-layered points----#
  lines(seq(0,length(conv_btt[[idx]])-1,by=10), 
        conv_btt[[idx]][seq(1,length(conv_btt[[idx]]),by=10)], 
        type = 'p', col=1, pch=1, cex=2)
  lines(seq(0,length(conv_mat[[idx]])-1,by=10), 
        conv_mat[[idx]][seq(1,length(conv_mat[[idx]]),by=10)], 
        type = 'p', col=2, pch=4, cex=2)
  legend('topright', legend = c('BTT','Matrix'), 
         lty=1, pch=c(1,4), col=1:2, cex=1.5)
  dev.off()
}

