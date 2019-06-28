## main-draft-dense_matrices.R
# Simulation for BTTSoftImpute with low-rank dense matrices
#
# Goal: Compare performances of the two methods:
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
### Note: Running time can be very high,
###       it highly depends on, esp., the rate of observed values
###       (parameter name: OBSRATE), so, change simulation parameters
###       including the parameters in the 'for' loops before running.
#
# Last modified at 2018.01.24. by Namgil Lee (Kangwon National University)
#
# < BTT for Missing Data Estimation >
# Copyright (C) 2018  Namgil Lee


## Read source code
rm(list=ls())
set.seed(11111)

#file.sources = list.files("../softImpute_BlockTT/", pattern="*.R", full.names = TRUE)
#sapply(file.sources,source,.GlobalEnv)



## Set Data Sizes
#
# 100 x J^N
#
ORDERS <- seq(1,3,by = 1)     # Not so large orderN, because full format is computed.
sizeI  <- 100                 # Total number of elements in full = I*J^N <= 100*4^6 = 409600.
SIZEJS = c(seq(4,20,by=4))
#sizeJ  <- 20     # J=4:  Considering real image data quantization.
rankRY <- 4     # RY=5: Number of singular values to estimate is not so small.
ttrankY <- 2    # R_n=4: TT-rank of V_Y is not large.

rankRX <- rankRY
myttrankmax <- 10     # R_n(Y) <= R_n(X) : TT-rank of V_X
mytol_df <- -1        # (unused)
mytol_dx <- 5e-4      # Small tol_dx  ---- sufficient accuracy;; 1e-3 / sqrt(5-1) = 5e-4
mytol_f  <- 0         # lambda=0  ==>  obj -> 0
maxiter <- 100000     # Sufficiently large, because we expect tol_df works...

#rateObs <- 3   #num_obs = orderN * rankR^2 * sizeJ  * rateObs
#NUMOBS <- rep(1000,length(ORDERS))
OBSRATE <- seq(0.2,0.9,by=0.1)

noisestd <- 0   #noise standard deviation
lambda = 0
num_rep = 5


time_btt = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS)))
time_mat = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS)))
cost_btt = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS)))
cost_mat = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS)))

mae_btt = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS))) #mean absolute error
mae_mat = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS)))
rmae_btt = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS))) #relative mean absolute error
rmae_mat = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS)))
rmae2_btt = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS))) #relative mean absolute error
rmae2_mat = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS)))

rmse_btt = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS))) #root mean squared error
rmse_mat = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS)))
rrmse_btt = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS))) #relative root mean squared error
rrmse_mat = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS)))
rrmse2_btt = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS))) #relative root mean squared error
rrmse2_mat = array(-1, c(num_rep, length(SIZEJS), length(OBSRATE), length(ORDERS)))



## Main iteration

for (ido in 1:length(ORDERS)) {
  orderN <- ORDERS[ido]
  #print(paste('----- ido:', ido, '-----'))
  for (idrate in 1:length(OBSRATE)) {
    rate_obs = OBSRATE[idrate]
    #print(paste('ido',ido,',rateobs',rate_obs))

    for (idj in 1:length(SIZEJS)) {
      sizeJ = SIZEJS[idj]
      print(paste('ido',ido,',sizeJ',sizeJ,',rateobs',rate_obs))

    for (idrep in 1:num_rep) {
      print(paste('idrep:', idrep))


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
      if (sizeI * sizeJ^orderN > 1e7) {
      } else {
        Y0full = Y0$u %*% diag(Y0$d) %*% t(bttFull(Y0$v))
        svdY0full = svd(Y0full)
        ## Check if 'bttFull()' is working correctly
        # print("singular values of Y0full: ")
        # print(svdY0full$d)
        # print("singular values of Y0 (original): ")
        # print(Y0$d)
      }

      ## Select observed values randomly from the original matrix Y0.
      # The number of observed values: (1-missP)*I*prod(J)
      # Y0_omega is in SparseMatrix format
      Y0_omega = select_observed_from_LR(x = Y0, rate = rate_obs)

      ##To check... if Y0_omega$v == Y0full(idobs)
      ##(<-- that is, selectObservedFromLR() is working well or not)
      # plot(sort(Y0_omega$v), main = 'Selected Values')
      # lines(sort(Y0full[attr(Y0_omega,'idobs')]), col='red')
      # legend('topright', legend = c('from LR', 'from Full'),
      #        lty = c(0,1), pch = c(1,-1))


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
      time_btt[idrep,idj,idrate,ido] = difftime(end_time, start_time, units = 'sec')

      Yres = sp_minus_lr_btt(Y_omega, Z_btt_als)
      obj = 0.5 * sum(Yres$v^2) + lambda * sum(Z_btt_als$d)
      cost_btt[idrep,idj,idrate,ido] = obj


      #****
      if (sizeI * sizeJ^orderN > 1e7) {
        print("****Too Large Size To Compute Full Matrix/Tensor*****")
      } else {
        Y0full = Y0$u %*% diag(Y0$d) %*% t(bttFull(Y0$v))  #Full. Computed once
        Zbres = Z_btt_als$u %*% diag(Z_btt_als$d) %*% t(bttFull(Z_btt_als$v))
        Zbres = Zbres - Y0full  #Residual
        mae_btt[idrep,idj,idrate,ido] = mean(abs(Zbres))
        rmae_btt[idrep,idj,idrate,ido] = mean(abs(Zbres))/mean(abs(Y0full))  # <- TO PAPER
        rmae2_btt[idrep,idj,idrate,ido] = mean(abs(Zbres)/abs(Y0full))
        rmse_btt[idrep,idj,idrate,ido] = sqrt(mean(Zbres^2))
        rrmse_btt[idrep,idj,idrate,ido] = sqrt(mean(Zbres^2))/sqrt(mean(Y0full^2))
        rrmse2_btt[idrep,idj,idrate,ido] = sqrt(mean((Zbres)^2/(Y0full)^2))
      }
      #****



      ## 2) Using standard thresholded SVD
      ##   We use our ALS algorithm but using matrix data, i.e.,
      #    Transform Y$j from matrix (array) into vector (scalars)
      # # require(softImpute)
      # # Ymat_omega1 = reduce_to_sparse_matrix(Y_omega)
      # # Ymat_omega = Incomplete(as.vector(Ymat_omega1$i),
      # #                         as.vector(Ymat_omega1$j),
      # #                         as.vector(Ymat_omega1$v))
      if (rankRY * sizeJ^orderN <= 2e7) {
        Ymat_omega = reduce_to_sparse_matrix(Y_omega)

        start_time = Sys.time()
        Z_mat_als = bttdSoftImpute(Y = Ymat_omega, lambda = 0,
                                   rank_Y = rankRX, ttrank_max = myttrankmax,
                                   tol_dx = mytol_dx, tol_df = mytol_df, tol_f = mytol_f,
                                   verbose = TRUE, method = 'als', maxit = maxiter)
        end_time = Sys.time()
        time_mat[idrep,idj,idrate,ido] = difftime(end_time, start_time, units = 'sec')

        Ymatres = sp_minus_lr_btt(Ymat_omega, Z_mat_als)
        obj = 0.5 * sum(Ymatres$v^2) + lambda * sum(Z_mat_als$d)
        cost_mat[idrep,idj,idrate,ido] = obj
      }



      #****
      #if (sizeI * sizeJ^orderN > 1e7)
      #  print("Too Large Size To Compute Full Matrix/Tensor")
      #Y0full = Y0$u %*% diag(Y0$d) %*% t(bttFull(Y0$v))  #Full. Computed once
      if (sizeI * sizeJ^orderN > 1e7) {
      } else {
        Zbres = Z_mat_als$u %*% diag(Z_mat_als$d) %*% t(bttFull(Z_mat_als$v))
        Zbres = Zbres - Y0full  #Residual
        mae_mat[idrep,idj,idrate,ido] = mean(abs(Zbres))
        rmae_mat[idrep,idj,idrate,ido] = mean(abs(Zbres))/mean(abs(Y0full))
        rmae2_mat[idrep,idj,idrate,ido] = mean(abs(Zbres)/abs(Y0full))
        rmse_mat[idrep,idj,idrate,ido] = sqrt(mean(Zbres^2))
        rrmse_mat[idrep,idj,idrate,ido] = sqrt(mean(Zbres^2))/sqrt(mean(Y0full^2))
        rrmse2_mat[idrep,idj,idrate,ido] = sqrt(mean((Zbres)^2/(Y0full)^2))
      }
      #****



    }#end for num_rep

  }#end for idrate
  }#end for sizeJ

  ## Save the results
  save(file = paste('resu_main-dense_matrices-forTable.Rdata' , sep=''),
       time_btt, cost_btt, time_mat, cost_mat,
       mae_btt,mae_mat,
       rmae_btt,rmae_mat,
       rmae2_btt,rmae2_mat,
       rmse_btt,rmse_mat,
       rrmse_btt,rrmse_mat,
       rrmse2_btt,rrmse2_mat
  )


}#end for ORDERS


## Print summary statistics : result tables has dimensions of [idrep,idj,idrate,ido]
## Manually change the indices of (idrate,ido) and print the results.
idrate=6; ido=2; print( summary(rrmse_btt[,,idrate,ido] - rrmse_mat[,,idrate,ido]) )


## Draw boxplots
idrate=6; ido=2;
if (idrate>=3 && ido==3) {
  ## For high-order, tensor is large, so,
  ## select a subset of three sizes; Jn=(4,8,12)
  dif_rrmse2 = dif_rrmse[,1:3,,]
  SIZEJS2 = SIZEJS[1:3]

} else {
  ## For lower-order, remove
  ## one of the box for simplicity
  dif_rrmse2 = dif_rrmse[,2:5,,]
  SIZEJS2 = SIZEJS[2:5]
}
mybox = boxplot(dif_rrmse2[,,idrate,ido])
#pdf(paste0(outdir,'/boxplot_Drrmse_rate', 10*OBSRATE[idrate],'_N',ORDERS[ido],'.pdf'))
#{
ymax = max(abs(c(mybox$stats[1,], mybox$stats[5,])))
boxplot(dif_rrmse2[,,idrate,ido], ylim = c(-ymax,ymax),
        names = SIZEJS2, cex.axis=2,
        xlab = expression('J'[n]), cex.lab=2,
        ylab = expression(paste(Delta,'RRMSE')))
abline(h=0,col='red',lty=3,lwd=2)
#}
#dev.off()

