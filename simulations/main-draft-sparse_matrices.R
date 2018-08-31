## main-draft-sparse_matrices.R
# Simulation for BTTSoftImpute with sparsely observed random Matrices
#
# Check how the observed values are estimated by our method
#
# Last modified at 2018.08.28. by Namgil Lee (Kangwon National University)
#
# < BTT for Missing Data Estimation >
# Copyright (C) 2018  Namgil Lee

rm(list=ls())
set.seed(33333)

file.sources = list.files("../softImpute_BlockTT/", pattern="*.R", full.names = TRUE)
sapply(file.sources,source,.GlobalEnv)

outdir = 'sparse_matrices'
if (!file.exists(outdir)) {
  dir.create(outdir)
}

orderN <- 2
sizeJ <- 20
RANKRS <- c(1:5) 
myttrankmax <- 10
mytol_x <- 1e-6
maxiter <- 10000


convergence_all = NULL
data_values_all = NULL
data_estims_all = NULL
time_cost_all = NULL



## Create a random matrix whose nonzero values are from a uniform distribution
data_values = runif(sizeJ, min = 0, max = 1)
Y_omega = list(i = 1:sizeJ, 
               j = matrix(data = rep(1:sizeJ, times = orderN), nrow = sizeJ), 
               v = data_values, 
               nrow = sizeJ, 
               ncol = rep(sizeJ, orderN), 
               dimnames = NULL)
data_values_all = data_values



for (idr in RANKRS) {
    rankR <- RANKRS[idr]
    
    print(paste('Rank:', rankR))



    ## Conduct tensor completion by using block-TT Decomposition (N. LEE)
    start_time = Sys.time()  
    Z_btt_als = bttdSoftImpute(Y = Y_omega, lambda = 0, 
                               rank_Y = rankR, ttrank_max = myttrankmax, tol_dx = mytol_x, 
                               verbose = FALSE, method = 'als', maxit = maxiter)
    end_time = Sys.time()  
    
    print(end_time - start_time)
    time_cost_all[idr] = end_time - start_time
    convergence_all[[idr]] = attr(Z_btt_als,'obj')




    ## Compare True value versus Estimated Value  
    my_est = Z_btt_als$u %*% diag(Z_btt_als$d, length(Z_btt_als$d)) %*% t(bttFull(Z_btt_als$v))
    id_obs = 1:sizeJ
    id_obs = id_obs + (id_obs-1)*sizeJ + (id_obs-1)*sizeJ^2
    data_observed = my_est[id_obs]
    data_estims_all[[idr]] = data_observed
    
    jpeg(paste(outdir, '/observed_values_R', rankR, '.jpg', sep=''))
    plot(sort(Y_omega$v), ylim = c(0,1), ylab = 'Observed values', xlab = 'Index', 
         cex.lab = 1.5, cex.main = 1.5)
    lines(sort(data_observed), col = 'red')#, pch = 4, type='b')
    legend('topleft', legend=c('Observed', 'Estimated'),
           col=c('black','red'), pch=c(1,4), lty=c(0,1))
    dev.off()
    
    pdf(paste(outdir, '/observed_values_R', rankR, '.pdf', sep='')) #width=10/2.54, height=6/2.54
    plot(sort(Y_omega$v), ylim = c(0,1), ylab = 'Observed values', xlab = 'Index', 
         cex.lab = 1.5, cex.main = 1.5)
    lines(sort(data_observed), col = 'red')#, pch = 4, type='b')
    legend('topleft', legend=c('Observed', 'Estimated'),
           col=c('black','red'), pch=c(1,-1), lty=c(0,1))
    dev.off()
    
    
    
}#end for RANKRS

## Plot convergence 
cost_max = max(convergence_all[[1]])

jpeg(paste(outdir, '/convergence.jpg', sep=''))
plot(c(0,3500), c(0, cost_max), type='n', 
     xlab = 'Iteration', ylab = 'Cost function', 
     main = 'Convergence of BTTSoftImpute method', 
     cex.lab=1.5, cex.main=1.5)
count <- min(5,length(RANKRS))
for (idr in 1:count) {
  lines(seq(0,length(convergence_all[[idr]])-1), 
        convergence_all[[idr]], type = 'l', col=idr, pch=idr)
}
for (idr in 1:count) {
  lines(seq(0,length(convergence_all[[idr]])-1,by=500), 
        convergence_all[[idr]][seq(1,length(convergence_all[[idr]]),by=500)], 
        type = 'p', col=idr, pch=idr)
}
legend('topright', legend = paste('R_X =', RANKRS[1:count]), 
       lty=1, pch=1:length(RANKRS), col=1:length(RANKRS))
dev.off()

pdf(paste(outdir, '/convergence.pdf', sep=''))
    #width=10/2.54, height=6/2.54)
plot(c(0,3500), c(0, cost_max), type='n', 
     xlab = 'Iteration', ylab = 'Cost function', 
     main = 'Convergence of BTTSoftImpute method', 
     cex.lab=1.5, cex.main = 1.5)
count <- min(5,length(RANKRS))
for (idr in 1:count) {
  lines(seq(0,length(convergence_all[[idr]])-1), 
        convergence_all[[idr]], type = 'l', col=idr, pch=idr)
}
for (idr in 1:count) {
  lines(seq(0,length(convergence_all[[idr]])-1,by=500), 
        convergence_all[[idr]][seq(1,length(convergence_all[[idr]]),by=500)], 
        type = 'p', col=idr, pch=idr)
}
legend('topright', legend = paste('R_X =', RANKRS[1:count]), 
       lty=1, pch=1:length(RANKRS), col=1:length(RANKRS))
dev.off()



