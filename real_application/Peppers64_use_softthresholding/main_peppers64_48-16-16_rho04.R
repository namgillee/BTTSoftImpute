## main_peppers64_48-16-16_rho04.R
# Real Application of BTTSoftImpute to Color Images: Peppers in resized 64x64x3
#
# Last modified at 2018.08.30. by Namgil Lee (Kangwon National University)
#
# < BTT for Missing Data >
# Copyright (C) 2018  Namgil Lee


library(tiff)
library(imager)

rm(list=ls())
file.sources = list.files("../../softImpute_BlockTT/", pattern="*.R", full.names = TRUE)
sapply(file.sources,source,.GlobalEnv)

set.seed(123)


datdir = '../data/'
filename = paste0(datdir, 'peppers.tiff')
outdir = './'
prefix = 'main03_peppers64_'


## Read Image

# grayscale image, 512 x 512 x 3
dat = readTIFF(filename) 

# Resize image
dim(dat) <- c(512, 512, 3, 1)
dat <- resize(dat, 2^6, 2^6) #256 x 256
dim(dat) <- c(2^6, 2^6, 3)

# # Save resized original image
# writeTIFF(dat, paste0(outdir, prefix, 'peppers64_orig.tiff'))


## Image Reconstruction

# Remove pixels at random 
II <- dim(dat)
num_elem <- prod(II)  #12288
prop_obs <- 0.4   #60% is missing
num_obs <- round(prop_obs * num_elem)  #4915
idx_obs <- sample(num_elem, num_obs)  #index of the observed values
idx_mis <- setdiff(1:num_elem, idx_obs)

# Save the image with only observed values,
# i.e., missing values have been removed and set at 0.
dat_rem <- dat
dat_rem[idx_mis] <- 0
# writeTIFF(dat_rem, paste0(outdir, prefix, 'peppers64_orig_rem.tiff'))

# Transform to a sparse tensor
# Original resized data has size 256x256x3 
# which will be changed to --> (2x2x2x2x2x2x2x2) x (2x2x2x2x2x2x2x2) x 3
SizeTo <- c(rep(2,6),rep(2,6),3)
vals <- dat[idx_obs]
idxs <- matrix(0, num_obs, length(SizeTo)) # (i, j1,j2,...,jN) coordinate

remainder <- idx_obs - 1
Qs <- cumprod(c(1, SizeTo[-length(SizeTo)]))
for (j in length(Qs):1) {
  divi <- Qs[j]
  quo <- floor(remainder/divi)  #quotient
  remainder <- remainder - divi * quo
  
  idxs[,j] <- quo + 1
}

# Permute sparse tensor 
# from size (2x2x2x2x2x2) x (2x2x2x2x2x2) x 3
# to size 3x(2x2)x...x(2x2)
idxs_permuted <- idxs[,c(13,1,7,2,8,3,9,4,10,5,11,6,12)]   ################

# Reshape sparse tensor 
# to size (3*16)x16x16x16
idxs_reshaped <- matrix(0,num_obs,3)

idxs_reshaped[,1] <- idxs_permuted[,1] + (idxs_permuted[,2]-1)*3 + 
  (idxs_permuted[,3]-1)*3*2^1 + 
  (idxs_permuted[,4]-1)*3*2^2 + (idxs_permuted[,5]-1)*3*2^3
idxs_reshaped[,2] <- idxs_permuted[,6] + (idxs_permuted[,7]-1)*2^1 + 
  (idxs_permuted[,8]-1)*2^2 + (idxs_permuted[,9]-1)*2^3
idxs_reshaped[,3] <- idxs_permuted[,10] + (idxs_permuted[,11]-1)*2^1 + 
  (idxs_permuted[,12]-1)*2^2 + (idxs_permuted[,13]-1)*2^3



Y_omega = list(i=idxs_reshaped[,1], 
               j=idxs_reshaped[,-1], 
               v=vals, 
               nrow=3*16, ncol=c(16,16), 
               dimnames=NULL)


## Conduct tensor completion by using block-TT(by Lee) and standard thresholded SVD(by Mazumder&Hastie)
## 1) Using block-TT-ALS
rankRX = 48             #<-- It should have been larger, but computation...
myttrankmax = 40        #<-- It could be larger, but computation... It will incur a regularization
mytol_f  <- 0         # lambda=0  ==>  obj -> 0
mytol_df <- 1e-5        # (unused) 1e-2
mytol_dx <- 1e-5      #

maxiter <- 100     # <---  small for test purpose


for (lambda in c(0.00, 0.01, 0.10, 0.50, 1.00, 5.00)) {  # (lambda in 1:10) {

  start_time = Sys.time()  
  Z_btt_als = bttdSoftImpute(Y = Y_omega, lambda = lambda, 
                             rank_Y = rankRX, ttrank_max = myttrankmax, 
                             tol_dx = mytol_dx, tol_df = mytol_df, tol_f = mytol_f,
                             verbose = TRUE, method = 'als', maxit = maxiter)
  end_time = Sys.time()  
  
  
  TIME_TAKEN = difftime(end_time, start_time, units = 'sec')
  
  Yres = sp_minus_lr_btt(Y_omega, Z_btt_als)
  obj = 0.5 * sum(Yres$v^2) + lambda * sum(Z_btt_als$d)
  COST_FUN = obj
  
  
  ## Reconstruct Images
  
  # X = USV'
  Zbtt = Z_btt_als$u %*% diag(Z_btt_als$d) %*% t(bttFull(Z_btt_als$v)) #48x16^3
  
  # Reshape
  dim(Zbtt) <- c(3,rep(2,4),rep(2,4),rep(2,4)) #17 dimensional
  
  # Permute to original 3D image
  #inverse permuation of : idxs[,c(13,1,7,2,8,3,9,4,10,5,11,6,12)]
  dat_lowrank <- aperm(Zbtt, c(2,4,6,8,10,12, 3,5,7,9,11,13, 1))
  dim(dat_lowrank) <- c(2^6, 2^6, 3)
  dat_lowrank[dat_lowrank>1] <- 1
  dat_lowrank[dat_lowrank<0] <- 0
  
  # # Save reconstructed image 
  #writeTIFF(dat_lowrank, paste0(outdir, prefix, 'recon00_rho', prop_obs, '.tiff'))
  
  #overritde observed value 
  #because Yhat = X + (Y_omega - X_omega)
  dat_lowrank[idx_obs] <- dat[idx_obs]  
  
  # Save reconstructed image 
  filestr = 
  writeTIFF(dat_lowrank, paste0(outdir, prefix, 
                                'recon_rho', 10* prop_obs, '_lam', 100*lambda, '.tiff'))
  
  
  #Yres
  print(paste('lambda:', lambda))
  print(TIME_TAKEN)
  print(COST_FUN)
  print(Z_btt_als$d)
  print(Z_btt_als$v$R)
  print(Z_btt_als$v$J)
}
