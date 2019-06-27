#' btt_times_mat
#'
#' This compute BlockTT times matrix,
#' and return BlockTT format
btt_times_mat <- function(tt, U)
{
  bk = tt$block
  K = tt$K

  cr = tt$G[[bk]]
  dim(cr) = c(tt$R[bk] * tt$J[bk] , K, tt$R[bk+1])
  cr = aperm(cr, c(1,3,2))
  dim(cr) = c(tt$R[bk] * tt$J[bk] * tt$R[bk+1], K)
  cr = cr %*% U
  dim(cr) = c(tt$R[bk] , tt$J[bk] , tt$R[bk+1], dim(U)[2])
  tt$G[[bk]] = aperm(cr, c(1,2,4,3))


  tt$K = dim(U)[2]
  return(tt)
}
