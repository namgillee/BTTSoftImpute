get_element_tt <- function(V, j)
{
  #get element from blockTT in the LRTT tensor V
  #Ignoring SPT in V
  #
  # V         : LRTT = SPT plus blockTT
  #             a list with arguments N, R, J, G, block, SP
  # V$G       : core tensor of blockTT
  # V$SP      : SPT with $i, $v, $nrow, $ncol, $dimnames
  # j         : vector of indices, with length(j) == V$N

  # Skip argument check

  N = V$N
  R = V$R

  if (N == 0) {
    warning("Input tensor is empty")
    return(0)
  } 
  
  if (N != length(j)) {
    warning("Index length is not equal to tensor dimension")
    return(0)
  }
  
  # if (is.matrix(j)) {
  #   # Each row of j may indicate an element,
  #   #but still this is a loop with each row,
  #   #skip.
  # }

  #
  # We consider that j is a vector
  #
  #if (!is.list(j)) {
  #  j = as.list(j)  #for convenience
  #}
  
  #
  # Main loop
  # compute V_blockTT[j,] + V_sparse[j,]
  # (1) V_blockTT[j,]
  #
  z = diag(1,R[1]) 
  for (n in 1:N) {
    cr = V$G[[n]]
    
    d_cr = dim(cr)
    if (length(d_cr)==4) {
      #Assume 4th order core
      cr = cr[,j[n],,]
    } else if (length(d_cr)==3) {
      #Assume 3rd order core
      cr = cr[,j[n],]
    } else {
      warning('In:get_element_tt::A core tensor does not fit.')
    }
    #if (n == V$block) {
    #  #Assume 4th order core
    #  cr = cr[,j[n],,]
    #} else { 
    #  #Assume 3rd order core
    #  cr = cr[,j[n],]
    #}
    dim(cr) <- c(R[n], length(cr)/R[n])
    
    z = z %*% cr
    dim(z) <- c(length(z)/R[n+1], R[n+1])
  }
  z = as.vector(z)
  
  #Ignoring SPT in V
  # # (2) V_sparse[j,]
  # #
  # idset = which( apply(V$SP[[1]], 1, function(x) all(x == j)) )
  # for (idnz in idset) {
  #   #spj1 = V$SP[[1]][idnz,]  # this is equal to j
  #   spv1 = V$SP[[2]][idnz,]
  #   z = z + spv1
  # }

  return(z)

}