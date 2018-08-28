spt_ttimes_mat_0 <- function(A, u)
{
  #transpose of sparse tensor (SP) and multiplying with a matrix 
  #In this version, we assume
  #
  # 
  # A : sparse tensor representing a matrix of size [M x L]
  #     $i (vector), $nrow, $v (matrix of size nnz x L)
  #     In general [SP]$i is a matrix, but 
  #     in this version assume A$i is a vector, or only use its first column 
  # u : a full matrix, size MxK
  # 
  # Output: A'u = an [L x K] matrix 
  
  iA1 = A$i
  if (is.matrix(iA1)) {
    warning('spt_ttimes_mat_0(): Input tensor A having multiple dimensions is not allowed')
  }

    
  # A$v : nnz x L matrix 
  # u[iA1,] : nnz x K matrix
  # out : L x K matrix
  #####out = array(0, c(dim(A$v)[2], dim(u)[2], dim(iA)[2]))
  out = t(A$v) %*% u[iA1,]  # Summation is also carried out  because iA1 is not matrix
  return(out)
  
}