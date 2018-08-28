reduce_to_sparse_matrix <- function(X)
{
  # Convert a sparse tensor X into the sparse matrix format
  #
  # Input X has attributes $i, $j, $v, $nrow, $ncol, $dimnames
  #
  # Output X also has the same attributes. 
  # But $i and $j are forced to be in vector instead of matrix, 
  # and $nrow and $ncol are forced to be scalars instead of vectors.
  
  M = length(X$nrow)
  N = length(X$ncol)
  
  if (M>1) {
    idata = X$i  #Consider that idata is a matrix
    cum_prod_rowsizes = cumprod(c(1,X$nrow[-M]))  #(1, J{1}, J{1}J{2}, ..., J{1}***J{N-1})
    idata = colSums( t(idata-1) *  cum_prod_rowsizes ) + 1 # SUM{ (j{n}-1)*J{1}J{2}*...*J{n-1} } + 1
    X$i = idata     #matrix(idata, nrow = length(idata))
    
    X$nrow = prod(X$nrow)
  } else if (is.array(X$i)) {
    X$i = as.vector(X$i)
  }
  
  if (N>1) {
    jdata = X$j  #Consider that jdata is a matrix
    cum_prod_colsizes = cumprod(c(1,X$ncol[-N]))  #(1, J{1}, J{1}J{2}, ..., J{1}***J{N-1})
    jdata = colSums( t(jdata-1) *  cum_prod_colsizes ) + 1 # SUM{ (j{n}-1)*J{1}J{2}*...*J{n-1} } + 1
    X$j = matrix(jdata, nrow = length(jdata))
    
    X$ncol = prod(X$ncol)
  } else if (is.array(X$j)) {
    X$j = as.vector(X$j)
  }
  
  #
  # X$dimnames  :  how to use it is not specified yet.
  #
  
  return(X)

}