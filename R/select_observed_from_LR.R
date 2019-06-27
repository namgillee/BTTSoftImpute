#' select_observed_from_LR
#'
#' Select observed values randomly from data x.
#' Here, 'randomly' means uniform random sampling.
#'
#' @param x A data matrix in LR (Low-rank (SVD)) format with attributes $u, $d, $v.
#'          $u is a matrix, x$v is a blockTT format.
#' @param numobs Number of observed values to select.
#' @param rate If numobs==NULL, then numobs = rate * numels(full(x)).
#' @return The returned value is a SparseMatrix format, with
#'         attributes $i, $j, $v, $nrow, $ncol, $dimnames.
#' @export
select_observed_from_LR <- function(x, numobs = NULL, rate = 1.0)
{
  retmat = NULL
  sizeI = nrow(x$u)
  J     = x$v$J
  N     = x$v$N
  sizeJ = prod(J)
  numels = sizeI * sizeJ

  if (is.null(sizeI)) {
    warning('Perhaps the input x$u is not a matrix.')
    return(retmat)#NULL
  } else if (is.null(sizeJ)) {
    warning('Perhaps the input x$v is not a blockTT.')
    return(retmat)#NULL
  }


  ## Select 'numobs' values out of 'numels' elements in the matrix data x
  if (is.null(numobs)) {
    numobs = floor(numels*rate)
  }
  if (numobs > numels) {
    warning('Number of observations is larger than number of elements.')
    numobs = numels
  }
  idobs = sample(numels, numobs)


  ## Extract the observed values
  idata = NULL
  jdata = NULL
  vdata = NULL
  for (i in 1:numobs) {
    #Compute (i,j)th element of the matrix data x = list(u,d,v), by
    #          x(i,j) = x$u[i,] %*% diag(x$d) %*% t(x$v[j,])
    #Note that the size of matrix data x is [sizeI x J[1] x ... x J[N]]
    #idobs-1 = i-1 + (j1-1)*I + (j2-1)*I*J1 + ... + (jN-1)*I*J1*...*J{N-1}
    id_vector_this = (idobs[i]-1) %/% cumprod(c(1, sizeI, J[-N]))   #quotient; of length N+1
    id_vector_this = id_vector_this %% c(sizeI, J) + 1    #residual+1; of length N+1

    idata = rbind(idata, id_vector_this[1])
    jdata = rbind(jdata, id_vector_this[-1])
    vdata = rbind(vdata, sum(x$u[id_vector_this[1],] * x$d * get_element_tt(x$v, id_vector_this[-1])))
  }


  # Return value
  retmat = list(i=idata, j=jdata, v=vdata, nrow=sizeI, ncol=J, dimnames=NULL)
  attr(retmat, 'idobs') = idobs
  return(retmat)
}
