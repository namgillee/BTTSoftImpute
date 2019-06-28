#' sp_minus_lr_btt
#'
#' This computes a sparse (SP) matrix minus
#' low-rank (LR) matrix operation.
#'
#' Last modified: 2017.08.25. by Namgil Lee
#'
#' @param A SParse (SP) format with elements i, j, v, nrow, ncol
#' @param Z LR format with elements u, d, and v.
#'        Z$v is in BlockTT format.
#' @export
sp_minus_lr_btt <- function(A, Z)
{
  irow = A$i
  jcol = A$j
  nnz = length(irow)
  if ( nnz == 0 )
    return(A)  #There is no observed element; nothing to subtract
  if ( any(Z$d < 0) ) {
    warning("There is a negative value in singular values")
    return(A)
  }
  if ( sum(Z$d) < .Machine$double.eps )
    return(A)   # because A-Z = A-0 = A

  #Change to matrix for convenience
  if ( !is.array(jcol) ||
       (length(jcol) == nrow(jcol)) ) {
    dim(jcol) = c(length(jcol), 1)
  }

  #Subtract Z from A
  if (!all(Z$v$G[[Z$v$block]] == 0)) {  # simple checking for zero Z$v
    for ( idnz in 1:nnz ) { # We know nnz>0
      m = irow[idnz]   #row index
      j = jcol[idnz, ] #column indices
      vz = sum(Z$u[m,] * Z$d *
               get_element_tt(Z$v, j))
      A$v[idnz] = A$v[idnz] - vz
    }
  }

  return(A)

}
