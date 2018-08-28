sp_times_btt <- function(spA, V)
{
  # Computes the matrix product  A %*% V 
  #
  # spA : sparse format with entries i, j, v, nrow, ncol, dimnames
  # V   : BlockTT with R, N, J, K, G, block
  #       It represents a [prod(V$J) x K] matrix
  
  bk = V$block
  
  U = sp_times_bttframe(spA, V, dir_rem = 0)
  
  cr = V$G[[bk]]
  cr = aperm(cr, c(1,2,4,3))
  dim(cr) <- c(V$R[bk] * V$J[bk] * V$R[bk+1], V$K)
  
  U = U %*% cr
  
  #***** Previous Old Code *******************************
  # irow = spA$i
  # jcol = spA$j
  # M = spA$nrow
  # nnz = length(irow)
  # N = V$N             # order (dimension) of block TT tensor V
  # K = V$K
  # # bk= V$block
  # # if (bk==0 || bk>N)
  # #   warning("V$block should be between 1 to N")
  # 
  # # Set jcol as a matrix
  # if ( !is.array(jcol) || (length(jcol) == nrow(jcol)) )
  #   dim(jcol) = c(length(jcol), 1)
  # 
  # # Incorrect dimension
  # if ( dim(jcol)[2] != N ) {
  #   warning("Dimension of TT is not equal to that of A")
  #   return(spA)
  # } 
  # 
  # # Compute U = A %*% V
  # U = matrix(0, M, K)
  # if (M==0 || nnz==0) 
  #   return(U)
  # 
  # 
  # if (!all(V$G[[N]] == 0)) { # A basic check for a zero tensor
  #   for (idnz in 1:nnz) {
  #     m = irow[idnz]      #row index
  #     j = jcol[idnz, ]    #column indices
  #     x = spA$v[idnz]
  #     vrow = get_element_tt(V, j)
  #     U[m,] = U[m,] + x * vrow
  #   }
  # }
  #**************************************************

  return(U)

}