#' @export
sp_times_bttframe <- function(spA, V, dir_rem = 0)
{
  # Computes the matrix product  A %*% V_{\neq}
  # Let n=V$block.
  # V_{\neq} is a frame matrix, which is either
  #
  # 1) dir_rem = 0               : V_{\neq n} of size  J{1}...J{N} x R{n}J{n}R{n+1}
  #
  # 2) dir_rem = 1               : V_{\neq n,n+1} of size  J{1}...J{N} x R{n}J{n}J{n+1}R{n+2}
  #
  # 3) dir_rem = -1              : V_{\neq n-1,n} of size  J{1}...J{N} x R{n-1}J{n-1}J{n}R{n+1}
  #
  # spA : sparse format with entries i, j, v, nrow, ncol, dimnames
  # V   : BlockTT with R, N, J, K, G, block
  #       It represents a [prod(V$J) x K] matrix

  irow = spA$i
  jcol = spA$j
  nnz = length(spA$v)
  N = V$N             # order (dimension) of block TT tensor V
  J = V$J
  R = V$R
  bk= V$block
  # if (bk==0 || bk>N)
  #   warning("V$block should be between 1 to N")

  # Set jcol as a matrix
  if ( !is.array(jcol) || (length(jcol) == nrow(jcol)) )
    dim(jcol) = c(length(jcol), 1)

  # Incorrect dimension
  if ( dim(jcol)[2] != N ) {
    warning("Dimension of TT is not equal to that of A")
    return(spA)
  }


  # TT-cores:  1...lend < {lend+1 ... rend-1} < rend...N

  lend = min(bk-1, bk-1 + dir_rem) #products of 1...lend TT cores
  rend = max(bk+1, bk+1 + dir_rem) #products of rend...N TT cores

  lend = max(lend,0)
  rend = min(rend,N+1)


  # Compute U = A %*% V_{\neq}
  U = array(0, c(spA$nrow, R[lend+1], prod(J[(lend+1):(rend-1)]), R[rend]))

  if (spA$nrow > 0 && nnz > 0 &&
      !all(V$G[[N]] == 0)) {    # A basic check for a zero tensor
        for (idnz in 1:nnz) {
          k = irow[idnz]      #row index
          j = jcol[idnz, ]    #column indices
          x = spA$v[idnz]
          ####vrow = get_element_tt(V, j)

          #Multiply from the left
          W = matrix(1, nrow = R[1], ncol = R[N+1])
          if (lend >= 1)
            for (m in 1:lend) {
              cr = V$G[[m]][,j[m],]
              dim(cr) = c(R[m],R[m+1])
              W = t(cr) %*% W
            }

          #Multiply from the right
          if (rend <= N)
            for (m in N:rend) {
              cr = V$G[[m]][,j[m],]
              dim(cr) = c(R[m],R[m+1])
              W = W %*% t(cr)
            }

          #fixed index, jn = j_n + (j_{n+1} -1)*J_n + ... + (j_r -1)*J_n*...*J_{r-1}
          jn = 1 + sum( (j[(lend+1):(rend-1)] - 1) *
                        cumprod( c(1,J[seq(lend+1,rend-1-1,length.out = rend-2-lend)]) ) )
          U[k,,jn,] = U[k,,jn,] + x * W

        }

  }

  dim(U) <- c(spA$nrow, prod(R[lend+1], J[(lend+1):(rend-1)], R[rend]))

  return(U)

}
