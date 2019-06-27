#' BTTSoftImpute
#'
#' BTTSoftImpute implements block tensor train decomposition for
#' missing data estimation. It extends the SoftImpute method
#' for matrix data imputation to higher-order tensor data.
#'
#' Last modified: 2018.01.15. by Namgil Lee (Kangwon National University)
#'
#' @param Y   A sparse matrix format with attributes
#'  $i, $j, $v, $nrow, $ncol, and $dimnames,
#'  which are [nnz]-vector, [nnz x N]-matrix, [nnz]-vector, scalar,
#'  [N]-vector, and list, respectively.
#' @param rank_Y   The rank of the estimated low rank matrix, X.
#' @param ttrank_max   The maximal TT rank of right singular vectors, X$v.
#' @return X = list(u, d, v), v is in block TT format
#' @export
#
# < BTT for Missing Data Estimation >
# Copyright (C) 2018 Namgil Lee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

bttdSoftImpute <- function(Y, lambda = 0, rank_Y = 2, ttrank_max = 100,
                           tol_dx = 1e-5, tol_df = 1e-5, tol_f = 0, maxit = 100, verbose = FALSE,
                           X0 = NULL, method = 'ALS')
{

  this.call = match.call()

  #--------------------------------------------------
  # Input parameters
  #
  irow = Y$i        #row indices of nonzero values
  jcol = Y$j        #column indices of nonzero values
  In = Y$nrow       #sizes of rows
  Jn = Y$ncol       #sizes of columns
  N = length(Jn)

  #
  # Set rank_Y, ttrank_max
  #
  # if (rank_Y > (rmax <- min(prod(In),prod(Jn)) - 1)) {
  #   rank_Y  =  rmax
  #   print(paste('In bttdSoftImpute(): rank_Y should not exceed min(dim(Y))-1; Changed to',
  #               rmax))
  # }
  if (rank_Y > (rmax <- min(prod(In),prod(Jn)))) {
    rank_Y  =  rmax
    print(paste('In bttdSoftImpute(): rank_Y should not exceed min(dim(Y)); Changed to',
                rmax))
  }
  ttrank_max  = min(ttrank_max, rank_Y*prod(Jn))
  #-------------------------------------------------


  #-------------------------------------------------
  # Initialize : X = list(u, d, v),
  #                   where v is in block TT format
  #
  if (!is.null(X0)) {

    if (!all(match(c("u", "d", "v"), names(X0), 0) > 0))
      stop("X0 does not have components u, d and v")

    X = X0

    if (!identical(In, nrow(X$u)) || !identical(Jn, X$v$J))
      stop("Sizes of X0 does not match the sizes of input matrix A")


    #====
    # The matrix rank, RA0, of the matrix SVD, X,
    # should be always fixed, whereas the
    # TT-Ranks of X$v can be changed.
    #
    RA0 <- length(X$d)
    if ( RA0 > rank_Y ) {
      # In the case that the rank of X is larger than
      # the given rank value, truncate it.

      X$u = X$u[, seq(rank_Y), drop = FALSE]
      X$d = X$d[seq(rank_Y)]
      bk = X$v$block
      X$v$G[[bk]] <- X$v$G[[bk]][,,seq(rank_Y),, drop = FALSE]

    } else if ( RA0 < rank_Y ) {
      # In the case that the rank of X is smaller than
      # the given rank value,

      ra = rank_Y - RA0    ##Fill in to the rank by size ra.

      # Fill in X$u
      U = X$u
      Ua = matrix(rnorm(In * ra), In, ra)
      Ua = Ua - U %*% (t(U) %*% Ua)
      Ua = svd(Ua)$u
      X$u = cbind(U, Ua)

      # Fill in X$d
      X$d = c(X$d[seq(RA0)], rep(0,ra))

      # Fill in X$v
      bk = X$v$block
      ##Jn = X$v$J
      R = X$v$R
      K = X$v$K #RA0
      cr = X$v$G[[bk]]
      cr = aperm(cr, c(1,2,4,3))
      dim(cr) = c(R[bk]*Jn[bk]*R[bk+1], K)
      Va = matrix(rnorm(R[bk]*Jn[bk]*R[bk+1]*ra),
                  R[bk]*Jn[bk]*R[bk+1], ra)
      Va = Va - cr %*% (t(cr) %*% Va)   #####Assume that bk'th core tensor is orthogonalized
      Va = svd(Va)$u
      cr = cbind(cr, Va)  #<-- ??case rank_Y > nrow??
      dim(cr) = c(R[bk], Jn[bk], R[bk+1], ra)
      X$v$G[[bk]] = aperm(cr, c(1,2,4,3))
      X$v$K = rank_Y
    }
    #====
  } else {
    # Initialize by block TT.
    # We could also use CP for initialization
    # CP init and BlockTT init can be compared
    # next time

    R = c(rank_Y, ceiling(rank_Y/cumprod(Jn[-N])), 1)
    X = list(u = matrix(0,In,rank_Y),  #!!!    #svd(matrix( rnorm(In*rank_Y), In, rank_Y ))$u,
             d = rep(1,rank_Y),
             v = bttRand(Jn, N, R, direction = -1) )

  }
  #-------------------------------------------------


  #-------------------------------------------------
  # Iteration
  # Perform iterations inside of each methods
  X = switch (toupper(method),
    ALS = bttdSsvdImpute_ALS(Y=Y, lambda=lambda, ttrank_max=ttrank_max,
                             tol_dx=tol_dx, tol_df=tol_df, tol_f=tol_f,
                             maxit=maxit, verbose=verbose, X0=X),
    ALS_V01 = bttdSsvdImpute_ALS_v01(Y=Y, lambda=lambda, ttrank_max=ttrank_max,
                                     tol_dx=tol_dx, tol_df=tol_df, tol_f=tol_f,
                                     maxit=maxit, verbose=verbose, X0=X), ##tol_df not yet implemnetd
    SVD = bttdSsvdImpute_SVD(Y=Y, lambda=lambda, ttrank_max=ttrank_max,
                             tol_dx=tol_dx, tol_df=tol_df, tol_f=tol_f,
                             maxit=maxit, verbose=verbose, X0=X),
    MALS = bttdSsvdImpute_MALS(Y=Y, lambda=lambda, ttrank_max=ttrank_max,
                               tol_dx=tol_dx, tol_df=tol_df, tol_f=tol_f,
                               maxit=maxit, verbose=verbose, X0=X),
    {
      print(paste('Unknown method:', method))
    }
  )

  #-------------------------------------------------
  # # OLD -- Iteration
  # Ares = sp_minus_lr_btt(Y, X)
  # obj_trace = 0.5 * sum(Ares$v^2) + lambda * sum(X$d)
  # ratio <- 1
  # iter <- 0
  # while ((ratio > tol_dx) & (iter < maxit)) {
  #   iter <- iter + 1
  #   X_old = X
  #
  #   # ALS routine for computing
  #   # SVD of Afilled == { Ares + X }
  #   Afill = list(Ares=Ares, X=X)
  #
  #   if (toupper(method) == 'ALS') {
  #     X = bttSparseSVD_ALS(Afill, ttrank_max, tol_dx, lambda, maxit,
  #                        verbose = FALSE, X0 = X)
  #
  #   } else if (toupper(method) == 'SVD') {
  #     X = bttSparseSVD_SVD(Afill, ttrank_max, tol_dx, lambda, maxit,
  #                           verbose = FALSE, X0 = X)
  #
  #   } else if (toupper(method) == 'MALS') {
  #     X = bttSparseSVD_MALS(Afill, ttrank_max, tol_dx, lambda, maxit,
  #                           verbose = FALSE, X0 = X)
  #
  #   } else {
  #     warning('Unknown method')
  #   }
  #
  #   # Ares = A_omega - X_omega
  #   Ares = sp_minus_lr_btt(Y, X)
  #
  #   # Convergence trace and stopping
  #   #
  #   ratio = lrRatio(X_old$u, X_old$d, X_old$v, X$u, X$d, X$v)
  #
  #   # Convergence trace
  #   obj = 0.5 * sum(Ares$v^2) + lambda * sum(X$d)
  #   obj_trace = c(obj_trace, obj)
  #   if (verbose) {
  #
  #     cat(iter, ":", "obj", format(round(obj, 5)), "ratio",
  #         ratio, "\n")
  #   }
  #
  # }
  #---------------------------------------------

  #---------------------------------------------
  # Output results
  attr(X, "call") = this.call
  return(X)

}
