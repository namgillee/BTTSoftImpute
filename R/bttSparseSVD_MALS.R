#Low rank tensor approximation (SVD) based on block TT format
#
#We use Modified ALS (MALS) to optimize two TT cores(V_i, v_{i+1}) and left-singular-vectors(U), i.e.,
#   U -> V1 -> V2 -> ...
#
#Inputs
#  Afill = list(Ares, X) : SPLR format
#  Ares   : SP format, $i, $j, $v, $nrow, $ncol
#  X      : LR format, $u, $d, $v
#  X$v    : BlockTT
#  lambda : perform soft-thresholding
#
#Output
#  X      : LR format
#
#Caution: Perhaps this function has not been verified sufficiently..
#
#Last modified: 2018.01.02. by Namgil Lee (Kangwon National University)
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

bttSparseSVD_MALS <- function(Afill, Rmax_TT,
                             tol_u, lambda, maxit, verbose,
                             X0 = NULL, final.svd = FALSE)
{

  if (Afill$X$v$N == 1) {
    if (verbose) {
      print('In MALS: Running ALS instead...')
    }
    X = bttSparseSVD_ALS(Afill = Afill, ttrank_max = Rmax_TT,
                         tol_u = tol_u, lambda = lambda, maxit = maxit, verbose = verbose,
                         X0 = X0, final.svd = final.svd)
    return(X)
  }


  eps = 1e-5

  this.call = match.call()

  # Initialize for X = {u,d,v}
  if (!is.null(X0)) {
    u = X0$u
    d = X0$d
    v = X0$v
  }
  else {
    u = Afill$X$u
    d = Afill$X$d
    v = Afill$X$v
  }

  # Iteration  for updating u, v, d via ALS
  N = v$N
  K = v$K
  ratio <- 1
  iter <- 0
  while ((ratio > tol_u) & (iter < maxit)) {
    iter <- iter + 1
    u_old = u
    v_old = v
    d_old = d

    ##### Update u #####
    # u <- Ares %*% v + Z0 %*% v

    u <- sp_times_btt(Afill$Ares, v) +
         Afill$X$u %*% (Afill$X$d  *  bttDot(Afill$X$v, v))
    if (lambda > 0 &&  sum(d) > .Machine$double.eps)
      u = u %*% diag(d/(d + lambda))
    svdu <- svd(u)
    u <- svdu$u
    d <- svdu$d
    # We skip updating because the bk'th core of v will
    # be updated in the next iteration.
    #v = btt_times_mat(v , svdu$v)


    ##### One full-sweep #####
    #### Note that, the case of $bkstarrt == N$ has not been  ####
    #### implemented, which can cause an error at the moment. ####
    bkstart = v$block
    course = c(bkstart:(N-1))
    if (N>2) {
      course = c(course, (N-2):1)
    }
    if (bkstart>2) {
      course = c(course, 2:(bkstart-1))
    }
    for (i in 1:length(course)) {
      n = course[i]

      ## Update (n,n+1)'th TT-core by
      ##   svd( V_{neq}' * Afill' * U );
      ## where
      ##   the size of the matrix   [V_{neq}' * Afill' * U]
      ##   is of [RIJR x K], and
      ##   the matrix   [Afill = {Ares+X}] is the sum
      ##   of sparse part and low-rank part
      if (n >= v$block) {
        vFrame = get_bttframe(v,1) #V_{neq bk, bk+1}
      } else {
        # n < v$block
        vFrame = get_bttframe(v,-1) #V_{neq bk-1, bk}
      }

      VAU_sp = sp_times_btt(Afill$Ares, vFrame) #SPT
      VAU_sp = t(VAU_sp) %*% u
              #VAU_sp = spt_ttimes_mat_0(VAU_sp, u)

      VAU_lr = bttDot(vFrame, Afill$X$v) %*%
               (Afill$X$d * (t(Afill$X$u)%*%u))  #V_{neq}' V0 D0 U0' U ## [RJR x K] matrix

      VAU = VAU_sp + VAU_lr
      if (lambda > 0  &&  sum(d) > .Machine$double.eps)
        VAU = VAU %*% diag(d/(d + lambda))


      ## Factorize n'th TT-core and merge to the next (either end, rl, lr)
      ## cf: VAU <- svdu$u, factorize and update v$G[[n]]
      ##     d <- svdu$d
      if (i < length(course) && n < course[i+1]) {
        # lr sweep
        # Update G[[n]], G[[n+1]]     # by splitting supercore W[[n,n+1]]
        # The core  G[[n+1]] will have the mode 'K'.

        dim(VAU) = c(v$R[n]*v$J[n], v$J[n+1]*v$R[n+2]*K)
        svdu = svd(VAU)
        Rnew = min( myChop2(svdu$d, eps), Rmax_TT )

        VAU = svdu$u[,seq(Rnew)]
        dim(VAU) = c(v$R[n], v$J[n], Rnew)
        v$G[[n]] = VAU

        # Do not need to update (n+1)th TT-core yet.
        # it will be updated in the next iteration.
        #--
        cr2 = svdu$d[seq(Rnew)] * t(svdu$v[,seq(Rnew),drop=FALSE])
        dim(cr2) = c(Rnew * v$J[n+1] * v$R[n+2], K)
        svdv = svd(cr2)
        cr2 = svdv$u  ##size??
        dim(cr2) <- c(Rnew, v$J[n+1], v$R[n+2], K)
        cr2 = aperm(cr2, c(1,2,4,3))
        v$G[[n+1]] = cr2

        d = svdv$d
        u = u %*% svdv$v
        #--

        v$R[n+1] = Rnew
        v$block = n+1

      } else {
        # rl sweep
        # Update G[[n]], G[[n+1]]     # by splitting supercore W[[n,n+1]]
        # This time, G[[n]] will have the mode 'K'.

        VAU = t(VAU)
        dim(VAU) = c(K*v$R[n]*v$J[n], v$J[n+1]*v$R[n+2])
        svdu = svd(VAU)
        Rnew = min( myChop2(svdu$d, eps), Rmax_TT )

        VAU = t(svdu$v[,seq(Rnew),drop=FALSE])
        dim(VAU) = c(Rnew, v$J[n+1], v$R[n+2])
        v$G[[n+1]] = VAU

        #--
        cr2 = svdu$u[,seq(Rnew),drop=FALSE] %*% diag(svdu$d[seq(Rnew)])
        dim(cr2) = c(K, v$R[n] * v$J[n] * Rnew)
        svdv = svd(t(cr2))
        cr2 = svdv$u
        dim(cr2) = c(v$R[n], v$J[n], Rnew, K)
        cr2 = aperm(cr2, c(1,2,4,3))
        v$G[[n]] = cr2

        d = svdv$d
        u = u %*% svdv$v
        #--

        v$R[n+1] = Rnew
        v$block = n

      }
    }

    ratio = lrRatio(u_old, d_old, v_old, u, d, v)

  }
  if (iter == maxit && verbose)
    warning(paste("Convergence not achieved by", maxit,
                  "iterations"))
  if ((lambda > 0) && final.svd) {
    u <- sp_times_btt(Afill$Ares, v) +
      Afill$X$u %*% (Afill$X$d  *  bttDot(Afill$X$v, v))
    svdu <- svd(u)
    u <- svdu$u
    d <- svdu$d
    v = btt_times_mat(v , svdu$v)   ##Be sure svdu$v is square!!
    d = pmax(d - lambda, 0)
  }
  X = list(u = u, d = d, v = v)
  attr(X,"call") <- this.call
  attr(X,"lambda") <- lambda
  return(X)
}
