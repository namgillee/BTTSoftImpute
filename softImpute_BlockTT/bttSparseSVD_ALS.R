bttSparseSVD_ALS <- function(Afill, ttrank_max, 
                             tol_u, lambda, maxit, verbose, 
                             X0 = NULL, final.svd = FALSE)
{
  #Low rank tensor approximation (SVD) based on block TT format
  #
  #We use ALS to optimize each TT cores(V_i) and left-singular-vectors(U), i.e., 
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
  #
  #Caution: Perhaps this function has not been verified sufficiently..
  #
  #Last modified: 2017.08.25. by Namgil Lee (Kangwon National University)
  #
  # < Regularized Block Tensor Train Decomposition >
  # Copyright (C) 2018 Namgil Lee
  #
  ##################################################################################
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
  ##################################################################################

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
    bkstart = v$block
    course = c(bkstart:N)
    if (N>1) {
      course = c(course, (N-1):1)
    } 
    if (bkstart>1) {
      course = c(course, 2:bkstart)
    }
    for (i in 1:length(course)) {
      n = course[i]
      
      ## Update n'th TT-core by
      ##    svd( V_{neq}' * Afill' * U ); 
      ## Afill={Ares+X} is the sum of sparse part and low-rank part
      vFrame = get_bttframe(v) #V_{neq} in BlockTT
      
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
      if (i==length(course)) {
        # end of the loop
        svdu <- svd(VAU)
        
        VAU <- svdu$u
        dim(VAU) = c(v$R[n], v$J[n], v$R[n+1], K)  ##??Assume K did not change##
        VAU = aperm(VAU, c(1,2,4,3))
        v$G[[n]] = VAU
        
        d <- svdu$d
        
        u = u %*% svdu$v
        
      } else if (n < (n2 <- course[i+1])) { 
        # lr sweep
        dim(VAU) = c(v$R[n]*v$J[n], v$R[n+1]*K) 
        svdu = svd(VAU)
        Rnew = min( myChop2(svdu$d, eps), ttrank_max )
        
        VAU = svdu$u[,seq(Rnew)]
        dim(VAU) = c(v$R[n], v$J[n], Rnew)
        v$G[[n]] = VAU
          
        # Do not need to update d yet
        # AND
        # We skip updating the next core because it will 
        # be updated in the next iteration from the other cores
        #--
        # cr2 = v$G[[n2]]
        # sz2 = dim(cr2)
        # dim(cr2) = c(sz2[1], prod(sz2[-1]))
        # vs = svdu$v[,seq(Rnew)] %*% diag(svdu$d[seq(Rnew)])
        # dim(vs) = c(v$R[n+1], K*Rnew )
        # cr2 = t(vs) %*% cr2    #sz2[1]==R[n2]==R[n+1]
        # dim(cr2) = c(K, Rnew, sz2[-1])  #K, Rnew, J[n2], R[n2+1]
        # cr2 = aperm(cr2, c(2,3,1,4))
        # v$G[[n2]] = cr2
        
        v$R[n+1] = Rnew
        v$block = n2
        
      } else {
        # rl sweep
        n2 <- course[i+1]  # We know n>course[i+1]
        
        VAU = t(VAU)
        dim(VAU) = c(K*v$R[n], v$J[n]*v$R[n+1]) 
        svdu = svd(VAU)
        Rnew = min( myChop2(svdu$d, eps), ttrank_max )
        
        VAU = svdu$v[,seq(Rnew)]
        VAU = t(VAU)
        dim(VAU) = c(Rnew, v$J[n], v$R[n+1])
        v$G[[n]] = VAU
        
        # Do not need to update d yet
        # AND
        # We skip updating the next core because it will 
        # be updated in the next iteration from the other cores
        #--
        # cr2 = v$G[[n2]]
        # sz2 = dim(cr2)
        # dim(cr2) = c(prod(sz2[-length(sz2)]), sz2[length(sz2)])
        # vs = svdu$u[,seq(Rnew)] %*% diag(svdu$d[seq(Rnew)])
        # dim(vs) = c(K, v$R[n], Rnew )
        # vs = aperm(vs, c(2,1,3))
        # dim(vs) = c(v$R[n], K*Rnew)
        # cr2 = cr2 %*% vs    #sz2[length(sz2)]==R[n]
        # dim(cr2) = c(sz2[-length(sz2)], K, Rnew)  #R[n2], J[n2], K, Rnew
        # v$G[[n2]] = cr2
        
        v$R[n] = Rnew
        v$block = n2
        
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