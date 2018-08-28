bttdSsvdImpute_ALS <- function(Y, lambda, ttrank_max, tol_dx=-1, tol_df=-1, tol_f=0,
                               maxit, verbose, X0 = NULL, final.svd=FALSE) 
{
  #Sparse SVD Imputation based on block TT format
  #Use ALS to optimize each TT cores(V_i) and left-singular-vectors(U), i.e., 
  #   U -> V1 -> V2 -> ...
  #
  # *************************************************************************
  # * In the version 'als', the matrix Afill = {Ares + udv'} are updated at *
  # * every iteration of updating u and each core tensor V_n, i.e.,         *
  # * U -> (update) -> V1 -> (update) -> V2 -> ... -> V1 -> (update)        *
  # *************************************************************************
  #
  #For inputs, see bttdSoftImpute.R
  #Assume that X0 is not NULL but given. 
  #Output
  #  X      : LR format, with attribute 'obj'
  #
  #Last modified: 2018.1.15. by Namgil Lee (Kangwon National University)
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
    warning('X0 is not given.')
    return(NULL)
  }
  
  # Iteration  for updating u, v, d via ALS
  N = v$N
  K = v$K
  Ares = sp_minus_lr_btt(Y, list(u=u,d=d,v=v))
  obj_trace = 0.5 * sum(Ares$v^2) + lambda * sum(d)
  ratio <- 1
  iter <- 0
  while ((ratio > tol_dx) & (iter < maxit)) {
    iter <- iter + 1
    u_old = u
    v_old = v
    d_old = d 

    #----
    # Assume that Afill = {Ares + udv'}
    # Note that, in this process, u,d,v are updated 
    # during iterations, instead of being fixed.
    # If they were fixed, it means performing SVD 
    # of Afill strictly.
    #----
    
    # Update u
    # u <- Afill %*% v
    u <- sp_times_btt(Ares, v) + u %*% diag(d,length(d)) 
    if (lambda > 0)   #&&  sum(d) > .Machine$double.eps (d may have been initialized by zeros)
      u = u %*% diag(d/(d + lambda),length(d))
    # Orthogonalize
    svdu <- svd(u)
    u <- svdu$u
    d <- svdu$d
    v = btt_times_mat(v, svdu$v)     #This line is necessary only to compute 'Ares=Y-X' accurately.

    
    # Update v
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

      # Ares = A_omega - X_omega
      Ares = sp_minus_lr_btt(Y, list(u=u,d=d,v=v))
      
      
      # Update v$G{n}
      # v$G{n} <- V_{neq}' * Afill' * u
      #        =  V_{neq}' * Ares' * u + V_{neq}' * vdu' * u
      #        =  V_{neq}' * Ares' * u + V_n * d
      ####vFrame = get_bttframe(v) #V_{neq} in BTT
      ####VAU = sp_times_btt(Ares, vFrame) #V_{neq}' * Ares'
      VAU = sp_times_bttframe(Ares, v) # This line abbreviates: VAU=sp_times_btt(Ares,get_bttframe(v))
      VAU = t(VAU) %*% u    #(V_{neq}' * Ares') * u
      cr = aperm(v$G[[n]], c(1,2,4,3)) 
      dim(cr) = c(v$R[n] * v$J[n] * v$R[n+1], K)
      cr = cr %*% diag(d,length(d))   #V_n * d
      VAU = VAU + cr
      if (lambda > 0)  #  &&  sum(d) > .Machine$double.eps
        VAU = VAU %*% diag(d/(d + lambda),length(d))
      # Orthogonalize
      svdu <- svd(VAU)
      VAU <- svdu$u
      d <- svdu$d
      u = u %*% svdu$v
      
      
      ## Factorize n'th TT-core and merge to the next (either to the end, n-1, n+1)
      ## cf: VAU <- svdu$u, factorize and update v$G[[n]]
      ##     d <- svdu$d
      if (i==length(course)) {
        # End of the loop
        # Update V_n, d, u
        #svdu <- svd(VAU)
        #VAU <- svdu$u
        #d <- svdu$d
        #u = u %*% svdu$v
        dim(VAU) = c(v$R[n], v$J[n], v$R[n+1], K)  #Assume that K did not change#
        VAU = aperm(VAU, c(1,2,4,3))
        v$G[[n]] = VAU
        
        
      } else if (n < (n2 <- course[i+1])) { 
        # lr sweep
        # Update V_n, V_{n+1}
        dim(VAU) = c(v$R[n]*v$J[n], v$R[n+1]*K) 
        svdu = svd(VAU)
        Rnew = min( myChop2(svdu$d, eps), ttrank_max )
        
        VAU = svdu$u[,seq(Rnew)]
        dim(VAU) = c(v$R[n], v$J[n], Rnew)
        v$G[[n]] = VAU
          
        # Here, update V_{n+1} in order to compute Ares=Y-X accurately --
        cr2 = v$G[[n2]]
        sz2 = dim(cr2)
        dim(cr2) = c(sz2[1], prod(sz2[-1]))
        vs = svdu$v[,seq(Rnew),drop=FALSE] %*% diag(svdu$d[seq(Rnew)],Rnew)
        dim(vs) = c(v$R[n+1], K*Rnew )
        cr2 = t(vs) %*% cr2    #sz2[1]==R[n2]==R[n+1]
        dim(cr2) = c(K, Rnew, sz2[-1])  #K, Rnew, J[n2], R[n2+1]
        cr2 = aperm(cr2, c(2,3,1,4))
        v$G[[n2]] = cr2
        #----
        
        v$R[n+1] = Rnew
        v$block = n2
        
      } else {
        # rl sweep
        # Update V_n, V_{n-1}
        n2 <- course[i+1]  # We know n>course[i+1]
        
        VAU = t(VAU)
        dim(VAU) = c(K*v$R[n], v$J[n]*v$R[n+1]) 
        svdu = svd(VAU)
        Rnew = min( myChop2(svdu$d, eps), ttrank_max )
        
        VAU = svdu$v[,seq(Rnew)]
        VAU = t(VAU)
        dim(VAU) = c(Rnew, v$J[n], v$R[n+1])
        v$G[[n]] = VAU
        
        # Here, update V_{n+1} in order to compute Ares=Y-X accurately --
        cr2 = v$G[[n2]]
        sz2 = dim(cr2)
        dim(cr2) = c(prod(sz2[-length(sz2)]), sz2[length(sz2)])
        vs = svdu$u[,seq(Rnew),drop=FALSE] %*% diag(svdu$d[seq(Rnew)],Rnew)
        dim(vs) = c(K, v$R[n], Rnew )
        vs = aperm(vs, c(2,1,3))
        dim(vs) = c(v$R[n], K*Rnew)
        cr2 = cr2 %*% vs    #sz2[length(sz2)]==R[n]
        dim(cr2) = c(sz2[-length(sz2)], K, Rnew)  #R[n2], J[n2], K, Rnew
        v$G[[n2]] = cr2
        #--
        
        v$R[n] = Rnew
        v$block = n2
        
      }
    }
    
    # Ares = A_omega - X_omega
    Ares = sp_minus_lr_btt(Y, list(u=u,d=d,v=v))
       
    # Convergence trace and stopping
    ratio = lrRatio(u_old, d_old, v_old, u, d, v)
    obj = 0.5 * sum(Ares$v^2) + lambda * sum(d)
    obj_trace = c(obj_trace, obj)
    len_o = length(obj_trace)
    dif_o = abs(obj_trace[len_o] - obj_trace[len_o-1])/obj_trace[len_o]
    if (verbose)
      cat(iter, ":", "obj", obj, "obj_df", dif_o, "ratio", ratio, "\n")
    
    if (dif_o <= tol_df)
      break
    
    if (obj <= tol_f)
      break
    
    
  }
  if (iter == maxit && verbose) 
    warning(paste("Convergence not achieved by", maxit, 
                  "iterations"))
  if (lambda > 0 && final.svd) {
    #Let's update u again, because updating v_n may be is more complex. 
    u <- sp_times_btt(Ares, v) + u %*% diag(d,length(d)) 
    svdu <- svd(u)
    u <- svdu$u
    d <- svdu$d
    v = btt_times_mat(v, svdu$v)
    d = pmax(d - lambda, 0)  
    
    # Ares = A_omega - X_omega
    Ares = sp_minus_lr_btt(Y, list(u=u,d=d,v=v))
    obj = 0.5 * sum(Ares$v^2) + lambda * sum(d)
    obj_trace = c(obj_trace, obj)
    if (verbose)
      cat("final SVD:", "obj", format(round(obj, 5)), "\n")
  }  

  #rank_max = sum(d>0)
  X = list(u = u, d = d, v = v)
  attr(X,"call") <- this.call
  attr(X,"obj") <- obj_trace
  attr(X,"lambda") <- lambda
  return(X)
}