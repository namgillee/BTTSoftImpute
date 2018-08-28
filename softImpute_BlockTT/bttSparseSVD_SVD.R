bttSparseSVD_SVD <- function(Afill, ttrank_max, 
                             tol_u, lambda, maxit, verbose, 
                             X0 = NULL, final.svd = FALSE)
{
  #Low rank tensor approximation (SVD) based on block TT format.
  #
  #We use SVDS to optimize each core tensors simulaneously with 
  #left-singular-vectors(U), i.e., 
  #   (U,V1) -> (U,V2) -> ...
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
  #Last modified: 2017.12.31

  #require(rARPACK)
  
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
      
      ## Update U and V_n (n'th TT-core) by
      ##    [U, d, V_n]  =  svd( Afill * V_{neq} ), 
      ## where  Afill = {Ares + X} is the sum of sparse part and low-rank part
      vFrame = get_bttframe(v) #V_{neq} in BlockTT
      
      AV_sp = sp_times_btt(Afill$Ares, vFrame) #[I x RJR] full matrix

      AV_lr = Afill$X$u %*% (Afill$X$d * bttDot(Afill$X$v, vFrame) )

      AV = AV_sp + AV_lr     ##local matrix of size [I x RJR]
      
      
      ####svdu = svds(AV, K)     ##The K largest singular values
      svdu = svd(AV)
      
      u <- svdu$u[,seq(K), drop = FALSE]            ##Update the left singular vectors
      
      d <- pmax(svdu$d[seq(K)] - lambda, 0)      ##Update the singular values with soft thresholding
      
      cr <- svdu$v[,seq(K)]    ##V_n, of size [RJR x K]

            
      ## Factorize n'th TT-core and merge to the next (either end, rl, lr)
      ##
      if (i==length(course)) {
        # end of the loop
        
        dim(cr) = c(v$R[n], v$J[n], v$R[n+1], K)  ##Assume K did not change##
        cr = aperm(cr, c(1,2,4,3))
        v$G[[n]] = cr
        
        
      } else if (n < (n2 <- course[i+1])) { 
        # lr sweep
        
        dim(cr) = c(v$R[n]*v$J[n], v$R[n+1]*K) 
        svdv = svd(cr)
        Rnew = min( myChop2(svdv$d, eps), ttrank_max )
        
        cr = svdv$u[,seq(Rnew)]
        dim(cr) = c(v$R[n], v$J[n], Rnew)
        v$G[[n]] = cr    #G(n) is a third-order core tensor, G(n+1) is a 4th-order.

        # We skip updating the next core because it will 
        # be updated in the next iteration from the other cores
        #--
        # cr2 = v$G[[n2]]
        # sz2 = dim(cr2)
        # dim(cr2) = c(sz2[1], prod(sz2[-1]))
        # vs = svdv$v[,seq(Rnew)] %*% diag(svdv$d[seq(Rnew)])
        # dim(vs) = c(v$R[n+1], K*Rnew )
        # cr2 = t(vs) %*% cr2      #sz2[1]==R[n2]==R[n+1]
        # dim(cr2) = c(K, Rnew, sz2[-1])  #K, Rnew, J[n2], R[n2+1]
        # cr2 = aperm(cr2, c(2,3,1,4))
        # v$G[[n2]] = cr2
        
        v$R[n+1] = Rnew
        v$block = n2
        
      } else {
        # rl sweep
        
        n2 <- course[i+1]  # We know n>course[i+1]
        
        cr = t(cr)
        dim(cr) = c(K*v$R[n], v$J[n]*v$R[n+1]) 
        svdv = svd(cr)
        Rnew = min( myChop2(svdv$d, eps), ttrank_max )
        
        cr = svdv$v[,seq(Rnew)]
        cr = t(cr)
        dim(cr) = c(Rnew, v$J[n], v$R[n+1])
        v$G[[n]] = cr
        
        # Do not need to update d yet
        # AND
        # We skip updating the next core because it will 
        # be updated in the next iteration from the other cores
        #--
        # cr2 = v$G[[n2]]
        # sz2 = dim(cr2)
        # dim(cr2) = c(prod(sz2[-length(sz2)]), sz2[length(sz2)])
        # vs = svdv$u[,seq(Rnew)] %*% diag(svdv$d[seq(Rnew)])
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