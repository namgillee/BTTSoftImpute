# BTTFunctions
#
# Block Tensor Train - Related Functions


myChop2 <- function(d, ep)
{
  # ep is bet. 0 and 1

  if (ep <= 0)
  {
    r1 = length(d)
    return(r1)
  }

  dsqcum = cumsum(sort(d^2, decreasing = F))
  epnorm = ep * sqrt( dsqcum[length(dsqcum)] )
  ff = sum(dsqcum < epnorm^2) # number of small entries
  if (ff <= 0)
    r1 = length(d)
  else
    r1 = length(d)-ff

  return(max(1,r1)) # return
}
#
# ##################################
# ttTensor <- function(Y)
# {
#   if (is.null(Y)) {
#     myTT = list()
#     myTT$J = 0
#     myTT$N = 0
#     myTT$R = 0
#     myTT$G = vector("list",0)
#     myTT$block = 0
#     myTT$SP = list(i = matrix(0, 0, 0),
#                    v = matrix(0, 0, 0),
#                    nrow = 0,
#                    ncol = 0,
#                    dimnames = list())
#     return(myTT)
#   }
#
#   if (is.array(Y)) {
#     # Run TT-SVD
#     Jn = dim(Y)
#     N = length(Jn)    # Assume Jn is not NULL
#     Rn = rep(1, N+1)  # Assume N >= 1
#
#     # Initialize
#     K = max(Rn[1],Rn[N+1])
#     myTT = list()
#     myTT$J = Jn
#     myTT$N = N
#     myTT$R = Rn
#     myTT$G = vector("list",N)
#     myTT$block = N+1
#     myTT$SP = list(i = matrix(0, 0, N),
#                    v = matrix(0, 0, K),
#                    nrow = Jn,
#                    ncol = K,
#                    dimnames = list(list(), seq(K)))
#
#     # Update $R, $G from X
#     ep = .Machine$double.eps / sqrt(N-1)
#     X = Y   #### caution
#     n = 1
#     while (n <= N-1) {
#       #compute SVD
#       dim(X) <- c( Rn[n] * Jn[n],  length(X)/Rn[n]/Jn[n] )
#       svdX   <- svd(X)
#
#       #truncate rank
#       R2 <- myChop2(svdX$d, ep)
#       U1 <- svdX$u[,1:R2]
#       V1 <- svdX$v[,1:R2]
#       s1 <- svdX$d[1:R2]   ## Assume decreasing order in svdX$d
#
#       #update TT-core and TT-rank
#       dim(U1) <- c(Rn[n], Jn[n], R2)
#       myTT$G[[n]] = U1
#       myTT$R[n+1] = R2
#
#       #update for the next iterate
#       Rn[n+1] = R2
#       X <- s1 * t(V1)
#       n = n+1
#     }
#     dim(X) <- c(Rn[N], Jn[N], 1) ## (keep it 3rd order)
#     myTT$G[[N]] = X
#
#     return(myTT)
#   }
#
# }

##################################
#' @export
bttRand <- function(J, N, R, direction = 1)
{
  if ( length(J) == 1 ) {
    J = J * rep(1,N)
  }
  if ( length(R) == 1 ) {
    R = c(1, rep(R, N-1), 1)
  }
  if ( is.null(direction)) {
    direction = 1
  }

  # Initialize myTT
  # Set $block, $K, direction
  # $block is set at either 1 or N
  myTT = list()
  if (R[1]>1) {           ## This is the only case when $block==1 < N
    myTT$block = 1        ## In other cases, $block==N
    myTT$K = R[1]
    R[1] = 1
    if (direction>0) {
      warning('Perhaps lr qr is not suitable when R[1]>1. Changed to: rl')
      direction = -1
    }
  }
  else if (R[N+1]>1) {
    myTT$block = N
    myTT$K = R[N+1]
    R[N+1] = 1
    if (direction<0) {
      warning('Perhaps rl qr is not suitable when R[N+1]>1. Changed to: lr')
      direction = 1
    }
  }
  else {
    if (direction>0)
      myTT$block = N
    else
      myTT$block = 1
    myTT$K = 1
  }

  K = myTT$K
  myTT$J = J
  myTT$N = N
  myTT$R = R
  myTT$G = vector("list",N)



  # Update $R and $G
  # G[[block]] is 4th order, the others are 3rd order tensors
  if (direction>0) {
    # LR qr
    for (i in 1:N) {
      if (i < N) {  #i!= myTT$block
        cr1 = qr(matrix( rnorm(R[i]*J[i]*R[i+1]), R[i]*J[i] ))
        cr1 = qr.Q(cr1)

        R[i+1] = dim(cr1)[2]
        dim(cr1) = c(R[i], J[i], R[i+1])
        myTT$G[[i]] = cr1 # Update core
      }
      else { #i==N, R[i+1]==1
        cr1 = qr(matrix( rnorm(R[i]*J[i]*K*R[i+1]), R[i]*J[i] ))
        cr1 = qr.Q(cr1)

        K = dim(cr1)[2]
        dim(cr1) = c(R[i], J[i], K, 1)
        myTT$G[[i]] = cr1 # Update core
      }
    }
  } else {
    # RL lq
    for (i in seq(N,1)) {
      if (i>1) {  #i!= myTT$block
        cr1 = qr(matrix( rnorm(J[i]*R[i+1]*R[i]), J[i]*R[i+1] ))
        cr1 = qr.Q(cr1)
        cr1 = t(cr1)

        R[i] = dim(cr1)[1]
        dim(cr1) = c(R[i], J[i], R[i+1])
        myTT$G[[i]] = cr1 # Update core
      }
      else { #i==1, R[1]==1
        cr1 = qr(matrix( rnorm(J[i]*R[i+1]*R[i]*K), J[i]*R[i+1] ))
        cr1 = qr.Q(cr1)
        cr1 = t(cr1)

        K = dim(cr1)[1]
        dim(cr1) = c(1, K, J[i], R[i+1])
        cr1 = aperm(cr1, c(1,3,2,4))
        myTT$G[[i]] = cr1 # Update core
      }
    }
  }

  myTT$R = R
  return(myTT)
}

##################################
#' @export
bttDot <- function(tt1, tt2, do_qr=FALSE)
{
  # Dot product of two block TT tensors
  # It returns a matrix of size K(tt1) x K(tt2)
  # (During the process it computes a tensor of size
  #   R(1;tt1) x R(1;tt2) x K(tt1) x K(tt2) x R(N+1;tt1) x R(N+1;tt2)
  #   but we assume that R(1)=R(N+1)=1)

  N = tt1$N
  if (N != tt2$N) {
    warning('Dimensions of two input TT tensors are different')
    return(0)
  }
  J = tt2$J  #Assume tt2$J == tt1$J
  R1 = tt1$R
  R2 = tt2$R
  bk1 = tt1$block
  bk2 = tt2$block

  if ( do_qr ) {
    rv1 = bttQR(tt1)
    tt1 = rv1$Q
    rv1 = rv1$R

    rv2 = bttQR(tt2)
    tt2 = rv2$Q
    rv2 = rv2$R
  }

  K1 = R1[1]  #Be careful : K1 is not tt1$K, K2 is not tt2$K
  K2 = R2[1]
  p = diag(K1*K2)
  dim(p) = c( K1*K2*R1[1], R2[1] )

  for (i in 1:N) {
    cr1 = tt1$G[[i]]
    cr2 = tt2$G[[i]]

    if (i==bk1)
      K1n = dim(cr1)[3]
    else
      K1n = 1
    if (i==bk2)
      K2n = dim(cr2)[3]
    else
      K2n = 1

    dim(cr2) = c( R2[i], J[i]*K2n*R2[i+1] )

    p = p %*% cr2; # size k1*k2*r1-, J*K2n*r2+
    dim(p) = c( K1*K2, R1[i]*J[i], K2n*R2[i+1] )
    p = aperm(p, c(1, 3, 2))
    dim(p) = c( K1*K2*K2n*R2[i+1], R1[i]*J[i] )

    dim(cr1) = c( R1[i]*J[i], K1n*R1[i+1] )

    p = p %*% Conj(cr1) # size K1*K2*K2n*R2+, K1n*R1+
    dim(p) = c( K1, K2*K2n, R2[i+1], K1n, R1[i+1] )
    p = aperm(p, c(1,4,2,5,3))
    dim(p) = c( K1*K1n*K2*K2n*R1[i+1], R2[i+1] )

    K1 = K1*K1n
    K2 = K2*K2n
  }

  if (R1[N+1]==1 && R2[N+1]==1)
    dim(p) = c( K1, K2 ) #Remove R1[N+1], R2[N+1]
  else {
    dim(p) = c( K1, K2, R1[N+1], R2[N+1] )
    p = aperm(1,3,2,4)
    dim(p) = c( K1*R1[N+1], K2*R2[N+1] ) #Merge K1&R1[N+1], K2&R2[N+1]
  }

  if ( do_qr ) {
    p = t(rv1) %*% p %*% rv2
  }

  return(p)
}

##################################
bttQR <- function(myTT)
{
  bk = myTT$block
  N = myTT$N
  K = myTT$K
  R = myTT$R
  J = myTT$J

  # Left-to-right orthogonalization
  n = 1
  while (n<bk) {
    cr = myTT$G[[n]]
    dim(cr) = c(R[n]*J[n], R[n+1])

    #compute QR decomposition
    crR = qr(cr)
    cr = qr.Q(crR)
    crR = qr.R(crR)

    #update TT-rank
    R[n+1] = dim(crR)[1]

    #update current core
    dim(cr) = c(R[n], J[n], R[n+1])
    myTT$G[[n]] = cr

    #multiply crR to the next core
    cr2 = myTT$G[[n+1]]
    sz2 = dim(cr2)
    dim(cr2) = c(sz2[1], prod(sz2[-1]))
    cr2 = crR %*% cr2  #sz2[1] == R[n+1]
    dim(cr2) = c(R[n+1], sz2[-1])
    myTT$G[[n+1]] = cr2

    n = n+1
  }

  # Right-to-left orthogonalization
  n = N
  while (n>bk) {
    cr = myTT$G[[n]]
    dim(cr) = c(R[n], J[n]*R[n+1])
    cr = t(cr)  #size (J[n]*R[n+1], R[n])

    #compute QR decomposition
    crR = qr(cr)
    cr = qr.Q(crR)
    crR = qr.R(crR)

    #update TT-rank
    R[n] = dim(crR)[1]

    #update current core
    cr = t(cr)
    dim(cr) = c(R[n], J[n], R[n+1])
    myTT$G[[n]] = cr

    #multiply crR to the next core
    cr2 = myTT$G[[n-1]]
    sz2 = dim(cr2)
    dim(cr2) = c(prod(sz2[-length(sz2)]), sz2[length(sz2)])
    cr2 = cr2 %*% t(crR)  #sz2[4] == R[n]
    dim(cr2) = c(sz2[-length(sz2)], R[n])
    myTT$G[[n-1]] = cr2

    n = n-1
  }

  ## Orthogonalize bk'th TT-core ##
  cr = myTT$G[[bk]]
  dim(cr) = c(R[bk]*J[bk], K, R[bk+1])
  cr = aperm(cr, c(1,3,2))
  dim(cr) = c(R[bk]*J[bk]*R[bk+1], K)

  crR = qr(cr)
  cr = qr.Q(crR)
  crR = qr.R(crR) ### We will return this qrR matrix

  #if (K!=dim(crR)[1])
  #  warning("In bttDot(): column size myTT$K has changed")
  K = dim(crR)[1]

  #update current core
  dim(cr) = c(R[bk], J[bk], R[bk+1], K)
  cr = aperm(cr, c(1,2,4,3))
  myTT$G[[bk]] = cr


  ## Update myTT ##
  myTT$R = R
  myTT$K = K

  return( list(Q=myTT, R=crR) )
}

#' @export
bttFull <- function(V)
{
  # Convert blockt TT tensor V into block matrix
  # of size [J{1}*J{2}*...*J{N} x K]
  #
  # Inputs:
  #   V       : block TT format with attributes N, R, J, G, block

  N = V$N
  R = V$R
  J = V$J
  K = V$K
  bk = V$block

  if (N == 0) {
    warning("Input tensor is empty")
    return(matrix(0, nrow = 0, ncol = K))
  }

  ## Main loop
  z = diag(1,R[N+1])
  for (n in N:1) {

    cr = V$G[[n]] #RxJxR

    if (n == bk) {
      cr = aperm(cr, c(3,1,2,4)) #KxRxJxR
    }

    dim(cr) <- c(length(cr)/R[n+1], R[n+1])
    z = cr %*% z #KxRxJJ..JxRxK'

    if (n == bk) {
      dim(z) <- c(K, length(z)/K)
      z = t(z) #RxJJ..JxRxK
    }

    dim(z) <- c(R[n], length(z)/R[n])
  }

  # Finally, z has size R{1} x J{1}*J{2}*...*J{N}*R{N+1}*K
  dim(z) = c(R[1]*prod(J)*R[N+1], K)

  return(z)

}
