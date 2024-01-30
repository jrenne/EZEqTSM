solve.model <- function(Model){
  # This procedure solves for the s.d.f.
  theta <- (1 - Model$gamma)/(1 - 1/Model$psi)
  n <- dim(Model$Pi)[1]
  q <- dim(Model$nu)[1]

  pc.bar <- 5
  eps <- .0000001

  for(i in 1:10){

    for(j in 1:2){
      # when j==1, baseline,
      # when j==2, perturbation of pc.bar
      if(j==1){
        pc.bar.j <- pc.bar
      }else{
        pc.bar.j <- pc.bar + eps
      }

      kappa.1 <- exp(pc.bar.j)/(1 + exp(pc.bar.j))
      kappa.0 <- log(1 + exp(pc.bar.j)) - kappa.1 * pc.bar.j

      A.1x <- (1 - 1/Model$psi) * solve(diag(n) - kappa.1*t(Model$Pi)) %*% t(Model$Pi) %*% Model$Beta[,2]
      A.1z <- theta/2*solve(diag(q) - kappa.1 * t(Model$nu)) %*% t(Model$nu) %*%
        (Model$chi.1 %*% (solve(t(Model$Pi)) %*% A.1x)^2 +
           (1 - Model$gamma) * Model$Alpha[,2]/theta^2)

      A.1 <- matrix(c(A.1x,A.1z,rep(0,3)),ncol=1)

      beta.c.tilde <- matrix(c(Model$Beta[,2],Model$Alpha[,2],0,1,0),ncol=1)
      aux <- theta * kappa.1 * A.1 + (1 - Model$gamma) * beta.c.tilde
      AUX <- compute.Laplace.Transf.X(Model,aux)
      phi.0 <- AUX$phi.0

      A.0 <- 1/(theta * (1 - kappa.1)) * ( theta * log(Model$delta) + theta * kappa.0 +
                                             (1 - Model$gamma) * Model$mu[2] + phi.0)
      if(j==1){
        f.1 <- A.0 - pc.bar.j
      }else{
        f.2 <- A.0 - pc.bar.j
      }
    }
    df <- (f.2 - f.1)/eps

    pc.bar <- c(pc.bar - f.1/df)
  }

  Model$lambda <- (theta - 1)*kappa.1*A.1 - Model$gamma*beta.c.tilde

  phi.lambda <- compute.Laplace.Transf.X(Model,matrix(Model$lambda,ncol=1))

  Model$eta.0  <- - theta*log(Model$delta) + Model$gamma * Model$mu[2] +
    (1 - theta)*kappa.0 + (1-theta)*A.0*(kappa.1 - 1) - phi.lambda$phi.0
  Model$eta.1 <- (theta - 1)*A.1 - phi.lambda$phi.1

  Model.solved <- Model
  Model.solved$A.0 <- A.0
  Model.solved$A.1 <- A.1
  Model.solved$pc.bar <- pc.bar
  Model.solved$kappa.0 <- kappa.0
  Model.solved$kappa.1 <- kappa.1

  Model.solved$dev.in.GN <- abs(f.1)

  return(Model.solved)
}

compute.prices <- function(Model,W,vec.H){
  # vec.H is a vector indicating the considered maturities (nb.mat)
  # W is of dimension (n+q+3)xK, where K is the number of considered types of bonds
  # outputs:
  #    a.h: dimension K x nb.mat
  #    b.h: dimension (n + q + 3) x K x nb.mat (i.e. class = array)

  nb.mat <- length(vec.H)
  n <- dim(Model$Pi)[1]
  q <- dim(Model$nu)[1]

  K <- dim(W)[2]

  a <- matrix(0,K,1)
  b <- matrix(0,n+q+3,K)

  all.a <- matrix(NaN,K,nb.mat)
  all.b <- array(NaN,c(n+q+3,K,nb.mat))
  all.a.return <- matrix(NaN,K,nb.mat)
  all.b.return <- array(NaN,c(n+q+3,K,nb.mat))

  Lambda <- Model$lambda %*% matrix(1,1,K)

  count <- 0
  for( i in 1:max(vec.H)){
    phi.lambda <- compute.Laplace.Transf.X(Model,Lambda,X=NaN)

    beta.c.tilde <-  matrix(c(Model$Beta[,2],Model$Alpha[,2],0,1,0),n+q+3,K)
    beta.pi.tilde <- matrix(c(Model$Beta[,3],Model$Alpha[,3],0,0,1),n+q+3,K)

    aux <- Lambda + b + W - beta.pi.tilde

    phi.aux <- compute.Laplace.Transf.X(Model,aux,X=NaN)

    a <- a - Model$mu[3] - matrix(Model$eta.0,K,1) + phi.aux$phi.0 - phi.lambda$phi.0

    b <- phi.aux$phi.1 - phi.lambda$phi.1 - matrix(Model$eta.1,n+q+3,K)

    if(sum(i==vec.H)>0){
      count <- count + 1
      all.a[,count] <- a
      all.b[,,count] <- b
      all.a.return[,count] <- -1/vec.H[count] * a
      all.b.return[,,count] <- -1/vec.H[count] * b
    }
  }
  return(list(
    all.a = all.a,
    all.b = all.b,
    all.a.return = all.a.return,
    all.b.return = all.b.return
  ))
}

compute.Laplace.Transf.X <- function(Model,u,X=NaN){
  n <- dim(Model$Pi)[1]
  q <- dim(Model$nu)[1]
  K <- dim(u)[2]
  # u has to be of dimension (n + q + 3) x K (if ones wants to compute the Laplace Transform for K vectors)
  u.x  <- matrix(u[1:n,],nrow=n)
  u.z  <- matrix(u[(n+1):(n+q),],nrow=q)
  u.epsilon  <- matrix(u[(n+q+1):(n+q+3),],nrow=3)

  phi.0 <- .5 * t(u.x * u.x) %*% matrix(Model$chi.0,ncol=1) +
    .5 * t(u.z * u.z) %*% matrix(Model$sigma.w^2,ncol=1) +
    .5 * ((t(u.epsilon) %*% Model$Omega) * t(u.epsilon)) %*% matrix(1,3,1) # This is 1/2 u.eps' Omega u.eps

  phi.1 <- rbind(
    t(Model$Pi) %*% u.x,
    .5 * t(Model$nu) %*% Model$chi.1 %*% (u.x * u.x) + t(Model$nu) %*% u.z,
    matrix(0,3,K)
  )

  if(!is.na(X)){
    LT <- exp(phi.0 + t(phi.1) %*% X)
  }else{
    LT <- NaN
  }

  return(list(
    phi.0 = phi.0,
    phi.1 = phi.1,
    LT = LT
  ))
}

compute.cov.mat.X <- function(Model){
  # This function computes the conditional and marginal covariance matrices of X,
  #       where X = [x',z',eps']'

  n <- dim(Model$Pi)[1] # dimension of x
  q <- dim(Model$nu)[1] # dimension of z

  big.Pi <- matrix(0,n+q+3,n+q+3)
  big.Pi[1:n,1:n] <- Model$Pi
  big.Pi[(n+1):(n+q),(n+1):(n+q)] <- Model$nu

  big.Sigma <- matrix(0,n+q+3,n+q+3)
  big.Sigma[1:n,1:n] <- diag(Model$chi.0)
  big.Sigma[(n+1):(n+q),(n+1):(n+q)] <- diag((Model$sigma.w)^2)
  big.Sigma[(n+q+1):(n+q+3),(n+q+1):(n+q+3)] <- Model$Omega

  vec.uncond.cov.matrix <- solve(diag(n+q+3)%x%diag(n+q+3) - big.Pi %x% big.Pi) %*% c(big.Sigma)
  uncond.cov.matrix <- matrix(vec.uncond.cov.matrix,n+q+3,n+q+3)
  uncond.autocov.matrix <- big.Pi %*% uncond.cov.matrix # This is the covariance matrix of X_{t},X_{t-1}

  # Elements regarding the conditional covariance matrix of X_{t+1}:
  Gamma_0 <- diag(c(Model$chi.0,Model$sigma.w^2,rep(0,3)))
  Gamma_0[(n+q+1):(n+q+3),(n+q+1):(n+q+3)] <- Model$Omega
  Gamma_1 <- matrix(0,n+q+3,n+q+3)
  Gamma_1[1:n,(n+1):(n+q)] <- t(Model$chi.1) %*% Model$nu

  return(list(
    uncond.cov.matrix=uncond.cov.matrix,
    uncond.autocov.matrix=uncond.autocov.matrix,
    big.Pi=big.Pi,
    big.Sigma=big.Sigma,
    Gamma_0=Gamma_0,
    Gamma_1=Gamma_1))
}
