# =============================================================================
# Continuous Latent Variable Model (CLVM)
# =============================================================================

library(MASS)
library(truncnorm)
library(boa)
library(bayesm)


# =============================================================================
# Truncated Inverse-Gamma sampler
# =============================================================================

rtigamma <- function(s,r,L,U){
  g <- 0
  while(g < L | g > U) g <- rigamma(1,s,r)
  return(g)
}

# =============================================================================
# Continuous Latent Variable Model
# =============================================================================

CLVM <- function(Y,F,mc.A,mc.B,L,U,
                 m.b,c.b,m.n,c.n,Un,Ls,Us,
                 simulate,burnin,Break){

  # ---------------------------------------------------------------------------
  # Data and initial quantities
  # ---------------------------------------------------------------------------

  n <- length(Y)
  K <- min(Y)
  Ma <- max(Y)
  bound <- c(-Inf,K:Ma)

  mc <- mc.A+mc.B
  mc.t <- mean(mc)

  X <- Y
  Z <- W <- c()

  E.p0 <- mc/(mc+mc.t)
  E.p0[E.p0 < L] <- L
  E.p0[E.p0 > U] <- U


  # ---------------------------------------------------------------------------
  # Initial values of latent variables
  # ---------------------------------------------------------------------------

  for(i in 1:n) Z[i] <- ifelse(E.p0[i] >= .5,1,0)

  beta <- as.vector(solve(t(F)%*%F)%*%t(F)%*%Y)
  F.B <- as.vector(F%*%beta)
  FF <- t(F)%*%F

  nu <- 1
  sig2 <- var(Y)

  H <- X-F.B
  F.Bp <- qnorm(E.p0)


  # ---------------------------------------------------------------------------
  # Storage
  # ---------------------------------------------------------------------------

  D.beta <- c()
  D.nu <- c()
  D.sig <- c()
  D.X <- c()
  D.Z <- c()

  R1 <- 0
  sim <- 0

  start_time <- Sys.time()


  # ===========================================================================
  # Gibbs sampler
  # ===========================================================================

  for(R in 1:simulate){

    # -------------------------------------------------------------------------
    # Full conditional of Y* = X
    # -------------------------------------------------------------------------

    for(i in 1:n){
      if(Y[i] == 3)  X[i] <- rtruncnorm(1,-Inf,3,F.B[i]+nu*H[i]*Z[i],sqrt(sig2))
      if(Y[i] > 3)   X[i] <- rtruncnorm(1,Y[i]-1,Y[i],F.B[i]+nu*H[i]*Z[i],sqrt(sig2))
    }

    ep <- X-F.B


    # -------------------------------------------------------------------------
    # Full conditional of Z and H
    # -------------------------------------------------------------------------

    v.h <- sig2/(nu^2+sig2)
    mu.h <- v.h*(F.Bp+ep*nu/sig2)

    C <- sqrt(v.h)*exp(.5*mu.h^2/v.h-.5*F.Bp^2)
    P <- pnorm(mu.h/sqrt(v.h))
    Q <- pnorm(-F.Bp)

    Z <- rbinom(n,1,P*C/(P*C+Q))

    for(i in 1:n){
      if(Z[i] == 1) H[i] <- rtruncnorm(1,0,Inf,mu.h[i],sqrt(v.h))
      if(Z[i] == 0) H[i] <- rtruncnorm(1,-Inf,0,F.Bp[i],1)
    }


    # -------------------------------------------------------------------------
    # Full conditional of beta
    # -------------------------------------------------------------------------

    S <- solve(diag(ncol(F))/c.b+FF/sig2)
    M <- S%*%(t(F)%*%(X-nu*H*Z)/sig2+m.b/c.b)

    beta <- as.vector(mvrnorm(1,as.vector(M),S))

    F.B <- as.vector(F%*%beta)
    ep <- X-F.B


    # -------------------------------------------------------------------------
    # Full conditional of nu
    # -------------------------------------------------------------------------

    S <- (1/c.n+sum(Z*H^2)/sig2)^(-1)
    M <- S*(sum((Z*H)*ep)/sig2+m.n/c.n)

    nu <- rtruncnorm(1,0,Un,M,sqrt(S))


    # -------------------------------------------------------------------------
    # Full conditional of sigma^2
    # -------------------------------------------------------------------------

    xx.s <- .5*sum((ep-nu*H*Z)^2)
    sig2 <- rtigamma((n-1)/2,xx.s,Ls,Us)


    # -------------------------------------------------------------------------
    # Save posterior samples
    # -------------------------------------------------------------------------

    if(R > burnin){
      R1 <- R1+1

      if(R1 == Break){
        sim <- sim+1

        D.beta <- rbind(D.beta,beta)
        D.nu <- c(D.nu,nu)
        D.sig <- c(D.sig,sig2)
        D.X <- rbind(D.X,X)
        D.Z <- rbind(D.Z,Z)

        R1 <- 0
      }
    }
  }


  # ===========================================================================
  # Posterior summaries
  # ===========================================================================

  E.beta <- apply(D.beta,2,summary)
  V.beta <- apply(D.beta,2,sd)

  E.nu <- summary(D.nu)
  V.nu <- sd(D.nu)
  HPD.nu <- boa.hpd(D.nu,.05)

  E.sig2 <- summary(D.sig)
  V.sig2 <- sd(D.sig)
  HPD.sig2 <- boa.hpd(D.sig,.05)

  E.X <- apply(D.X,2,mean)
  E.Z <- apply(D.Z,2,mean)


  # ===========================================================================
  # Computational time
  # ===========================================================================

  end_time <- Sys.time()
  time <- end_time-start_time

  cat("time ====>",time,"\n")


  # ===========================================================================
  # Return posterior samples and summaries
  # ===========================================================================

  return(list(
    beta=D.beta,
    nu=D.nu,
    sig2=D.sig,
    X=D.X,
    Z=D.Z,
    E.beta=E.beta,
    V.beta=V.beta,
    E.nu=E.nu,
    V.nu=V.nu,
    HPD.nu=HPD.nu,
    E.sig2=E.sig2,
    V.sig2=V.sig2,
    HPD.sig2=HPD.sig2,
    E.X=E.X,
    E.Z=E.Z
  ))
}