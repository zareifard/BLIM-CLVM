# =============================================================================
# Binary Latent Indicator Model (BLIM)
# =============================================================================

library(MASS)
library(truncnorm)
library(boa)
library(bayesm)

# =============================================================================
# Truncated Gamma sampler
# =============================================================================

rtigamma <- function(s,r,L,U){
  g <- 0
  while(g < L | g > U) g <- rigamma(1,s,r)
  return(g)
}


# =============================================================================
# P(a < X < b), X ~ N(m,s^2)
# =============================================================================

Pab <- function(a,b,m,s){
  pn <- pnorm(b,m,s)-pnorm(a,m,s)
  return(pn)
}

# =============================================================================

BLIM <- function(Y,F,mc.A,mc.B,L,U,v.p0=0.2,
                 m.b,c.b,m.bz,c.bz,Up.z,Ls,Us,
                 simulate,burnin,Break){

  # ---------------------------------------------------------------------------
  # Data and initial quantities
  # ---------------------------------------------------------------------------

  n <- length(Y)
  Lev <- nrow(as.matrix(table(Y)))
  K <- min(Y); Ma <- max(Y)
  bound <- c(-Inf,K:Ma)

  mc <- mc.A+mc.B; mc.t <- mean(mc)

  X <- Y
  Z <- W <- c()

  E.p0 <- mc/(mc+mc.t)
  E.p0[E.p0 < L] <- L
  E.p0[E.p0 > U] <- U


  # ---------------------------------------------------------------------------
  # Initial values of latent indicators and Beta prior
  # ---------------------------------------------------------------------------

  for(i in 1:n) Z[i] <- ifelse(E.p0[i] >= .5,1,0)

  al <- (E.p0*(1-E.p0)/v.p0-1)*E.p0
  be <- (1-E.p0)/E.p0*al

  if(any(al <= 0) || any(be <= 0))
    stop("Invalid Beta prior parameters. Check L, U and v.p0.")


  # ---------------------------------------------------------------------------
  # Initial values
  # ---------------------------------------------------------------------------

  beta <- as.vector(solve(t(F)%*%F)%*%t(F)%*%Y)
  F.B <- as.vector(F%*%beta); FF <- t(F)%*%F

  beta.z <- 1
  sig2 <- var(Y)


  # ---------------------------------------------------------------------------
  # Storage
  # ---------------------------------------------------------------------------

  D.beta.z <- c(); D.beta <- c(); D.sig <- c()
  D.X <- c(); D.Z <- c()

  R1 <- 0; sim <- 0
  start_time <- Sys.time()


  # ===========================================================================
  # Gibbs sampler
  # ===========================================================================

  for(R in 1:simulate){

    # -------------------------------------------------------------------------
    # Full conditional distributions of Z and W
    # -------------------------------------------------------------------------

    for(i in 1:n){
      for(j in K:Ma){
        if(Y[i] == j){
          P <- Pab(bound[j-K+1],bound[j-K+2],F.B[i]+beta.z,sqrt(sig2))
          Q <- Pab(bound[j-K+1],bound[j-K+2],F.B[i],sqrt(sig2))
        }
      }

      pw <- P*al[i]/(P*al[i]+Q*be[i])

      W[i] <- ifelse(rbinom(1,1,pw) == 1,
                     rbeta(1,al[i]+1,be[i]),
                     rbeta(1,al[i],be[i]+1))

      Z[i] <- rbinom(1,1,P*W[i]/(P*W[i]+Q*(1-W[i])))
    }


    # -------------------------------------------------------------------------
    # Full conditional distribution of X
    # -------------------------------------------------------------------------

    for(i in 1:n){
      for(j in K:Ma){
        if(Y[i] == j)
          X[i] <- rtruncnorm(1,bound[j-K+1],bound[j-K+2],
                             F.B[i]+beta.z*Z[i],sqrt(sig2))
      }
    }


    # -------------------------------------------------------------------------
    # Full conditional distribution of beta
    # -------------------------------------------------------------------------

    S <- solve(diag(ncol(F))/c.b+FF/sig2)
    M <- S%*%(t(F)%*%(X-beta.z*Z)/sig2+m.b/c.b)
    beta <- as.vector(mvrnorm(1,as.vector(M),S))

    F.B <- as.vector(F%*%beta)


    # -------------------------------------------------------------------------
    # Full conditional distribution of beta.z (alpha)
    # -------------------------------------------------------------------------

    S <- as.numeric((1/c.bz+sum(Z)/sig2)^(-1))
    M <- S%*%(t(Z)%*%(X-F.B)/sig2+m.bz/c.bz)

    beta.z <- rtruncnorm(1,0,Up.z,M,sqrt(S))


    # -------------------------------------------------------------------------
    # Full conditional distribution of sigma^2
    # -------------------------------------------------------------------------

    xx.s <- .5*sum((X-F.B-beta.z*Z)^2)
    sig2 <- rtigamma((n-1)/2,xx.s,Ls,Us)


    # -------------------------------------------------------------------------
    # Save posterior samples
    # -------------------------------------------------------------------------

    if(R > burnin){
      R1 <- R1+1

      if(R1 == Break){
        sim <- sim+1

        D.beta <- rbind(D.beta,beta)
        D.beta.z <- c(D.beta.z,beta.z)
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

  E.beta.z <- summary(D.beta.z)
  V.beta.z <- sd(D.beta.z)
  HPD.beta.z <- boa.hpd(D.beta.z,.05)

  E.sig2 <- summary(D.sig)
  V.sig2 <- sd(D.sig)
  HPD.sig2 <- boa.hpd(D.sig,.05)

  E.X <- apply(D.X,2,mean)
  E.Z <- apply(D.Z,2,mean)


  # ---------------------------------------------------------------------------
  # Computational time
  # ---------------------------------------------------------------------------

  end_time <- Sys.time()
  time <- end_time-start_time

  cat("time ====>",time,"\n")


  # ---------------------------------------------------------------------------
  # Return posterior samples and summaries
  # ---------------------------------------------------------------------------

  return(list(
    beta=D.beta,
    alpha=D.beta.z,
    sig2=D.sig,
    X=D.X,
    Z=D.Z,
    E.beta=E.beta,
    V.beta=V.beta,
    E.beta.z=E.beta.z,
    V.beta.z=V.beta.z,
    HPD.beta.z=HPD.beta.z,
    E.sig2=E.sig2,
    V.sig2=V.sig2,
    HPD.sig2=HPD.sig2,
    E.X=E.X,
    E.Z=E.Z
  ))
}