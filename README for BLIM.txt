
Binary Latent Indicator Model (BLIM)
========================================

This repository provides an R implementation of the Binary Latent Indicator Model (BLIM) for modeling interaction frequencies in
ChIA-PET data.


FUNCTION
========

BLIM <- function(Y,F,mc.A,mc.B,L,U,v.p0=0.2,
                 m.b,c.b,m.bz,c.bz,Up.z,Ls,Us,
                 simulate,burnin,Break)


ARGUMENTS
=========

Y
Observed interaction frequencies y_i, i = 1,...,n.

F
Design matrix of covariates. The first column is typically an intercept.

mc.A
Marginal PET counts for anchor A_i.

mc.B
Marginal PET counts for anchor B_i.

L
Lower bound for the prior mean of the mixture weights.

U
Upper bound for the prior mean of the mixture weights.

v.p0
Prior variance used to determine the Beta prior parameters of the mixture
weights omega_i. 

m.b
Prior mean of beta. 

c.b
Prior variance parameter for beta. It should be determined according to
prior information. Large values correspond to diffuse priors.

m.bz
Prior mean of alpha. 

c.bz
Prior variance parameter of alpha. It should be determined according to
prior information. Large values correspond to diffuse priors.

Up.z
Upper truncation point for alpha. 

Ls
Lower truncation bound for sigma^2. 

Us
Upper truncation bound for sigma^2. 

simulate
Total number of MCMC iterations.

burnin
Number of initial MCMC iterations discarded as burn-in.

Break
Number of iterations between two successive saved posterior samples.
