# Continuous Latent Variable Model (CLVM)

This repository provides an R implementation of the Continuous Latent
Variable Model (CLVM) for modeling interaction frequencies in ChIA-PET
data.

# FUNCTION

CLVM(Y,F,mc.A,mc.B,L,U,
     m.b,c.b,m.n,c.n,Un,Ls,Us,
     simulate,burnin,Break)

# ARGUMENTS

Y
Observed interaction frequencies y_i, i = 1,...,n.

F
Design matrix of covariates. The first column is typically an intercept.

mc.A
Marginal PET counts for anchor A_i.

mc.B
Marginal PET counts for anchor B_i.

L
Lower bound for the prior mean of the mixture/signal probability.

U
Upper bound for the prior mean of the mixture/signal probability.

m.b
Prior mean vector for beta.
The value should be selected according to available prior information.

c.b
Prior variance parameter for beta.
The value should be selected according to prior information. Large values
correspond to more diffuse priors.

m.n
Prior mean of nu.
The prior for nu is truncated to the positive domain. 

c.n
Prior variance parameter of nu. Large values correspond to a more diffuse prior.

Un
Upper truncation bound for the prior of nu.

Ls
Lower truncation bound for sigma^2.

Us
Upper truncation bound for sigma^2.

simulate
Total number of Gibbs sampling iterations.

burnin
Number of initial iterations discarded as burn-in.

Break
Number of iterations between two successive saved posterior samples.

