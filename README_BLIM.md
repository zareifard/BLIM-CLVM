# BLIM

Bayesian latent-variable model for ordinal observations, with a regression component for the latent response and an item-level probability component. This README documents the current implementation behavior described for the repository; it does not imply that the short example is sufficient for convergence assessment.

## Dependencies

The implementation uses `MASS`, `truncnorm`, `boa`, and `bayesm` (including `rigamma`). Install the packages if needed:

```r
install.packages(c("MASS", "truncnorm", "boa", "bayesm"))
```

Load the function by sourcing its R file (adjust the path to your repository layout):

```r
source("BLIM.R")
```

## Function

```r
BLIM(Y, F, mc.A, mc.B, v.p0, L, U, m.b, c.b, m.bz, c.bz,
     Up.z, Ls, Us, simulate, burnin, Break)
```

Use named arguments when calling the function. `Break` is the spacing, in iterations, between saved draws after burn-in; it is not an independent MCMC chain or a convergence diagnostic.

| Argument | Meaning and requirements |
|---|---|
| `Y` | Numeric/integer-valued response vector of length `n`. Current category/indexing logic is designed for integer categories; ensure the observed support is consecutive from `K=min(Y)` through `Ma=max(Y)` (no missing categories). |
| `F` | Numeric design matrix with `n` rows and `p` columns, full column rank. The code uses `solve(t(F) %*% F)`, so rank deficiency or incompatible dimensions will fail. Include an intercept column explicitly if desired. |
| `mc.A`, `mc.B` | Numeric vectors of length `n`; their entries and resulting denominator `mc + mc.t` must be positive. |
| `v.p0` | Scalar variance parameter for the beta probability prior. For every observation, use a value satisfying `0 < v.p0 < E.p0 * (1-E.p0)` after clamping, so the derived beta shape parameters `al` and `be` are positive. Choose `L` and `U` compatibly. |
| `L`, `U` | Lower and upper bounds used to clamp `E.p0`, with `L <= U`; practical use requires bounds compatible with valid beta parameters. |
| `m.b` | Numeric prior-mean vector of length `p` for regression coefficients. |
| `c.b` | Scalar diagonal prior variance parameter, used through `I/c.b`; use a positive value. |
| `m.bz` | Scalar prior mean for `z`/the latent scale parameter. |
| `c.bz` | Scalar diagonal prior variance parameter for `z`, used through `I/c.bz`; use a positive value. |
| `Up.z` | Upper truncation bound for `z`; its sampling interval is `[0, Up.z]`. Use `Up.z > 0`. |
| `Ls`, `Us` | Lower and upper bounds for `sigma2`; sampling is constrained to `[Ls, Us]`. Use compatible bounds, typically `0 < Ls < Us`. |
| `simulate` | Total number of iterations. Must be an integer greater than `burnin` for draws to be saved. |
| `burnin` | Number of initial iterations discarded. The implementation does not validate controls. |
| `Break` | Positive integer spacing between saved iterations after burn-in. The implementation does not validate controls. |

The implementation computes `n = length(Y)`, `mc = mc.A + mc.B`, `mc.t = mean(mc)`, and `E.p0 = mc/(mc + mc.t)`, then clamps `E.p0` to `[L,U]`. It constructs `bound = c(-Inf, K:Ma)` and indexes it using `Y`. Accordingly, use integer-valued consecutive categories and check the support before fitting.

## Reproducible example

This is a runnable synthetic-data call example. Its deliberately short chain is for demonstrating inputs and outputs only, not for inference or convergence claims.

```r
set.seed(20250308)
n <- 60
x <- as.numeric(scale(rnorm(n)))
F <- cbind(intercept = 1, x = x)
Y <- sample(3:8, size = n, replace = TRUE)
mc.A <- runif(n, 50, 150)
mc.B <- runif(n, 50, 150)

fit <- BLIM(
  Y = Y, F = F,
  mc.A = mc.A, mc.B = mc.B,
  v.p0 = 0.05, L = 0.05, U = 0.95,
  m.b = rep(0, ncol(F)), c.b = 100,
  m.bz = 0, c.bz = 100, Up.z = 10,
  Ls = 0.01, Us = 100,
  simulate = 200, burnin = 50, Break = 5
)
str(fit)
```

The generated `mc.A` and `mc.B` are positive, and their resulting `E.p0` values are close to 0.5; `v.p0=0.05` gives positive beta shape parameters for this example.

## Outputs

The function returns the following quantities (inspect `str(fit)` for the exact list/container structure in the sourced version):

| Field | Description |
|---|---|
| `beta` | Stored draws for regression coefficients (`D.beta`). |
| `alpha` | Stored draws for the scalar probability/latent parameter. |
| `sig2` | Stored draws for `sigma2`. |
| `X`, `Z` | Stored latent-variable draws, as matrices. |
| `E.beta` | Column-wise `summary` of `D.beta` (`apply(D.beta, 2, summary)`); this is a summary matrix, not merely a posterior mean. |
| `V.beta` | Column-wise standard deviations of coefficient draws. |
| `E.beta.z` | `summary(alpha)`. |
| `V.beta.z` | Standard deviation of `alpha`. |
| `HPD.beta.z` | 95% highest-posterior-density interval from `boa.hpd(..., 0.05)`. |
| `E.sig2` | `summary(sig2)`. |
| `V.sig2` | Standard deviation of `sig2`. |
| `HPD.sig2` | 95% HPD interval from `boa.hpd(..., 0.05)`. |
| `E.X`, `E.Z` | Column means of the stored `X` and `Z` matrices. |

## MCMC storage and caveats

After `burnin`, the implementation stores a draw every `Break` iterations. With integer controls, `simulate > burnin`, and positive `Break`, the expected number of saved draws is `floor((simulate-burnin)/Break)`. The code does not validate these controls, so provide sensible integer values yourself.

The function computes the beta-prior shapes from the clamped `E.p0` and `v.p0`; invalid shapes (`al <= 0` or `be <= 0`) are rejected by the code. Ensure `F` is full rank, all vector/matrix dimensions agree, category support is consecutive, and variance/bound parameters are valid. The example only demonstrates a call. Choose chain length based on the application and assess convergence, mixing, and sensitivity before interpreting posterior summaries.
