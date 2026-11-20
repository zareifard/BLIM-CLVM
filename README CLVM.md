# Continuous Latent Variable Model (CLVM)

CLVM is an R implementation of a Bayesian continuous latent-variable model for interaction frequencies in ChIA-PET data. The model is applicable to both intra- and inter-chromosomal interactions; the GM12878 example below focuses on **intra-chromosomal interactions**.

## Installation and package loading

Install the required R packages once, then load them in each R session:

```r
install.packages(c("MASS", "truncnorm", "boa", "bayesm"))
library(MASS)
library(truncnorm)
library(boa)
library(bayesm)
```

From the repository root, source the implementation:

```r
source("CLVM.R")
```

## Function

```r
CLVM <- function(Y,F,mc.A,mc.B,L,U,m.b,c.b,m.n,c.n,Un,Ls,Us,simulate,burnin,Break)
```

There are no default arguments; supply every argument when calling the function.

| Argument | Description |
|---|---|
| `Y` | Vector of observed interaction frequencies, `y_i`. In this implementation, the sampler updates `X` for observations with `Y == 3` or `Y > 3`; use integer values of at least 3. |
| `F` | Covariate design matrix with one row per observation. Its first column is usually an intercept. It must have full column rank. |
| `mc.A`, `mc.B` | Marginal PET counts for anchors A and B, respectively. |
| `L`, `U` | Lower and upper bounds used to clamp the initial/prior signal probability calculated from the marginal counts. |
| `m.b`, `c.b` | Prior mean vector for regression coefficients `beta` and scalar covariance scale, respectively. `c.b` is used in the prior covariance `I/c.b`. |
| `m.n`, `c.n` | Prior mean and scalar variance scale for the positive parameter `nu`; `c.n` is used in the prior variance `1/c.n`. |
| `Un` | Upper truncation bound for `nu`; the sampler draws `nu` on `[0, Un]`. |
| `Ls`, `Us` | Lower and upper truncation bounds for `sigma^2`. |
| `simulate` | Total number of Gibbs sampling iterations. |
| `burnin` | Number of initial iterations discarded as burn-in. |
| `Break` | Interval between saved posterior draws after burn-in. |

The implementation requires compatible vector/matrix dimensions, a full-rank `F`, positive `c.b` and `c.n`, and sensible positive truncation bounds for `sigma^2`. The response support is expected to be consecutive; the current sampler has no update branch for `Y < 3`.

## Example: GM12878 intra-chromosomal data

The dataset contains 1,000 GM12878 records. It is whitespace-delimited, has no header, and contains exactly nine columns:

1. Observed interaction frequency, `y_i`
2. Start of anchor A, `f_i^{x_A}`
3. End of anchor A, `f_i^{y_A}`
4. Start of anchor B, `f_i^{x_B}`
5. End of anchor B, `f_i^{y_B}`
6. Mappability, `f_i^{map}`
7. GC percentage, `f_i^{gc}`
8. Marginal count `mc(A_i)`
9. Marginal count `mc(B_i)`

**Repository data note:** `GM.txt` is not present in the current GitHub `main` tree. Add/upload the dataset to the repository root first. The raw URL below works only after the file has been uploaded there:

```r
url <- "https://raw.githubusercontent.com/zareifard/BLIM-CLVM/main/GM.txt"
tmp <- tempfile(fileext = ".txt")
download.file(url, destfile = tmp, mode = "wb")
D <- read.table(tmp, header = FALSE)
```

Construct the response, design matrix, and marginal counts:

```r
Y <- D[[1]]
F <- cbind(rep(1, length(Y)), as.matrix(D[, 2:7]))
mc.A <- D[[8]]
mc.B <- D[[9]]
beta <- as.vector(solve(t(F) %*% F) %*% t(F) %*% Y)
```

Run CLVM using the example settings:

```r
fit <- CLVM(
  Y = Y, F = F, mc.A = mc.A, mc.B = mc.B,
  L = 0.3, U = 0.7,
  m.b = beta, c.b = 10^4,
  m.n = 1, c.n = 10, Un = 10,
  Ls = 0.01, Us = 3 * var(Y),
  simulate = 100, burnin = 10, Break = 2
)
```

## Returned output

The returned list includes posterior draws for `beta`, `nu`, `sig2`, `X`, and `Z`, as well as summaries `E.beta`, `V.beta`, `E.nu`, `V.nu`, `HPD.nu`, `E.sig2`, `V.sig2`, `HPD.sig2`, `E.X`, and `E.Z`. HPD intervals are calculated at the 0.05 level by `boa.hpd`.

## Reproducibility

Set and record an R random seed with `set.seed(...)` before running the sampler if reproducible simulation output is required. Record the input data, package versions, and MCMC settings with the analysis. The data-download example depends on `GM.txt` first being added to the repository root.
