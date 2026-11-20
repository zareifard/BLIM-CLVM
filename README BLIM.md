# Binary Latent Indicator Model (BLIM)

BLIM is an R implementation of a Bayesian model for interaction frequencies in ChIA-PET data. The model can be applied to both intra- and inter-chromosomal interactions; the example analysis in this README focuses on **intra-chromosomal interactions**.

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
source("BLIM.R")
```

## Function

```r
BLIM <- function(Y,F,mc.A,mc.B,L,U,v.p0=0.2,m.b,c.b,m.bz,c.bz,Up.z,Ls,Us,simulate,burnin,Break)
```

Only `v.p0` has a default; all other arguments are required.

| Argument | Description |
|---|---|
| `Y` | Vector of observed interaction frequencies, `y_i`. |
| `F` | Design matrix of covariates, including an intercept column. |
| `mc.A`, `mc.B` | Marginal PET counts for anchors A and B. |
| `L`, `U` | Lower and upper bounds for the prior mean of the mixture weights. |
| `v.p0` | Prior variance used to determine the Beta prior parameters for the mixture weights. Defaults to `0.2`. |
| `m.b`, `c.b` | Prior mean and variance parameter for `beta`. |
| `m.bz`, `c.bz` | Prior mean and variance parameter for `alpha`. |
| `Up.z` | Upper truncation point for `alpha`. |
| `Ls`, `Us` | Lower and upper truncation bounds for `sigma^2`. |
| `simulate` | Total number of MCMC iterations. |
| `burnin` | Number of initial iterations discarded as burn-in. |
| `Break` | Number of iterations between successive saved posterior samples. |

## Example: GM12878 intra-chromosomal data

The example dataset contains 1,000 GM12878 records. It is whitespace-delimited, has no header, and contains exactly nine columns:

1. Observed interaction frequency, `y_i`
2. Start of anchor A, `f_i^{x_A}`
3. End of anchor A, `f_i^{y_A}`
4. Start of anchor B, `f_i^{x_B}`
5. End of anchor B, `f_i^{y_B}`
6. Mappability score, `f_i^{map}`
7. GC content percentage, `f_i^{gc}`
8. `mc(A_i)`
9. `mc(B_i)`

**Repository data note:** `GM.txt` is not yet included in the current GitHub `main` tree. Add/upload `GM.txt` to the repository root before running the download code below. The raw-file URL will work only after that file has been added:

```r
url <- "https://raw.githubusercontent.com/zareifard/BLIM-CLVM/main/GM.txt"
tmp <- tempfile(fileext = ".txt")
download.file(url, destfile = tmp, mode = "wb")
D <- read.table(tmp, header = FALSE)
```

Construct the response, design matrix, anchor counts, and starting prior mean as follows:

```r
Y <- D[[1]]
n <- length(Y)
F <- cbind(rep(1,n), as.matrix(D[,2:7]))
mc.A <- D[[8]]
mc.B <- D[[9]]
beta <- as.vector(solve(t(F)%*%F)%*%t(F)%*%Y)
```

Run BLIM with the example settings:

```r
fit <- BLIM(
  Y = Y, F = F, mc.A = mc.A, mc.B = mc.B,
  L = 0.3, U = 0.7, v.p0 = 0.2,
  m.b = beta, c.b = 10^4,
  m.bz = 2, c.bz = 10,
  Up.z = 100, Ls = 0.01, Us = 3*var(Y),
  simulate = 100, burnin = 10, Break = 2
)
```

## Returned output

The returned object contains posterior samples for `beta`, `alpha`, `sig2`, `X`, and `Z`, along with summaries including `E.beta`, `V.beta`, `E.beta.z`, `V.beta.z`, `HPD.beta.z`, `E.sig2`, `V.sig2`, `HPD.sig2`, `E.X`, and `E.Z`. These fields summarize the corresponding posterior samples; the reported HPD calls use a 0.05 level.

## Reproducibility

For reproducible simulation output, set and record an R random seed with `set.seed(...)` before calling `BLIM`. Keep the input data, package versions, and MCMC settings together with the analysis. The example download requires `GM.txt` to have first been uploaded to the repository root as noted above.
