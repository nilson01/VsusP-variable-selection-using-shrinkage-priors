pkgname <- "VsusP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "VsusP-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('VsusP')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("OptimalHbi")
### * OptimalHbi

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: OptimalHbi
### Title: Variable selection using shrinkage priors :: OptimalHbi
### Aliases: OptimalHbi

### ** Examples


n <- 100
p <- 20
X <- matrix(rnorm(n * p), n, p)
beta <- exp(rnorm(p))
Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
df <- data.frame(X, Y)
rv.hs <- bayesreg::bayesreg(Y ~ ., df, "gaussian", "horseshoe+", 200, 100)

Beta <- t(rv.hs$beta)
lower <- 0
upper <- 1
l <- 5
S2Mbeta <- Sequential2MeansBeta(Beta, lower, upper, l)

bi <- S2Mbeta$b.i
Hbi <- S2Mbeta$H.b.i
OptimalHbi(bi, Hbi)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("OptimalHbi", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("S2MVarSelection")
### * S2MVarSelection

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: S2MVarSelection
### Title: Variable selection using shrinkage priors :: S2MVarSelection
### Aliases: S2MVarSelection

### ** Examples


n <- 100
p <- 20
X <- matrix(rnorm(n * p), n, p)
beta <- exp(rnorm(p))
Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
df <- data.frame(X, Y)
# Fit a model using gaussian horseshoe+ for 200 samples
# # recommended n.samples is 5000 and burning is 2000
rv.hs <- bayesreg::bayesreg(Y ~ ., df, "gaussian", "horseshoe+", 200, 100)

Beta <- rv.hs$beta
H <- 12
impVariablesGLM <- S2MVarSelection(Beta, H)
impVariablesGLM




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("S2MVarSelection", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Sequential2Means")
### * Sequential2Means

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Sequential2Means
### Title: Variable selection using shrinkage priors :: Sequential2Means
### Aliases: Sequential2Means

### ** Examples

# -----------------------------------------------------------------
# Example 1: Gaussian Model and Horseshoe prior
n <- 100
p <- 20
X <- matrix(rnorm(n * p), n, p)
beta <- exp(rnorm(p))
Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
b.i <- seq(0, 1, 0.05)

# Sequential2Means with horseshoe+ using gibbs sampling
# recommended n.samples is 5000 and burning is 2000
S2M <- Sequential2Means(X, Y, b.i, "horseshoe+", 200, 100)
Beta <- S2M$Beta
H.b.i <- S2M$H.b.i

# -----------------------------------------------------------------
# Example 2: Gaussian Model and ridge prior

n <- 100
p <- 20
X <- matrix(rnorm(n * p), n, p)
beta <- exp(rnorm(p))
Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
b.i <- seq(0, 1, 0.05)
# Sequential2Means with ridge regression using gibbs sampling
# recommended n.samples is 5000 and burning is 2000
S2M <- Sequential2Means(X, Y, b.i, "ridge", 200, 100)
Beta <- S2M$Beta
H.b.i <- S2M$H.b.i




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Sequential2Means", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Sequential2MeansBeta")
### * Sequential2MeansBeta

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Sequential2MeansBeta
### Title: Variable selection using shrinkage prior :: Sequential2MeansBeta
### Aliases: Sequential2MeansBeta

### ** Examples


# -----------------------------------------------------------------
# Example 1: Gaussian Model and Horseshoe prior

n <- 100
p <- 20
X <- matrix(rnorm(n * p), n, p)
beta <- exp(rnorm(p))
Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
df <- data.frame(X, Y)

# beta samples for gaussian model using horseshow prior and gibbs sampling
rv.hs <- bayesreg::bayesreg(Y ~ ., df, "gaussian", "horseshoe+", 200, 100)

Beta <- t(rv.hs$beta)
lower <- 0
upper <- 1
l <- 20
S2Mbeta <- Sequential2MeansBeta(Beta, lower, upper, l)
H.b.i <- S2Mbeta$H.b.i

# -----------------------------------------------------------------
# Example 2: normal model and lasso prior

#' n <- 100
p <- 20
X <- matrix(rnorm(n * p), n, p)
beta <- exp(rnorm(p))
Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
df <- data.frame(X, Y)
rv.hs <- bayesreg::bayesreg(Y ~ ., df, "normal", "lasso", 200, 100)

Beta <- t(rv.hs$beta)
lower <- 0
upper <- 1
l <- 15
S2Mbeta <- Sequential2MeansBeta(Beta, lower, upper, l)
H.b.i <- S2Mbeta$H.b.i




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Sequential2MeansBeta", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
