## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>")

## ----eval = FALSE-------------------------------------------------------------
#  # Plain installation
#  devtools::install_github("nilson01/VsusP")
#  
#  # For installation with vignette
#  devtools::install_github("nilson01/VsusP", build_vignettes = TRUE)
#  
#  # load the library
#  library(VsusP)
#  

## ----example-setup------------------------------------------------------------
# Assuming MASS is installed, if not uncomment the next line
if (!requireNamespace("MASS", quietly = TRUE)) {
    install.packages("MASS")
}
library(MASS)


# Set seed for reproducibility
set.seed(20221208)

# Simulate data
sim.XY <- function(n, p, beta) {
  X <- matrix(rnorm(n * p), n, p)
  Y <- X %*% beta + rnorm(n)
  return(list(X = X, Y = Y, beta = beta))
}

n <- 100
p <- 20
beta <- exp(rnorm(p))
data <- sim.XY(n, p, beta)

## ----s2m-apply----------------------------------------------------------------

b.i <- seq(0, 1, by = 0.05)
S2M <- VsusP::Sequential2Means(X = data$X, Y = as.vector(data$Y), b.i = b.i, prior = "horseshoe+", n.samples = 5000, burnin = 2000)
Beta <- S2M$Beta
H.b.i <- S2M$H.b.i

## ----s2m-results--------------------------------------------------------------
VsusP::OptimalHbi(bi = b.i, Hbi = H.b.i)

optimal_b_i <- b.i[which.min(S2M$H.b.i)]
optimal_Hbi <- min(S2M$H.b.i)

cat("Optimal b.i: \n", optimal_b_i, "\n")
cat("Optimal H.b.i: \n", optimal_Hbi, "\n")


## ----s2m-results important Variables------------------------------------------
H <- 13
# Variable selection
impVariablesGLM <- VsusP::S2MVarSelection(Beta, H)
impVariablesGLM

## ----sequential2means-example-------------------------------------------------
set.seed(20221208)
n <- 100
p <- 20
X <- matrix(rnorm(n * p), n, p)
Y <- X %*% exp(rnorm(p)) + rnorm(n)
b.i <- seq(0, 1, length.out = 20)
result <- VsusP::Sequential2Means(X, as.vector(Y), b.i, prior = "horseshoe+", n.samples = 5000, burnin = 2000)

cat("Beta: \n")
print(result$Beta[1:5, ])
cat("\n\n")

cat("H.b.i: \n")
print(result$H.b.i)


## ----optimal-hbi-example------------------------------------------------------
VsusP::OptimalHbi(b.i, result$H.b.i)

## ----s2mvar-selection-example-------------------------------------------------
H <- 17 
significantVariables <- VsusP::S2MVarSelection(result$Beta, H)
cat("Significant Variables: \n")
print(significantVariables)


## ----sequential2meansbeta-example---------------------------------------------
Beta <- result$Beta
lower <- 0
upper <- 1
l <- 20

S2MBeta = VsusP::Sequential2MeansBeta(Beta, lower, upper, l)

cat("p : \n")
print(S2MBeta$p)
cat("\n\n")


cat("bi : \n")
print(S2MBeta$b.i)
cat("\n\n")


cat("Hbi : \n")
print(S2MBeta$H.b.i)

## ----Simulation---------------------------------------------------------------

library(MASS)
library(VsusP)

set.seed(12345) 

n <- 50
p <- 300
r <- 10
rho <- 0.95

Sigma <- diag(1, p)
block_size <- 5
num_blocks <- p / block_size

for (b in 0:(num_blocks - 1)) {
  block_start <- b * block_size + 1
  block_end <- min(p, block_start + block_size - 1)
  Sigma[block_start:block_end, block_start:block_end] <- rho
  diag(Sigma[block_start:block_end, block_start:block_end]) <- 1
}

Sigma <- Sigma + diag(0.001, p)  # Ensure positive definiteness

X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
beta_true <- c(6, rep(3, 4), rep(1, 5), rep(0, p - 10))
Y <- as.vector(X %*% beta_true + rnorm(n))

var_y <- var(Y)
b_i_range <- seq(0.5 * var_y, 10 * var_y, length.out = 20)  # Check the scale of var(Y)

results <- Sequential2Means(X = X, Y = Y, b.i = b_i_range, prior = "horseshoe+", n.samples = 5000, burnin = 2000)

plot(b_i_range, results$H.b.i, type = 'b', xlab = "Tuning parameter b.i", ylab = "Estimated number of signals (H.b.i)",
     main = "Plot of b.i vs H.b.i for Simulation 2")

optimal_b_i <- b_i_range[which.min(results$H.b.i)]
optimal_Hbi <- min(results$H.b.i)

cat("Optimal b.i:", optimal_b_i, "\n")
cat("Optimal H.b.i:", optimal_Hbi, "\n")

H <- optimal_Hbi
# Variable selection
Beta <- results$Beta

impVariablesGLM <- S2MVarSelection(Beta, H)
impVariablesGLM



