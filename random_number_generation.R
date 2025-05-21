## setup
library("checkmate")

# --- Random Variable Generation ---------------------------------------------------------------

# Implementation of RNG with 42 as starting value
lehmerRNG <- function(x.0, a = 7^5, m = 2^31 - 1, n = 1000) {
  assertInt(x.0, lower = 1)
  assertInt(a)
  assertInt(m, lower = 1)
  assertInt(n, lower = 1)
  X <- integer(n + 1)
  U <- numeric(n)
  X[[1]] <- x.0
  for (i in seq_len(n)) {
    X[[i + 1]] <- (a * X[[i]]) %% m
    U[[i]] <- X[[i + 1]] / m
  }
  U
}

# Test the function with start value 42
resRng <- lehmerRNG(42)
par(mfrow = c(1, 2))
plot(resRng)
F.resRng <- ecdf(resRng)
plot(F.resRng, main = "")

## Inverse-transformation
# X ~ exp(lambda = 2)
# Our inverse cdf is:
f.inverse.exp <- function(u, lambda = 2) {
  assertNumber(u, lower = 0, upper = 1)
  assertCount(lambda)
  - log(1 - u) / lambda
}

X.exp <- vapply(resRng, f.inverse.exp, numeric(1))

# Histogramm of result
par(mfrow = c(1, 1))
hist(X.exp,
     breaks = 30,
     probability = TRUE,
     xlab = "x",
     main = ""
     )
curve(dexp(x, rate = 2), add = TRUE, lwd = 2, col = "red")

# 1000 sims of Kolmogorov-Smirnov test: sample against Exp(2)

sims <- 1000
p.values.ks <- numeric(sims)

for (i in seq_len(sims)) {
  U <- lehmerRNG(42)
  X.exp.sim <- vapply(U, f.inverse.exp, numeric(1))
  test <- ks.test(X.exp.sim, "pexp", rate = 2)
  p.values.ks[[i]] <- test$p.value
}

round(mean(p.values.ks), 4)

# 1000 sims of chi-square test: sample against Exp(2)

p.values.chisq <- numeric(sims)

for (i in seq_len(sims)) {
  U <- lehmerRNG(42)
  X.exp.sim <- vapply(U, f.inverse.exp, numeric(1))
  brks <- seq(0, max(X.exp.sim), l = 11)
  freqs <- hist(X.exp.sim, breaks = brks, plot = FALSE)$counts
  # expexted probalities
  p <- numeric(length(freqs))
  for (j in seq_len(length(p))) {
    p[[j]] <- pexp(brks[[j + 1]], rate = 2) - pexp(brks[[j]], rate = 2)
  }
  test <- suppressWarnings(chisq.test(x = freqs, p = p, rescale.p = TRUE))
  p.values.chisq[[i]] <- test$p.value
}

round(mean(p.values.chisq), 4)

# --- Rejection Sampling ---------------------------------------------------------------

# Function that does rejection sampling
rejectionSampler <- function(n.samples, a, verbose = TRUE) {
  assertInt(n.samples, lower = 1)
  assertNumber(a, lower = 0)
  stopifnot(a != 0)
  assertFlag(verbose)
  
  samples <- numeric(n.samples)
  max.dens <- exp(0)
  M <- max.dens * (2 * a)
  cat("M: ", M, "\n")
  iters <- 0
  i <- 0
  
  while(i < n.samples) {
    iters <- iters + 1
    y <- runif(1, -a, a)
    u <- runif(1)
    target.dens <- exp(-y^2)
    prob.dens <- 1 / (2 * a)
    
    if (u <= target.dens /( M * prob.dens)) {
      i <- i + 1
      samples[[i]] <- y
    }
  }
  
  if (verbose) {
    cat("Acceptance rate:", n.samples/iters, "\n")
  }
  return(samples)
}

# Generate samples
set.seed(1234)
a.values <- c(0.1, 1, 10)
n <- 10000
results.rej.samp <- list()

for (a in a.values) {
  samp <- rejectionSampler(n.samples = n,
                            a         = a,
                            verbose   = FALSE)
  results.rej.samp[[paste0("a=", a)]] <- samp
}

# Plot results
par(mfrow = c(1, 3))
plot(density(results.rej.samp$`a=0.1`), main = "Density of Samples (a = 0.1)",
     xlab = "x", ylab = "Density", xlim = c(-2.5, 2.5))
plot(density(results.rej.samp$`a=1`), main = "Density of Samples (a = 1)",
     xlab = "x", ylab = "Density", xlim = c(-2.5, 2.5))
plot(density(results.rej.samp$`a=10`), main = "Density of Samples (a = 10)",
     xlab = "x", ylab = "Density", xlim = c(-2.5, 2.5))

# --- Importance Sampling -------------------------------------------------------------------

# We use importance sampling now to estimate the integral from -a to a of
# h(x) = cos(x) / (1 + x^2)

a <- 1
samples <- rejectionSampler(n, a)

# Function h
h <- function(x) {
  cos(x) / (1 + x^2)
}

# Plot h
par(mfrow = c(1, 1))
curve(h, from = -a, to = a, col = "blue", lwd = 2,
      ylab = "h(x)", xlab = "x", main = "Function to Estimate")
# next to the target densit
curve(exp(-x^2), from = -a, to = a, col = "red", lwd = 2,
      add = TRUE)
legend("bottomright", legend = c("h(x)", "exp(-x^2)"),
       col = c("blue", "red"), lty = c(1, 1), lwd = c(2, 2))

f.unnorm <- function(x) {
  exp(-x^2)
}
norm.constant <- integrate(f.unnorm, -a, a)$value
f.normed <- function(x) {
  f.unnorm(x) / norm.constant
}

# Look at scaling
curve(f.unnorm, from = -a, to = a, col = "red", lwd = 2,
      ylab = "f(x)", xlab = "x", main = "Target Density",
      ylim = c(0, 1.2))
curve(f.normed, from = -a, to = a, col = "blue", lwd = 2,
      add = TRUE)

# Compute weights, for a
weights <- h(samples) / f.normed(samples)

# Importance sampling estimate of the integral, with SE
estimate <- mean(weights)
std_error <- sd(weights) / sqrt(n) # standard error of the mean

# Display result
cat("Estimate of the integral:", estimate, "\n")
cat("Estimated standard error:", std_error, "\n")
true_I <- integrate(h, lower = -a, upper = a)
cat("True value of the integral:", true_I$value, "\n")