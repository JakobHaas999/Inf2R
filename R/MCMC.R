## setup
library("checkmate")

# --- Markov Chains ---------------------------------------------------------------
# Markov Chain about the weather on the next day

# Matrix P:
P.weather <- matrix(c(
  1/3, 1/2, 1/6,
  1/6, 2/3, 1/6,
  1/6, 1/2, 1/3
), nrow = 3, byrow = TRUE)

# Labels for P
labels <- c("S", "C", "R")
rownames(P.weather) <- labels
colnames(P.weather) <- labels

# Function to do MCMC
MCMC <- function(n, P, start, seed = NULL) {
  # Input checks
  assertInt(n, lower = 1)
  assertMatrix(P, nrows = ncol(P), mode = "numeric", any.missing = FALSE)
  assertTRUE(
    all.equal(unname(apply(P, 1, sum)), rep(1, nrow(P)))
  )
  if (!is.null(names(P))) {
    assertTRUE(all(rownames(P) == colnames(P)))
  }
  assert(checkString(start),
         checkInt(start, lower = 1, upper = nrow(P)))
  if (is.character(start)) {
    assertTRUE(start %in% rownames(P))
  }
  assertInt(seed, null.ok = TRUE)
  
  chain <- character(n)
  chain[[1]] <- start
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  for(i in seq_len(n)[-1]) {
    prior <- chain[[i - 1]]
    chain[[i]] <- sample(colnames(P), size = 1, prob = P.weather[prior, ])
  }
  prop.table(table(chain))
}

# Simulation of our Markov Chain
# Parameter
n <- 1000
start.weather <- "C"
seed <- 1922

MCMC(n, P = P.weather, start = start.weather, seed = seed)[labels]

# --- Metropolis-Hastings-algorithm ------------------------------------------------

# Implement the Metropolis-Hastings-algorithm for Y = {x in X | sum x_i > 10}
# with X_i ~ Exp(1)

## Helper functions for the algorithm

# Target density
# Use a function to evaluate the logarithm of the target density
target <- function(x) {
  if (sum(x) > 10) {
    # use log=TRUE to avoid numerical stability issues for low probabilities
    # sum because log space, otherwise we would multiply
    sum(dexp(x, rate = 1, log = TRUE))
  } else {
    -Inf
  }
}

# Proposal density
proposal <- function(x, y, sigma = 1) {
  # Again use log = TRUE
  sum(dnorm(y, x, sd = sigma, log = TRUE))
}

# acceptance probability
acceptance <- function(x, y) {
  min(1,
      exp(
        (target(y) - target(x)) +
         proposal(x, y) - proposal(y, x)
      ))
}

# Implement the algorithm
metropolis_hastings <- function(target, proposal, acceptance,
                                x.0, n.samples = 100000, sigma = 1) {
  res <- matrix(NA, nrow = n.samples, ncol = length(x.0))
  x.k <- x.0
  accepted <- 0
  for (i in seq_len(n.samples)) {
    ystar <- rnorm(length(x.k), mean = x.k, sd = sigma)
    alpha <- acceptance(x.k, ystar)
    if (alpha == 1 || runif(1) < alpha) {
      x.k <- ystar
      accepted <- accepted + 1
    }
    res[i, ] <- x.k
  }
  cat("Acceptance rate was: ", accepted / n.samples)
  res
}

# Perform the algorithm
samples<- metropolis_hastings(target, proposal, acceptance, x.0 = c(5.01, 5.01))

# Plot results

burn.in <- 1000  # burn in time

par(mfrow = c(1, 3))
plot(samples[burn.in:2500, 1], samples[burn.in:2500, 2], type = "l",
     main = "Plot of the sample",
     xlab = "X1", ylab = "X2")

plot(rowSums(samples)[0:10000], type = "l",
     main = "Sum S of the two variables over time",
     xlab = "Iteration", ylab = "S")

plot(density(rowSums(samples[(burn.in+1):nrow(samples), ])),
     main = paste("Distribution of the sum of the", 2,
                  "variables"),
     xlim = c(9, 17))

# Autocorrelation for the sum of the both variables
par(mfrow = c(1, 1))
acf(rowSums(samples), main = "Autocorrelation for the sum of the 2 variables")
# Thin samples to make the sample i.i.d
samples.iid <- samples[seq_len(nrow(samples)) %% 50 == 0, ]
acf(rowSums(samples.iid), main = "Autocorrelation for the sum of the 2 variables")

# --- Gibbs-Metropolis-Hastings ---------------------------------------------------------------

### Define new helper functions for the new algorithm

# New target
target_gibbs <- function(x.j, x) {
  if (sum(x) > 10) {
    # again log = TRUE
    dexp(x.j, rate = 1, log = TRUE)
  } else {
    -Inf
  }
}

# Univariate proposal distribution
proposal_gibbs <- function(x.j, y.j, sigma = 1) {
  dnorm(y.j, mean = x.j, sd = sigma, log = TRUE)
}

acceptance_gibbs <- function(x.j, y.j, x.new, x.current) {
  log.alpha <- target_gibbs(x.j, x.new) - target_gibbs(y.j, x.current) +
    proposal_gibbs(y.j, x.j) - proposal_gibbs(x.j, y.j)
  min(1, exp(log.alpha))
}


# Implement gibbs sampling
gibbs_sampling <- function(target, proposal, acceptance,
                           x.0, n.samples = 100000, sigma = 1) {
  res <- matrix(NA, nrow = n.samples, ncol = length(x.0))
  accepted <- 0
  x.k <- x.0
  
  for (i in seq_len(n.samples)) {
    for (j in seq_along(x.k)) {
      proposed <- rnorm(1, mean = x.k[[j]], sd = sigma)
      proposed.x <- x.k
      proposed.x[[j]] <- proposed
      alpha <- acceptance(proposed, x.k[[j]], proposed.x, x.k)
      if (alpha == 1 || runif(1) < alpha) {
        x.k[[j]] <- proposed
        accepted <- accepted + 1
      }
    }
    res[i, ] <- x.k
  }
  cat("Acceptance rate was:", accepted / (n.samples * length(x.0)), "\n")
  res
}

# Perform the algorithm
samples.gibbs <- gibbs_sampling(target_gibbs, proposal_gibbs, acceptance_gibbs,
                                c(5.01, 5.01))

# Plot the result
plot(samples.gibbs[burn.in:(burn.in + 2500), 1], samples.gibbs[burn.in:(burn.in + 2500), 2],
     type = "l", xlab = "X1", ylab = "X2")


# --- Metropolis-Hastings-algorithm with a discrete random variable ------------------------------
# TODO