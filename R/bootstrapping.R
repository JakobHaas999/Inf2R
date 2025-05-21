## setup
library("checkmate")
library("latex2exp")

# --- Bootstrapping ---------------------------------------------------------------

# Sample X_1, ..., X_120 with X_i ~ 0.6 N(0, 1) + 0.4 N(5, 2^2)

# Function that samples from this mixture
sample_from_mixture <- function(n, seed = NULL) {
  assertInt(n, lower = 1)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  replicate(n, {
    u <- runif(1)
    if (u <= 0.6) {
      rnorm(1, 0, 1)
    } else {
      rnorm(1, 5, 2)
    }
  })
}

# Perform the function on our sample
seed <- 123  # seed
n.samples <- 120 # sample size

sample.mixture <- sample_from_mixture(n = n.samples, seed = seed)

# Compute IQR
iqr.sample <- IQR(sample.mixture)
cat(sprintf("The IQR of our sample is %s", round(iqr.sample, 3)))

# Function for bootstrapping

bootstrap_from_sample <- function(smp, B, seed = NULL) {
  assertNumeric(smp, any.missing = FALSE)
  assertInt(B, lower = 1)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  replicate(B, {
    bootstrap <- smp[sample.int(length(smp), length(smp), replace = TRUE)]
    IQR(bootstrap)
  })
}

# Perform the function on "sample.mixture"
seed <- 1234  # seed
B <- 500  # bootstrap samples

iqr.bootstrap <- bootstrap_from_sample(sample.mixture, B = B, seed = seed)

# Compute SE and Skewness

se.iqr <- sd(iqr.bootstrap)
skew.iqr <- moments::skewness(iqr.bootstrap)

cat(sprintf("The SE of our IQR sample is %s",
            round(se.iqr, 3)))
cat(sprintf("The skewness of our IQR sample is %s",
            round(skew.iqr, 3)))

# Plot a histogram of "iqr.bootstrap" with overlaid kernel density estimate
hist(iqr.bootstrap, probability = TRUE, ylim = c(0, 1),
     xlab = TeX("$\\{I_b^*\\}$"))
lines(density(iqr.bootstrap), col = "red")

# Calculate Confidence Intervalls
CI.bootstrap <- quantile(iqr.bootstrap, probs = c(0.025, 0.975))
CI.normal.approx <- c(
  lower.bound = mean(iqr.bootstrap) - 1.96 * se.iqr,
  upper.bound = mean(iqr.bootstrap) + 1.96 * se.iqr
)

# Calculate widths
width.CI.boot <- CI.bootstrap[[2]] - CI.bootstrap[[1]]
width.CI.norm.approx <- CI.normal.approx[[2]] - CI.normal.approx[[1]]

cat("Width of the 95% bootstrap CI is: ", width.CI.boot)
cat("Width of the 95% normal-approximation CI is: ", width.CI.norm.approx)



