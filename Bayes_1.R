

# load data
load("Data/Ex_3_2.RData")
x_t <- data$x_t
pi_it <- data$pi_it
pi_new <- data$pi_new


# Prior props of the four models
prior.probs <- c(M1 = 0.3, M2 = 0.35, M3 = 0.2, M4 = 0.15)

# Calculate posterior props
post.probs.unnormed <- numeric(length(prior.probs))
names(post.probs.unnormed) <- names(prior.probs)

for (i in seq_along(post.probs.unnormed)) {
  likelihood <- prod(
    (pi_it[i, ]^x_t) * (1 - pi_it[i, ])^(1 - x_t)
  )
  post.probs.unnormed[[i]] <- likelihood * prior.probs[[i]]
}

# Norm post
post.probs <- post.probs.unnormed / sum(post.probs.unnormed)

# Get the most plausibel model

most.plausibel.model <- post.probs[which.max(post.probs.unnormed)]
cat(sprintf("The most plausible model is '%s' with a posterior probability of %s",
            names(most.plausibel.model), round(most.plausibel.model[[1]], 3)), "\n")

# Predictive posterior propability
pi.hat <- sum(pi_new * post.probs)
cat("Predictive posterior propability is: ", round(pi.hat, 3), "\n")



## Now failures occured on 20 days

# Plotting distribution of each model’s predictions
par(mfrow = c(2, 2))
for (i in 1:4) {
  hist(pi_it[i, ], main = sprintf("Model: %s", i), xlab = "")
}
# -> Model 1: low probality for failures
# -> Model 2: neutral
# -> Model 3: low values in the the intervall [0.2, 0.4]
# -> Model 4: high probality for failures

# Calculate everything for the new circumstances
set.seed(1111)
x_t_new <- sample(c(rep(1, 20), rep(0, 10)))

post.probs.new <- numeric(length(prior.probs))
names(post.probs.new) <- names(prior.probs)

for (i in seq_along(post.probs.new)) {
  likelihood <- prod(
    (pi_it[i, ]^x_t_new) * (1 - pi_it[i, ])^(1 - x_t_new)
  )
  post.probs.new[[i]] <- likelihood * prior.probs[[i]]
}
# Norm probs
post.probs.new <- post.probs.new / sum(post.probs.new)

most.plausibel.model.new <- post.probs.new[which.max(post.probs.new)]
pi.hat.new <- sum(post.probs.new * pi_new)

cat(sprintf("The most plausible model is now '%s' with probality %s",
            names(most.plausibel.model.new), round(most.plausibel.model.new[[1]], 3)), "\n")
cat("New predictive posterior propability is: ", round(pi.hat.new, 3))



target <- function(x, mu) {
  dnorm(x, mu, 1)
}

proposal <- function(x, mu_t, sigma2) {
  rnorm(x, mu_t, sqrt(sigma2))
}

mu_0 <- 1
n_sim <- 1000
res <- numeric(n_sim)
sigma2 <- 3