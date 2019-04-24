# Example 1: n is a power of 2
# ----------------------------
# Create the baseline mean function. (The "spikes" function is used
# as an example here.)
n <- 2^9
t <- 1:n/n
spike.f <- function (x) (0.75 * exp(-500 * (x - 0.23)^2) +
  1.5 * exp(-2000 * (x - 0.33)^2) + 3 * exp(-8000 * (x - 0.47)^2) +
  2.25 * exp(-16000 * (x - 0.69)^2) + 0.5 * exp(-32000 * (x - 0.83)^2))
mu.s <- spike.f(t)

# Gaussian case
# -------------
# Scale the signal to be between 0.2 and 0.8
mu.t <- (1 + mu.s)/5

# Create the baseline variance function. (The function V2 from Cai &
# Wang (2008) is used here.)
var.fn <- (1e-04 + 4 * (exp(-550 * (t - 0.2)^2) +
                        exp(-200 * (t - 0.5)^2) +
                        exp(-950 * (t - 0.8)^2)))/1.35

# Set the signal-to-noise ratio.
rsnr    <- sqrt(5)
sigma.t <- sqrt(var.fn)/mean(sqrt(var.fn)) * sd(mu.t)/rsnr^2

# Simulate an example dataset.
X.s <- rnorm(n,mu.t,sigma.t)

# Run smash (Gaussian version is run since observations are not
# counts).
print(system.time(
  mu.est <- smash(X.s,ashparam = list(control = list(verbose = TRUE)))))

# Plot the true mean function as well as the estimated one.
plot(mu.t,type = "l")
lines(mu.est,col = 2)

# Example 2: n is not a power of 2
# --------------------------------
n    <- 1000
t    <- 1:n/n
mu.s <- spike.f(t)

# Scale the signal to be between 0.2 and 0.8.
mu.t <- (1 + mu.s)/5

# Create the baseline variance function.
var.fn  <- (1e-04 + 4 * (exp(-550 * (t - 0.2)^2) +
                         exp(-200 * (t - 0.5)^2) +
                         exp(-950 * (t - 0.8)^2)))/1.35
sigma.t <- sqrt(var.fn)/mean(sqrt(var.fn)) * sd(mu.t)/rsnr^2

# Simulate an example dataset.
X.s <- rnorm(n,mu.t,sigma.t)

# Run smash.
print(system.time(
  mu.est <- smash(X.s,ashparam = list(control = list(verbose = TRUE)))))

# Plot the true mean function as well as the estimated one.
plot(mu.t,type = "l")
lines(mu.est,col = 2)

print(sessionInfo())
