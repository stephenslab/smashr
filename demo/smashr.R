# Illustration of Gaussian and Poisson denoising using the "smash"
# function.
library(smashr)

# Create the baseline mean function. Here we use the "spikes" function.
n = 2^9
t = 1:n/n
spike.f = function(x) (0.75 * exp(-500 * (x - 0.23)^2) +
            1.5  * exp(-2000 * (x - 0.33)^2) +
	    3    * exp(-8000 * (x - 0.47)^2) + 
            2.25 * exp(-16000 * (x - 0.69)^2) +
	    0.5  * exp(-32000 * (x - 0.83)^2))
mu.s = spike.f(t)

# GAUSSIAN CASE
# -------------
# Scale the signal to be between 0.2 and 0.8.
mu.t = (1 + mu.s)/5
plot(mu.t, type = "l")

# Create the baseline variance function. (The function V2 from Cai &
# Wang 2008 is used here.)
var.fn = (1e-04 + 4 * (exp(-550 * (t - 0.2)^2) +
  exp(-200 * (t - 0.5)^2) + exp(-950 * (t - 0.8)^2)))/1.35
plot(var.fn, type = "l")

# Set the signal-to-noise ratio.
rsnr = sqrt(5)
sigma.t = sqrt(var.fn)/mean(sqrt(var.fn)) * sd(mu.t)/rsnr^2

# Simulate an example dataset.
X.s = rnorm(n, mu.t, sigma.t)

# Run smash on these data.
mu.est <- smash(X.s, "gaus")

# Plot the estimated mean function (red) against the ground-truth (black).
plot(mu.t, type = "l")
lines(mu.est, col = 2)

# POISSON CASE
# ------------
# Scale the signal to be non-zero and to have a low average intensity.
mu.t = 0.01 + mu.s

# Simulate an example dataset.
X.s = rpois(n, mu.t)

# Run smash.
mu.est = smash(X.s, "poiss")

# Plot the estimated mean function (red) against the ground-truth (black).
plot(mu.t, type = "l")
lines(mu.est, col = 2) 
