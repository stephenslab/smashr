This repository contains an R package for performing Gaussian denoising (also known as nonparametric regression) with heteroskedastic errors, as well as Poisson denoising.

To install the package run the following lines
```
install.packages("devtools")
library(devtools)
install_github("zrxing/smashr")
```

To use the package run 
```
library(smash)
```
and for help using the package run
```
?smash
```

#Background
Nonparametric regression is often used when one does not wish to assume any particular relationship (eg. linear) between a predictor and a response. Commonly seen examples come from astronomy, finance and biology. 

For the Gaussian case, the statistical problem is of the form Y_i=\mu_i+\epsilon_i, i=1,...,n, where \mu is the underlying mean function and assumed to be "spatially structured", and \epsilon_i's are independent Gaussian noise with mean 0 and variance \sigma_i^2. The goal is to recover \mu as accurately as possible given Y.

For the Poisson case, the model is similar, given by Y_i\sim Pois(\mu_i), i=1,...,n, where \mu is the underlying intensity for the Poisson process, and also assumed to be "spatially structured". The goal is to recover \mu as accurately as possible given Y.

The key idea behind the method is via wavelet transform (eg. Donoho & Johnstone (1994)) for the Gaussian case, and an extremely similar multiscale decomposition (Kolaczyk (1999)) for the Poisson case. Instead of traditional ways to perform shrinkage/thresholding however, both versions of smash make use of the Adaptive Shrinkage procedure introduced by Stephens (2015). Other details can be found in the companion paper Xing and Stephens (in preparation).

#Demonstration
The method can be used straight out of the box using the default settings, and should work for most reasonable scenarios (eg the noise does not completely overwhelm the signal). An example is shown below:

```
# Create the baseline mean function (The "spikes" function is used as an example here)
n = 2^9
t = 1:n/n
spike.f = function(x) (0.75 * exp(-500 * (x - 0.23)^2) + 1.5 * exp(-2000 * (x - 0.33)^2) + 3 * exp(-8000 * (x - 0.47)^2) + 
    2.25 * exp(-16000 * (x - 0.69)^2) + 0.5 * exp(-32000 * (x - 0.83)^2))
mu.s = spike.f(t)

# Gaussian case
# Scale the signal to be between 0.2 and 0.8
mu.t = (1 + mu.s)/5
plot(mu.t, type = "l")
# Create the baseline variance function (The function V2 from Cai & Wang (2008) is used here)
var.fn = (1e-04 + 4 * (exp(-550 * (t - 0.2)^2) + exp(-200 * (t - 0.5)^2) + exp(-950 * (t - 0.8)^2)))/1.35
plot(var.fn, type = "l")
# Set the signal-to-noise ratio
rsnr = sqrt(5)
sigma.t = sqrt(var.fn)/mean(sqrt(var.fn)) * sd(mu.t)/rsnr^2
# Simulate an example dataset
X.s = rnorm(n, mu.t, sigma.t)
# Run smash
mu.est <- smash(X.s, "gaus")
# Plot the true mean function as well as the estimated one
plot(mu.t, type = "l")
lines(mu.est$mu.est, col = 2)

# Poisson case
# Scale the signal to be non-zero and to have a low average intensity
mu.t = 0.01 + mu.s
# Simulate an example dataset
X.s = rpois(n, mu.t)
# Run smash
mu.est = smash(X.s, "poiss")
# Plot the true mean function as well as the estimated one
plot(mu.t, type = "l")
lines(mu.est, col = 2) 
```



