#' @useDynLib smashr
#' 
#' @title smashr: Smoothing using Adaptive SHrinkage in R
#' 
#' @description This package performs nonparametric regression on
#' univariate Poisson or Gaussian data using multi-scale methods.  For
#' the Poisson case, the data $x$ is a vector, with $x_j \sim
#' Poi(\mu_j)$ where the mean vector $\mu$ is to be estimated.  For
#' the Gaussian case, the data $x$ are a vector with $x_j \sim
#' N(\mu_j, \sigma^2_j)$. Where the mean vector $\mu$ and variance
#' vector $\sigma^2$ are to be estimated.  The primary assumption is
#' that $\mu$ is spatially structured, so $\mu_j - \mu_{j+1}$ will
#' often be small (that is, roughly, $\mu$ is smooth). Also $\sigma$
#' is spatially structured in the Gaussian case (or, optionally,
#' $\sigma$ is constant, not depending on $j$).
#' 
#' @details The function \code{\link{smash}} provides a minimal
#' interface to perform simple smoothing.  It is actually a wrapper to
#' \code{\link{smash.gaus}} and \code{\link{smash.poiss}} which
#' provide more options for advanced use.  The only required input is
#' a vector of length 2^J for some integer J.  Other options include
#' the possibility of returning the posterior variances, specifying a
#' wavelet basis (default is Haar, which performs well in general due
#' to the fact that smash uses the translation-invariant transform)
#' 
#' @author Matthew Stephens and Zhengrong Xing
#' 
#' @examples
#' # Create the baseline mean function (The "spikes" function is used as an example here)
#' n = 2^9
#' t = 1:n/n
#' spike.f = function(x) (0.75 * exp(-500 * (x - 0.23)^2) + 1.5 * exp(-2000 * (x - 0.33)^2) + 3 * exp(-8000 * (x - 0.47)^2) + 
#'     2.25 * exp(-16000 * (x - 0.69)^2) + 0.5 * exp(-32000 * (x - 0.83)^2))
#' mu.s = spike.f(t)
#' 
#' # Gaussian case
#' # Scale the signal to be between 0.2 and 0.8
#' mu.t = (1 + mu.s)/5
#' plot(mu.t, type = "l")
#' # Create the baseline variance function (The function V2 from Cai & Wang (2008) is used here)
#' var.fn = (1e-04 + 4 * (exp(-550 * (t - 0.2)^2) + exp(-200 * (t - 0.5)^2) + exp(-950 * (t - 0.8)^2)))/1.35
#' plot(var.fn, type = "l")
#' # Set the signal-to-noise ratio
#' rsnr = sqrt(5)
#' sigma.t = sqrt(var.fn)/mean(sqrt(var.fn)) * sd(mu.t)/rsnr^2
#' # Simulate an example dataset
#' X.s = rnorm(n, mu.t, sigma.t)
#' # Run smash
#' mu.est <- smash(X.s, "gaus")
#' # Plot the true mean function as well as the estimated one
#' plot(mu.t, type = "l")
#' lines(mu.est, col = 2)
#' 
#' # Poisson case
#' # Scale the signal to be non-zero and to have a low average intensity
#' mu.t = 0.01 + mu.s
#' # Simulate an example dataset
#' X.s = rpois(n, mu.t)
#' # Run smash
#' mu.est = smash(X.s, "poiss")
#' # Plot the true mean function as well as the estimated one
#' plot(mu.t, type = "l")
#' lines(mu.est, col = 2) 
#' @docType package
#' @name smashr
NULL
