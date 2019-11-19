#' @title Estimate the underlying mean or intensity function from
#'   Gaussian or Poisson data, respectively.
#'
#' @description This is a wrapper function for
#'   \code{\link{smash.gaus}} or \code{\link{smash.poiss}} as
#'   appropriate. For details see \code{\link{smash.gaus}} and
#'   \code{\link{smash.poiss}}.
#'
#' @details Performs nonparametric regression on univariate Poisson or
#'   Gaussian data using wavelets. For the Poisson case, the data are
#'   assumed to be i.i.d. from an underlying inhomogeneous mean
#'   function that is "smooth". Similarly for the Gaussian case, the
#'   data are assumed to be independent with an underlying smooth mean
#'   function.  In the Gaussian case, the variances are allowed vary,
#'   but are similarly "spatially structured" as with the mean
#'   function. The functions \code{ashsmooth.gaus} and
#'   \code{ashsmooth.pois} perform smoothing for Gaussian and Poisson
#'   data respectively. The only required input is a vector of length
#'   2^J for some integer J. Other options include the possibility of
#'   returning the posterior variances, specifying a wavelet basis
#'   (default is Haar, which performs well in general due to the fact
#'   that we used the translation-invariant version).
#'
#' @param x A vector of observations. Reflection is done automatically
#'   if length of \code{x} is not a power of 2.
#'
#' @param model Specifies the model (Gaussian or Poisson). Can be
#'   NULL, in which case the Poisson model is assumed if x consists of
#'   integers, and the Gaussian model is assumed otherwise. One of
#'   'gaus' or 'poiss' can also be specified to fit a specific model.
#'
#' @return See \code{smash.gaus} or \code{smash.poiss} for details.
#'
#' @examples
#'
#' # First Gaussian example
#' # ----------------------
#' # In this first example, the length of the signal is a power of 2.
#' #
#' # Create the baseline mean function. (The "spikes" function is used
#' # as an example here.)
#' set.seed(2)
#' n <- 2^9
#' t <- 1:n/n
#' spike.f <- function (x) (0.75 * exp(-500 * (x - 0.23)^2) +
#'   1.5 * exp(-2000 * (x - 0.33)^2) + 3 * exp(-8000 * (x - 0.47)^2) +
#'   2.25 * exp(-16000 * (x - 0.69)^2) + 0.5 * exp(-32000 * (x - 0.83)^2))
#' mu.s <- spike.f(t)
#'
#' # Scale the signal to be between 0.2 and 0.8
#' mu.t <- (1 + mu.s)/5
#' plot(mu.t,type = "l")
#'
#' # Create the baseline variance function. (The function V2 from Cai &
#' # Wang (2008) is used here.)
#' var.fn <- (1e-04 + 4 * (exp(-550 * (t - 0.2)^2) +
#'                         exp(-200 * (t - 0.5)^2) +
#'                         exp(-950 * (t - 0.8)^2)))/1.35
#' plot(var.fn,type = "l")
#'
#' # Set the signal-to-noise ratio.
#' rsnr    <- sqrt(5)
#' sigma.t <- sqrt(var.fn)/mean(sqrt(var.fn)) * sd(mu.t)/rsnr^2
#'
#' # Simulate an example dataset.
#' X.s <- rnorm(n,mu.t,sigma.t)
#'
#' # Run smash (Gaussian version is run since observations are not
#' # counts).
#' mu.est <- smash(X.s)
#'
#' # Plot the true mean function as well as the estimated one.
#' plot(mu.t,type = "l")
#' lines(mu.est,col = 2)
#'
#' # First Poisson example
#' # ---------------------
#' # Scale the signal to be non-zero and to have a low average intensity.
#' mu.t <- 0.01 + mu.s
#'
#' # Simulate an example dataset.
#' X.s <- rpois(n,mu.t)
#'
#' # Run smash (the Poisson version is run since observations are counts).
#' mu.est <- smash(X.s)
#'
#' # Plot the true mean function as well as the estimated one.
#' plot(mu.t,type = "l")
#' lines(mu.est,col = 2)
#'
#' # Second Gaussian example
#' # -----------------------
#' # In this second example, we demonstrate that smash also works even
#' # when the signal length (n) is not a power of 2.
#' n    <- 1000
#' t    <- 1:n/n
#' mu.s <- spike.f(t)
#'
#' # Scale the signal to be between 0.2 and 0.8.
#' mu.t <- (1 + mu.s)/5
#'
#' # Create the baseline variance function.
#' var.fn  <- (1e-04 + 4 * (exp(-550 * (t - 0.2)^2) +
#'                          exp(-200 * (t - 0.5)^2) +
#'                          exp(-950 * (t - 0.8)^2)))/1.35
#' sigma.t <- sqrt(var.fn)/mean(sqrt(var.fn)) * sd(mu.t)/rsnr^2
#'
#' # Simulate an example dataset.
#' X.s <- rnorm(n,mu.t,sigma.t)
#'
#' # Run smash.
#' mu.est <- smash(X.s)
#'
#' # Plot the true mean function as well as the estimated one.
#' plot(mu.t,type = "l")
#' lines(mu.est,col = 2)
#'
#' # Second Poisson example
#' # ----------------------
#' # The Poisson version of smash also works with signals that are not
#' # exactly of length 2^J for some integer J.
#' #
#' # Scale the signal to be non-zero and to have a low average intensity.
#' mu.t <- 0.01 + mu.s
#'
#' # Simulate an example dataset
#' X.s <- rpois(n,mu.t)
#'
#' # Run smash (Poisson version is run since observations are counts).
#' mu.est <- smash(X.s)
#'
#' # Plot the true mean function as well as the estimated one.
#' plot(mu.t,type = "l")
#' lines(mu.est,col = 2)
#'
#' @export
#'
smash = function (x, model = NULL, ...) {
  if(!is.null(model)){
    if (!(model == "gaus" | model == "poiss")) {
      stop("Error: model must be NULL or one of 'gaus' or 'poiss'")
    }
  }

  if (is.null(model)){
    if (!isTRUE(all.equal(trunc(x),x))) {
      model = "gaus"
    } else {
      model = "poiss"
    }
  }

  if (model == "gaus") {
    return(smash.gaus(x, ...))
  }

  if (model == "poiss") {
    return(smash.poiss(x, ...))
  }
}

# @description Computes the non-decimated wavelet transform matrix for
#   a given basis.
# @param n The sample size. Must be a power of 2.
# @param filter.number Specifies the type of wavelet basis used.
# @param family Specifies the type of wavelet basis used.
# @return The NDWT matrix for the specified basis, with the entries squared.
ndwt.mat = function (n, filter.number, family) {
  J = log2(n)
  X = diag(rep(1, n))
  W = matrix(0, n * J, n)
  W = apply(X, 1, wd.D, filter.number = filter.number, family = family,
            type = "station")
  return(list(W2 = W^2))
}

# A wrapper function for code in mu.smooth and var.smooth.
#
#' @importFrom ashr ash
shrink.wc = function (wc, wc.var.sqrt, ashparam, jash, df, SGD) {
    if (jash == FALSE) {
      zdat.ash = withCallingHandlers(do.call(ash,
        c(list(betahat=wc, sebetahat=wc.var.sqrt), ashparam)))
    } else {
      zdat.ash = jasha(wc, wc.var.sqrt, df = df, SGD = SGD)
    }
    return(zdat.ash)
}

# Returns "mu.est" if posterior variances are not computed, and a list
# with elements "mu.est" and "mu.est.var" otherwise.
#
#' @importFrom stats dnorm
#' @importFrom wavethresh accessD
#' @importFrom wavethresh putD
#' @importFrom wavethresh AvBasis
#' @importFrom wavethresh convert
#' @importFrom ashr calc_loglik
#' @importFrom ashr get_fitted_g
#' @importFrom ashr set_data
#' @importFrom ashr get_pm
#' @importFrom ashr get_psd
#'
mu.smooth = function (wc, data.var, basis, tsum, Wl, return.loglr,
                      post.var, ashparam, J, n) {
    wmean = matrix(0, J, n)
    wvar = matrix(0, J, n)
    if (basis[[1]] == "haar") {
        y = wc
        vtable = cxxtitable(data.var)$sumtable
        logLR.scale = c()
        for (j in 0:(J - 1)) {
            ind.nnull = (vtable[j + 2, ] != 0)
            zdat.ash = shrink.wc(y[j + 2, ind.nnull], sqrt(vtable[j + 2,
                                 ind.nnull]), ashparam, jash = FALSE,
                                 df = NULL, SGD = FALSE)
            wmean[j + 1, ind.nnull] = get_pm(zdat.ash)/2
            wmean[j + 1, !ind.nnull] = 0
            if (return.loglr == TRUE) {
                spins = 2^(j + 1)
                logLR.temp =
                  calc_loglik(get_fitted_g(zdat.ash),
                    set_data(y[j + 2, ind.nnull],
                      sqrt(vtable[j + 2, ind.nnull]),NULL,0)) -
                      sum(dnorm(y[j + 2, ind.nnull], 0,
                                sqrt(vtable[j + 2, ind.nnull]), log = TRUE))
                logLR.scale[j + 1] = logLR.temp/spins
            }
            if (post.var == TRUE) {
                wvar[j + 1, ind.nnull] = get_psd(zdat.ash)^2/4
                wvar[j + 1, !ind.nnull] = 0
            }
        }
        wwmean = -wmean
        mu.est = cxxreverse_gwave(tsum, wmean, wwmean)
        if (return.loglr == TRUE) {
            logLR = sum(logLR.scale)
        }
        if (post.var == TRUE) {
            wwvar = wvar
            mu.est.var = cxxreverse_gvwave(0, wvar, wwvar)
        }
    } else {
        x.w = wc
        # Diagonal of W*V*W'.
        x.w.v = apply((rep(1, n * J) %o% data.var) * Wl$W2, 1, sum)
        x.pm = rep(0, n)
        x.w.v.s = rep(0, n * J)
        logLR.scale = 0
        for (j in 0:(J - 1)) {
            index = (((J - 1) - j) * n + 1):((J - j) * n)
            x.w.j = wavethresh::accessD(x.w, j)
            x.w.v.j = x.w.v[index]
            ind.nnull = (x.w.v.j != 0)
            zdat.ash = shrink.wc(x.w.j[ind.nnull],
              sqrt(x.w.v.j[ind.nnull]), ashparam, jash = FALSE,
                df = NULL, SGD = FALSE)
            x.pm[ind.nnull] = get_pm(zdat.ash)
            x.pm[!ind.nnull] = 0
            x.w = wavethresh::putD(x.w, j, x.pm)
            if (return.loglr == TRUE) {
                spins = 2^(J - j)
                logLR.temp = calc_loglik(get_fitted_g(zdat.ash),
                  set_data(x.w.j[ind.nnull],
                           sqrt(x.w.v.j[ind.nnull]), NULL, 0)) -
                    sum(dnorm(x.w.j[ind.nnull], 0, sqrt(x.w.v.j[ind.nnull]),
                              log = TRUE))
                logLR.scale[j + 1] = logLR.temp/spins
            }
            if (post.var == TRUE) {
              x.w.v.s[index[ind.nnull]] = get_psd(zdat.ash)^2
              x.w.v.s[index[!ind.nnull]] = 0
            }
        }
        mu.est = wavethresh::AvBasis(wavethresh::convert(x.w))
        if (return.loglr == TRUE) {
            logLR = sum(logLR.scale)
        }
        if (post.var == TRUE) {
            mv.wd = wd.var(rep(0, n), filter.number = basis$filter.number,
                           family = basis$family, type = "station")
            mv.wd$D = x.w.v.s
            mu.est.var = AvBasis.var(convert.var(mv.wd))
        }
    }
    if (return.loglr == TRUE & post.var == TRUE) {
        return(list(mu.est = mu.est, mu.est.var = mu.est.var, logLR = logLR))
    } else if(return.loglr == TRUE & post.var == FALSE) {
        return(list(mu.est = mu.est, logLR = logLR))
    } else if(return.loglr == FALSE & post.var == TRUE) {
        return(list(mu.est = mu.est, mu.est.var = mu.est.var))
    } else {
        return(mu.est)
    }
}

# Returns "var.est" if posterior variances are not computed, and a
# list with elements 'var.est' and 'var.est.var' otherwise.
#
#' @importFrom ashr get_pm
#' @importFrom ashr get_psd
#' @importFrom wavethresh wd
var.smooth = function (data, data.var, x.var.ini, basis, v.basis, Wl,
                       filter.number, family, post.var, ashparam, jash,
                       weight, J, n, SGD) {
    wmean = matrix(0, J, n)
    wvar = matrix(0, J, n)
    if (basis[[1]] == "haar" | v.basis == FALSE) {
        vtable = cxxtitable(data.var)$sumtable
        vdtable = cxxtitable(data)$difftable
        for (j in 0:(J - 1)) {
            ind.nnull = (vtable[j + 2, ] != 0)
            zdat.ash = shrink.wc(vdtable[j + 2, ind.nnull],
                                 sqrt(vtable[j + 2, ind.nnull]),
                                ashparam, jash = jash,
                                 df = min(50, 2^(j + 1)), SGD = SGD)
            wmean[j + 1, ind.nnull] = get_pm(zdat.ash)/2
            wmean[j + 1, !ind.nnull] = 0
            if ((sum(is.na(wmean[j + 1, ])) > 0) & (SGD == TRUE)) {
                zdat.ash = shrink.wc(vdtable[j + 2, ind.nnull],
                                     sqrt(vtable[j + 2, ind.nnull]),
                                     ashparam, jash = jash,
                  df = min(50, 2^(j + 1)), SGD = FALSE)
                wmean[j + 1, ind.nnull] = get_pm(zdat.ash)/2
                wmean[j + 1, !ind.nnull] = 0
            }
            if (post.var == TRUE) {
                wvar[j + 1, ind.nnull] = get_psd(zdat.ash)^2/4
                wvar[j + 1, !ind.nnull] = 0
            }
        }
        wwmean = -wmean
        var.est = cxxreverse_gwave(weight * sum(data) + (1 - weight) *
                                   sum(x.var.ini), wmean, wwmean)
        if (post.var == TRUE) {
            wwvar = wvar
            var.est.var = cxxreverse_gvwave(0, wvar, wwvar)
        }
    } else {
        x.w = wavethresh::wd(data, filter.number = filter.number,
                             family = family, type = "station")

        # Diagonal of W*V*W'.
        x.w.v = apply((rep(1, n * J) %o% data.var) * Wl$W2, 1, sum)
        x.pm = rep(0, n)
        x.w.v.s = rep(0, n * J)
        for (j in 0:(J - 1)) {
            index = (((J - 1) - j) * n + 1):((J - j) * n)
            x.w.j = accessD(x.w, j)
            x.w.v.j = x.w.v[index]
            ind.nnull = (x.w.v.j != 0)
            zdat.ash = shrink.wc(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]),
                                 ashparam, jash = jash,
                df = min(50, 2^(j + 1)), SGD = SGD)
            x.pm[ind.nnull] = get_pm(zdat.ash)
            x.pm[!ind.nnull] = 0
            if ((sum(is.na(x.pm)) > 0) & (SGD == TRUE)) {
                zdat.ash =
                   shrink.wc(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]),
                             ashparam, jash = jash,
                  df = min(50, 2^(j + 1)), SGD = FALSE)
                x.pm[ind.nnull] = get_pm(zdat.ash)
                x.pm[!ind.nnull] = 0
            }
            x.w = putD(x.w, j, x.pm)
            if (post.var == TRUE) {
                x.w.v.s[index[ind.nnull]] = get_psd(zdat.ash)^2
                x.w.v.s[index[!ind.nnull]] = 0
            }
        }
        var.est = AvBasis(convert(x.w))
        if (post.var == TRUE) {
            mv.wd = wd.var(rep(0, n), filter.number = basis$filter.number,
                           family = basis$family, type = "station")
            mv.wd$D = x.w.v.s
            var.est.var = AvBasis.var(convert.var(mv.wd))
        }
    }
    if (post.var == TRUE) {
        return(list(var.est = var.est, var.est.var = var.est.var))
    } else {
        return(var.est)
    }
}

# Set default ash parameters. By default, ashparam$df = NULL,
# ashparam$mixsd = NULL and ashparam$g = NULL.
#
#' @importFrom utils modifyList
setAshParam.gaus = function (ashparam) {

  if (!is.list(ashparam))
    stop("Error: invalid parameter 'ashparam'")
  ashparam.default = list(pointmass = TRUE, prior = "nullbiased",
                          gridmult = 2, mixcompdist = "normal",
                          nullweight = 10, outputlevel = 2, fixg = FALSE)
  ashparam = modifyList(ashparam.default, ashparam)
  if (!is.null(ashparam[["g"]]))
    stop(paste("Error: ash parameter 'g' can only be NULL; if you want",
               "to specify ash parameter 'g' use multiseq arguments",
               "'fitted.g' and/or 'fitted.g.intercept'"))

  if (!((is.null(ashparam[["mixsd"]])) |
        (is.numeric(ashparam[["mixsd"]]) &
         (length(ashparam[["mixsd"]]) < 2))))
    stop(paste("Error: invalid parameter 'mixsd', 'mixsd' must be",
               "null or a numeric vector of length >=2"))
  if (!((ashparam[["prior"]] == "nullbiased") |
        (ashparam[["prior"]] == "uniform") |
        is.numeric(ashparam[["prior"]])))
    stop(paste("Error: invalid parameter 'prior', 'prior' can be",
               "a number or 'nullbiased' or 'uniform'"))
  return(ashparam)
}

#' @title Estimate underlying mean function from noisy Gaussian data.
#'
#' @description This function takes a data vector as input performs
#'   signal denoising using wavelet decomposition and an adaptive
#'   shrinkage prior on the wavelet parameters. The data are assumed to
#'   be (mostly) independent and Gaussian, but not necessarily
#'   identically distributed.
#'
#' @details We assume that the data come from the model \eqn{Y_t =
#'   \mu_t + \epsilon_t} for \eqn{t=1,...,T}, where \eqn{\mu_t} is an
#'   underlying mean, assumed to be spatially structured (or treated as
#'   points sampled from a smooth continous function), and
#'   \eqn{\epsilon_t \sim N(0, \sigma_t)}, and are independent. Smash
#'   provides estimates of \eqn{\mu_t} and \eqn{\sigma_t^2} (and their
#'   posterior variances if desired).
#'
#' @param x A vector of observations. Reflection is done automatically
#'   if length of \code{x} is not a power of 2.
#'
#' @param sigma A vector of standard deviations. Can be provided if
#'   known or estimated beforehand.
#'
#' @param v.est Boolean indicating if variance estimation should be
#'   performed instead.
#'
#' @param joint Boolean indicating if results of mean and variance
#'   estimation should be returned together.
#'
#' @param v.basis Boolean indicating if the same wavelet basis should
#'   be used for variance estimation as mean estimation. If false,
#'   defaults to Haar basis for variance estimation (this is much faster
#'   than other bases).
#'
#' @param post.var Boolean indicating if the posterior variance should
#'   be returned for the mean and/or variance estiamtes.
#'
#' @param family Choice of wavelet basis to be used, as in
#'   \code{wavethresh}.
#'
#' @param filter.number Choice of wavelet basis to be used, as in
#'   \code{wavethresh}.
#'
#' @param return.loglr Boolean indicating if a logLR should be returned.
#'
#' @param jash Indicates if the prior from method JASH should be
#'   used. This will often provide slightly better variance estimates
#'   (especially for nonsmooth variance functions), at the cost of
#'   computational efficiency. Defaults to FALSE.
#'
#' @param SGD Boolean indicating if stochastic gradient descent should
#'   be used in the EM. Only applicable if jash=TRUE.
#'
#' @param weight Optional parameter used in estimating overall
#'   variance. Only works for Haar basis. Defaults to 0.5. Setting this
#'   to 1 might improve variance estimation slightly.
#'
#' @param min.var The minimum positive value to be set if the
#'   variance estimates are non-positive.
#'
#' @param ashparam A list of parameters to be passed to \code{ash};
#'   default values are set by function \code{setAshParam.gaus}.
#'
#' @param homoskedastic indicates whether to assume constant variance
#'   (if v.est is true)
#'
#' @param reflect A logical indicating if the signals should be
#'   reflected.
#'
#' @return \code{smash.gaus} returns the following by default:
#'
#'   \item{mu.res}{A list with the mean estimate, its posterior variance
#'   if \code{post.var} is TRUE, the logLR if \code{return.loglr} is
#'   TRUE, or a vector of mean estimates if neither \code{post.var} nor
#'   \code{return.loglr} are TRUE.}
#'
#'   If \code{v.est} is TRUE, then \code{smash.gaus} returns the
#'   following:
#'
#'   \item{var.res}{A list with the variance estimate, its posterior
#'   variance if \code{post.var} is TRUE, or a vector of variance
#'   estimates if \code{post.var} is \code{FALSE} In addition, if
#'   \code{joint} is TRUE, then both \code{mu.res} and \code{var.res}
#'   are returned.}
#'
#' @examples
#'
#' n=2^10
#' t=1:n/n
#' spike.f = function(x) (0.75*exp(-500*(x-0.23)^2) +
#'   1.5*exp(-2000*(x-0.33)^2) + 3*exp(-8000*(x-0.47)^2) +
#'   2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
#' mu.s = spike.f(t)
#'
#' # Gaussian case
#' mu.t = (1+mu.s)/5
#' plot(mu.t,type='l')
#' var.fn = (0.0001 + 4*(exp(-550*(t-0.2)^2) + exp(-200*(t-0.5)^2) +
#'   exp(-950*(t-0.8)^2)))/1.35
#' plot(var.fn,type='l')
#' rsnr=sqrt(5)
#' sigma.t=sqrt(var.fn)/mean(sqrt(var.fn))*sd(mu.t)/rsnr^2
#' X.s=rnorm(n,mu.t,sigma.t)
#' mu.est=smash.gaus(X.s)
#' plot(mu.t,type='l')
#' lines(mu.est,col=2)
#'
#' @importFrom wavethresh wd
#'
#' @export
#'
smash.gaus = function (x, sigma = NULL, v.est = FALSE, joint = FALSE,
                       v.basis = FALSE, post.var = FALSE, filter.number = 1,
                       family = "DaubExPhase", return.loglr = FALSE,
                       jash = FALSE, SGD = TRUE, weight = 0.5,
                       min.var = 1e-08, ashparam = list(),
                       homoskedastic = FALSE, reflect=FALSE) {

   # Check inputs x.
   if (!(is.numeric(x) & length(x) > 1))
     stop("Argument \"x\" should be a numeric vector with more than one element")

   # Check input sigma.
   if (!is.null(sigma))
     if (!(is.numeric(sigma) & (length(sigma) == 1 | length(sigma) == length(x))))
       stop("Argument \"sigma\" should be NULL or a scalar or a numeric vector of the same length as \"x\"")

   if (reflect | !ispowerof2(length(x))) {

     # Reflect variances.
     if(!is.null(sigma) & length(sigma) > 1) {
       reflect.sigma = reflect(sigma)
       sigma = reflect.sigma$x
     }

     # Reflect data.
     reflect.res = reflect(x)
     idx         = reflect.res$idx
     x           = reflect.res$x

   } else {
      idx = 1:length(x)
    }
    J = log2(length(x))
    n = length(x)

   if (length(sigma) == 1) {
      sigma = rep(sigma, n)
   }

    if(!v.est & homoskedastic){
        stop("Error: can't set homoskedastic TRUE if not estimating variance")
    }
    if(v.est & homoskedastic){
      sigma = sd_estimate_gasser_etal(x)
      v.est = FALSE;
    }

     ashparam = setAshParam.gaus(ashparam)


    if (v.est) {
        weight = 1
    } else {
        weight = 0.5
    }

    if (filter.number == 1 & family == "DaubExPhase") {
        basis = "haar"
    } else {
        basis = list(family = family, filter.number = filter.number)
    }
    if (post.var == TRUE & v.est == TRUE & jash == TRUE) {
        stop(paste("Error: Posterior variances for variance estimate",
                   "not returned for method JASH"))
    }
    if (joint == TRUE) {
        v.est = TRUE
    }

    tsum = sum(x)
    x.w.d = cxxtitable(x)$difftable
    Wl = NULL
    if (basis[[1]] != "haar") {
        x.w.d = wavethresh::wd(x, filter.number = filter.number,
                               family = family, type = "station")
        Wl = ndwt.mat(n, filter.number = filter.number, family = family)
    }

    if (is.null(sigma)) {
        var.est1.ini = (x - lshift(x))^2/2
        var.est2.ini = (rshift(x) - x)^2/2
        var.est.ini = (var.est1.ini + var.est2.ini)/2
        mu.est = mu.smooth(x.w.d, var.est.ini, basis, tsum, Wl, FALSE,
                           FALSE, ashparam, J, n)
        var.est = (x - mu.est)^2
        var.var.est = 2/3 * var.est^2
        var.est = var.smooth(var.est, var.var.est, var.est.ini, basis,
                             v.basis, Wl, filter.number, family, FALSE,
                             ashparam, jash, weight, J, n, SGD = SGD)
        var.est[var.est <= 0] = 1e-08
        sigma = sqrt(var.est)
    }

    ashparam.mean = ashparam
    ashparam.mean$gridmult = 64
    mu.res = mu.smooth(x.w.d, sigma^2, basis, tsum, Wl, return.loglr,
                       post.var, ashparam.mean, J, n)

    if (!v.est) {

      # Truncate the reflected x back to the original length m.
      if (post.var == TRUE) {
        mu.res$mu.est = mu.res$mu.est[idx]
        mu.res$mu.est.var = mu.res$mu.est.var[idx]
        return(mu.res)
      }else if(return.loglr == TRUE){
        mu.res$mu.est = mu.res$mu.est[idx]
        return(mu.res)
      }else{
        return(mu.res[idx])
      }


    } else {
        if (return.loglr == FALSE & post.var == FALSE) {
            mu.est = mu.res
        } else {
            mu.est = mu.res$mu.est
        }
        var.est = (x - mu.est)^2
        var.var.est = 2/3 * var.est^2
        var.res = var.smooth(var.est, var.var.est, 0, basis, v.basis, Wl,
                             filter.number, family, post.var, ashparam,
                             jash, 1, J, n, SGD = SGD)
        if (post.var == FALSE) {
            var.res[var.res <= 0] = min.var

            # Truncate the reflected x back to the original length m.
            var.res = var.res[idx]
            mu.res = mu.res[idx]
        } else {
            var.res$var.est[var.res$var.est <= 0] = min.var
            # Truncate the reflected x back to the original length m.
            mu.res$mu.est       = mu.res$mu.est[idx]
            mu.res$mu.est.var   = mu.res$mu.est.var[idx]
            var.res$var.est     = var.res$var.est[idx]
            var.res$var.est.var = var.res$var.est.var[idx]
        }
        if (joint == FALSE) {
            return(var.res = var.res)
        } else {
            return(list(mu.res = mu.res, var.res = var.res))
        }
    }
}

#' @title Estimate homoskedastic standard deviation from nonparamatric
#'    regression.
#'
#' @param x The data.
#'
#' @return An estimate of the standard deviation
#'
#' @details Uses formula (3) from Brown and Levine (2007), Annals of
#'   Statistics, who attribute it to Gasser et al.
#'
#' @export
#'
sd_estimate_gasser_etal = function(x){
  n = length(x)
  sqrt(2/(3 * (n - 2)) * sum((1/2 * x[1:(n - 2)] - x[2:(n - 1)] +
                              1/2 * x[3:n])^2))
}

#' @title TI thresholding with heteroskedastic errors.
#'
#' @param x The data. Should be a vector of length a power of 2.
#'
#' @param sigma The standard deviation function. Can be provided if
#'   known or estimated beforehand.
#'
#' @param method The method to estimate the variance function. Can be
#'   'rmad' for running MAD as described in Gao (1997), or 'smash'.
#'
#' @param filter.number The wavelet basis to be used.
#'
#' @param family The wavelet basis to be used.
#'
#' @param min.level The primary resolution level.
#'
#' @return returns a vector of mean estimates
#'
#' @details The 'rmad' option effectively implements the procedure
#'   described in Gao (1997), while the 'smash' option first estimates
#'   the variance function using package \code{smash} and then performs
#'   thresholding given this variance function.
#'
#' @references Gao, Hong-Ye (1997) Wavelet shrinkage estimates for
#'   heteroscedastic regression models. MathSoft, Inc.
#'
#' @examples
#'
#' n=2^10
#' t=1:n/n
#' spike.f = function(x) (0.75*exp(-500*(x-0.23)^2) +
#'   1.5*exp(-2000*(x-0.33)^2) +
#'   3*exp(-8000*(x-0.47)^2) +
#'   2.25*exp(-16000*(x-0.69)^2) +
#'   0.5*exp(-32000*(x-0.83)^2))
#' mu.s = spike.f(t)
#'
#' # Gaussian case
#' # -------------
#' mu.t=(1+mu.s)/5
#' plot(mu.t,type='l')
#' var.fn = (0.0001+4*(exp(-550*(t-0.2)^2) +
#'   exp(-200*(t-0.5)^2) + exp(-950*(t-0.8)^2)))/1.35
#' plot(var.fn,type='l')
#' rsnr=sqrt(5)
#' sigma.t=sqrt(var.fn)/mean(sqrt(var.fn))*sd(mu.t)/rsnr^2
#' X.s=rnorm(n,mu.t,sigma.t)
#' mu.est.rmad<-ti.thresh(X.s,method='rmad')
#' mu.est.smash<-ti.thresh(X.s,method='smash')
#' plot(mu.t,type='l')
#' lines(mu.est.rmad,col=2)
#' lines(mu.est.smash,col=4)
#'
#' @importFrom wavethresh wd
#' @importFrom caTools runmad
#'
#' @export
#'
ti.thresh = function (x, sigma = NULL, method = "smash", filter.number = 1,
                      family = "DaubExPhase", min.level = 3,
                      ashparam = list()) {
    n = length(x)
    J = log2(n)
    if (length(sigma) == 1)
        sigma = rep(sigma, n)

    if (is.null(sigma) & method == "rmad") {
        lambda.thresh = 1
    } else {
        lambda.thresh = 2
    }

    if (is.null(sigma)) {
        if (method == "smash") {
            sigma = sqrt(smash.gaus(x, v.est = TRUE, v.basis = TRUE,
                                    filter.number = filter.number,
                                    family = family, ashparam = ashparam,
                                    weight = 1))
        } else if (method == "rmad") {
            x.w = wavethresh::wd(x, filter.number = filter.number,
                                 family = family, type = "station")
            win.size = round(n/10)
            odd.boo = (win.size%%2 == 1)
            win.size = win.size + (1 - odd.boo)
            sigma = runmad(accessD(x.w, J - 1), win.size,
                           endrule = "func")
        } else {
            stop("Error: Method not recognized")
        }
    }

    if (filter.number == 1 & family == "DaubExPhase") {
        tsum = sum(x)
        vdtable = cxxtitable(x)$difftable
        vtable = cxxtitable(sigma^2)$sumtable
        wmean = threshold.haar(vdtable, vtable, lambda.thresh, levels = 0:(J - 1 - min.level), type = "hard")
        wwmean = -wmean
        mu.est = cxxreverse_gwave(tsum, wmean, wwmean)
    } else {
        x.w = wavethresh::wd(x, filter.number = filter.number,
                             family = family, type = "station")
        Wl = ndwt.mat(n, filter.number = filter.number, family = family)

        # Diagonal of W*V*W'.
        x.w.v = apply((rep(1, n * J) %o% (sigma^2)) * Wl$W2, 1, sum)
        x.w.t = threshold.var(x.w, x.w.v, lambda.thresh,
                              levels = (min.level):(J - 1),
                              type = "hard")
        mu.est = AvBasis(convert(x.w.t))
    }
    return(mu.est)
}

#' @title Reflect and extend a vector.
#'
#' @description Reflects a vector if it has length a power of 2;
#'   otherwise extends the vector to have length a power of 2 and then
#'   reflects it.
#'
#' @param x An n-vector.
#'
#' @return A list with two list elements: \code{"x"} containing the
#'   reflected signal; and \code{"idx"} containing the indices of the
#'   original signal.
#'
#' @export
#'
reflect <- function (x) {
    n = length(x)
    J = log2(n)
    if ((J%%1) == 0) {

        # if J is an integer, i.e. n is a power of 2.
        x = c(x, x[n:1])
        return(list(x=x, idx = 1:n))
    } else {
        n.ext = 2^ceiling(J)
        lnum = round((n.ext - n)/2)
        rnum = n.ext - n - lnum
        if (lnum == 0) {
            x.lmir = NULL
        } else {
            x.lmir = x[lnum:1]
        }
        if (rnum == 0) {
            x.rmir = NULL
        } else {
            x.rmir = x[n:(n - rnum + 1)]
        }
        x.ini = c(x.lmir, x, x.rmir)
        x.mir = x.ini[n.ext:1]
        x = c(x.ini, x.mir)
        return(list(x = x, idx = (lnum + 1):(lnum + n)))
    }
}

# Compute posterior mean and var for log(p), log(q), log(p0/p1) and
# log(q0/q1) This function returns posterior means and variances of
# log(p), log(q), log(p0/p1) and log(q0/q1) as lp, lq, lpratio and
# lqratio, respectively, where p is the probability of going left and
# q=1-p.
#
# Input log indicates if mean estimation is to be performed on the log
# scale.
#
# The return value is a list with elements lp.mean, lp.var, lq.mean
# and lq.var.
compute.res <- function (alpha, log) {
    if (log == TRUE) {

        # Find mean and variance of log(q).
        lp = ff.moments(alpha$mean, alpha$var)
        lq = ff.moments(-alpha$mean, alpha$var)
        return(list(lp.mean = lp$mean, lp.var = lp$var,
                    lq.mean = lq$mean, lq.var = lq$var))
    } else {
        lp = ff.moments_exp(alpha$mean, alpha$var)
        lq = ff.moments_exp(-alpha$mean, alpha$var)
        return(list(lp.mean = lp$mean, lp.var = lp$meansq,
                    lq.mean = lq$mean, lq.var = lq$meansq))
    }
}

# For each resolution, performs shrinkage and returns the posterior
# means and variances in matrix form.
#
#' @importFrom data.table rbindlist
#' @importFrom ashr ash
#' @importFrom ashr get_pm
#' @importFrom ashr get_psd
#'
getlist.res = function (res, j, n, zdat, log, shrink, ashparam) {
  ind = ((j - 1) * n + 1):(j * n)
  if (shrink == TRUE) {

    # Apply ash to vector of intercept estimates and SEs.
    zdat.ash = withCallingHandlers(do.call(ash,
        c(list(betahat = zdat[1, ind], sebetahat = zdat[2, ind]), ashparam)))

    # Find mean and variance of alpha.
    alpha.mv = list(mean = get_pm(zdat.ash),
      var = get_psd(zdat.ash)^2)
    } else {
      alpha.mv = list(mean = fill.nas(zdat[1, ind]),
          var = fill.nas(zdat[2, ind])^2)
  }
  res.j = compute.res(alpha.mv, log)
  res = rbindlist(list(res, res.j))
  return(res)
}

# Reconstructs the signal in data space.
#
# @param ls estimated total intensity
# @param res the matrix of wavelet proportions/probabilities
# @param log bool, indicating if signal reconstruction should be in
#   log space.
# @param n length of data
# @param J log2(n)
#
# Return reconstructed signal in original data space.
recons.mv = function (ls, res, log, n, J) {
    if (log == TRUE) {

        # Reconstructs estimate from the 'wavelet' space on the log level.
        est.mean = cxxreverse_pwave(log(ls), matrix(res$lp.mean, J, n,
                     byrow = TRUE), matrix(res$lq.mean, J, n, byrow = TRUE))
        est.var = cxxreverse_pwave(0, matrix(res$lp.var, J, n, byrow = TRUE),
                     matrix(res$lq.var, J, n, byrow = TRUE))
    } else {

        # Reconstruction on non-log level.
        est.mean = exp(cxxreverse_pwave(log(ls), log(matrix(res$lp.mean, J,
            n, byrow = TRUE)), log(matrix(res$lq.mean, J, n, byrow = TRUE))))
        est.ms = exp(cxxreverse_pwave(2 * log(ls),
            log(matrix(res$lp.var, J, n, byrow = TRUE)),
            log(matrix(res$lq.var, J, n, byrow = TRUE))))
        est.var = pmax(est.ms - est.mean^2, 0)
    }
    return(list(est.mean = est.mean, est.var = est.var))
}

# Set default ash parameters. By default, ashparam$df = NULL,
# ashparam$mixsd = NULL and ashparam$g = NULL.
#
#' @importFrom utils modifyList
setAshParam.poiss = function (ashparam) {
    if (!is.list(ashparam))
        stop("Error: invalid parameter 'ashparam'")
    ashparam.default = list(pointmass = TRUE, prior = "nullbiased",
                            gridmult = 2, mixcompdist = "normal",
                            nullweight = 10, outputlevel = 2, fixg = FALSE)
    ashparam = modifyList(ashparam.default, ashparam)
    if (!is.null(ashparam[["g"]]))
      stop(paste("Error: ash parameter 'g' can only be NULL; if you want",
                 "to specify ash parameter 'g' use multiseq arguments",
                 "'fitted.g' and/or 'fitted.g.intercept'"))

    if(!((is.null(ashparam[["mixsd"]])) |
         (is.numeric(ashparam[["mixsd"]])
          & (length(ashparam[["mixsd"]]) < 2))))
      stop(paste("Error: invalid parameter 'mixsd', 'mixsd' must be null",
                 "or a numeric vector of length >=2"))
    if(!((ashparam[["prior"]] == "nullbiased") |
         (ashparam[["prior"]] == "uniform") |
         is.numeric(ashparam[["prior"]])))
      stop(paste("Error: invalid parameter 'prior', 'prior' can be a",
                 "number or 'nullbiased' or 'uniform'"))
    return(ashparam)
}

# Set default glm.approx parameters.
#
#' @importFrom utils modifyList
setGlmApproxParam = function (glm.approx.param) {
  if(!is.list(glm.approx.param))
    stop("Error: invalid parameter 'glm.approx.param'")
  glm.approx.param.default =
    list(minobs = 1, pseudocounts = 0.5, all = FALSE, center = FALSE,
         forcebin = TRUE, repara = TRUE, lm.approx = FALSE, disp = "add")
  glm.approx.param = modifyList(glm.approx.param.default, glm.approx.param)

  if(glm.approx.param[["minobs"]]%%1 != 0|glm.approx.param[["minobs"]] <1 )
    stop("Error: minobs must be an integer larger than or equal to 1")
  if(!(glm.approx.param[["disp"]] %in% c("add", "mult")))
    stop("Error: parameter disp must be 'add' or 'mult'")
  if(!(is.numeric(glm.approx.param[["pseudocounts"]]) &
       glm.approx.param[["pseudocounts"]] > 0))
    stop(paste("Error: invalid parameter 'pseudocounts',",
               "'pseudocounts' must be a positive number"))

  return(glm.approx.param)
}

#' @title Estimate the underlying intensity for Poisson counts.
#'
#' @description Main smoothing procedure for Poisson data. Takes a
#'   univariate inhomogeneous Poisson process and estimates its mean
#'   intensity.
#'
#' @details We assume that the data come from the model \eqn{Y_t \sim
#'   Pois(\mu_t)} for \eqn{t=1,...,T}, where \eqn{\mu_t} is the
#'   underlying intensity, assumed to be spatially structured (or
#'   treated as points sampled from a smooth continous function). The
#'   \eqn{Y_t} are assumed to be independent. Smash provides estimates
#'   of \eqn{\mu_t} (and its posterior variance if desired).
#'
#' @param x A vector of Poisson counts (reflection is done
#'   automatically if length of \code{x} is not a power of 2).
#'
#' @param post.var Boolean, indicates if the posterior variance should
#'   be returned.
#'
#' @param log bool, determines if smoothed signal is returned on log
#'   scale or not
#'
#' @param reflect A logical indicating if the signals should be
#'   reflected.
#'
#' @param glm.approx.param A list of parameters to be passed to
#'   \code{glm.approx}; default values are set by function
#'   \code{setGlmApproxParam}.
#'
#' @param ashparam A list of parameters to be passed to \code{ash};
#'   default values are set by function \code{setAshParam.poiss}.
#'
#' @param cxx bool, indicates if C++ code should be used to create TI
#'   tables.
#'
#' @param lev integer from 0 to J-1, indicating primary level of
#'   resolution. Should NOT be used (ie shrinkage is performed at all
#'   resolutions) unless there is good reason to do so.
#'
#' @return \code{smash.poiss} returns the mean estimate by default,
#'   with the posterior variance as an additional component if
#'   \code{post.var} is TRUE.
#'
#' @examples
#'
#' n=2^10
#' t=1:n/n
#' spike.f = function(x) (0.75*exp(-500*(x-0.23)^2) +
#'   1.5*exp(-2000*(x-0.33)^2) + 3*exp(-8000*(x-0.47)^2) +
#'   2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
#' mu.s=spike.f(t)
#' mu.t=0.01+mu.s
#' X.s=rpois(n,mu.t)
#' mu.est=smash.poiss(X.s)
#' plot(mu.t,type='l')
#' lines(mu.est,col=2)
#'
#' @export
#'
smash.poiss = function (x, post.var = FALSE, log = FALSE, reflect = FALSE,
                        glm.approx.param = list(), ashparam = list(),
                        cxx = TRUE, lev = 0) {
    if (is.matrix(x)) {
        if (nrow(x) == 1) {

            # Change matrix x to vector.
            x = as.vector(x)
        } else {
            stop("x cannot have multiple rows")
        }
    }

    ashparam = setAshParam.poiss(ashparam)
    glm.approx.param = setGlmApproxParam(glm.approx.param)

    if (!is.numeric(x))
        stop("Error: invalid parameter 'x': 'x' must be numeric")
    if (!is.logical(reflect))
        stop("Error: invalid parameter 'reflect', 'reflect' must be bool")

    J = log2(length(x))

    # If ncol(x) is not a power of 2, reflect x.
    if ((J %% 1) != 0)
        {
            reflect = TRUE
        }

    # Reflect signal.
    if (reflect | !ispowerof2(length(x))) {
      reflect.res     = reflect(x)
      reflect.indices = reflect.res$idx
      x               = reflect.res$x
    }

    n = length(x)
    J = log2(n)

    # Create the parent TI table for each signal, and put into rows of
    # matrix y.
    ls = sum(x)

    if (!cxx) {
        y = as.vector(t(ParentTItable(x)$parent))
    } else {
        y = cxxSParentTItable(x)
    }

    zdat = withCallingHandlers(do.call(glm.approx, c(list(x = y, g = NULL),
                                       glm.approx.param)))

    res = list()

    # Loop through resolutions, smoothing each resolution separately
    for (j in 1:(J - lev)) {
        res = getlist.res(res, j, n, zdat, log, TRUE, ashparam)
    }

    # Do not smooth for coarser levels, with everything the same as
    # above but using the estimate and its variance as the posterior
    # mean and variance ie flat prior.
    if (lev != 0) {
        for (j in (J - lev + 1):J) {
            res = getlist.res(res, j, n, zdat, log, FALSE, ashparam)
        }
    }
    recons = recons.mv(ls, res, log, n, J)
    if (reflect==TRUE | floor(log2(length(x)))!=ceiling(log2(length(x)))) {
        recons$est.mean = recons$est.mean[reflect.indices]
        recons$est.var = recons$est.var[reflect.indices]
    }
    if (post.var == FALSE) {

        # If post.var = TRUE then only estimate (and not variance) is
        # returned.
        return(est = recons$est.mean)
    } else {
        return(list(est = recons$est.mean, var = recons$est.var))
    }
}

# Fill NAs with zeros (or other value).
fill.nas <- function (x, t = 0) {
  x[is.na(x)] <- t
  return(x)
}

# Check that x is a (non-negative) power of 2.
ispowerof2 <- function (x)
  x >= 1 & 2^ceiling(log2(x)) == x
