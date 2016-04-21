# Disclaimer: Functions 'interleave', 'lsift', 'rshift', 'simplify', 'setAshParam.gaus' and 'setAshParam.poiss' in
# this file are adapted from the same functions defined in package 'Multiseq' under the GPL license, authored by
# Ester Pantaleo, Heejung Shim, Matthew Stephens and Zhengrong Xing. 


#' @title Estimate the underlying mean or intensity function from Gaussian or Poisson data respectively
#'
#' @description This is a wrapper function for \code{\link{smash.gaus}} or \code{\link{smash.poiss}} as appropriate. For details see \code{\link{smash.gaus}} and \code{\link{smash.poiss}}
#'
#' @details Performs nonparametric regression on univariate Poisson or Gaussian data using wavelets. For the Poisson case, the data are assumed to be i.i.d. 
#' from an underlying inhomogeneous mean function that is "smooth". Similarly for the Gaussian case, the data are assumed to be independent with an underlying smooth mean function. 
#' In the Gaussian case, the variances are allowed vary, but are similarly "spatially structured" as with the mean function. The functions \code{ashsmooth.gaus} and \code{ashsmooth.pois} 
#' perform smoothing for Gaussian and Poisson data respectively. The only required input is a vector of length 2^J for some integer J.
#' Other options include the possibility of returning the posterior variances, specifying a wavelet basis (default is Haar, which performs well in general due to the fact that we used the translation-invariant version).
#' @param x: a vector of observations
#' @param model: specifies the model (Gaussian or Poisson). Can be NULL, in which case the Poisson model is assumed if x
#' consists of integers, and the Gaussian model is assumed otherwise. One of 'gaus' or 'poiss' can also be specified to fit a specific model. 
#' 
#' @return See \code{smash.gaus} or \code{smash.poiss} for details.
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
#' # Run smash (Gaussian version is run since observations are not counts)
#' mu.est <- smash(X.s)
#' # Plot the true mean function as well as the estimated one
#' plot(mu.t, type = "l")
#' lines(mu.est$mu.est, col = 2)
#' 
#' # Poisson case
#' # Scale the signal to be non-zero and to have a low average intensity
#' mu.t = 0.01 + mu.s
#' # Simulate an example dataset
#' X.s = rpois(n, mu.t)
#' # Run smash (Poisson version is run since observations are counts)
#' mu.est = smash(X.s)
#' # Plot the true mean function as well as the estimated one
#' plot(mu.t, type = "l")
#' lines(mu.est, col = 2) 
#' 
#' @export
smash = function(x, model = NULL, ...){
  if(!is.null(model)){
    if(!(model == "gaus" | model == "poiss")){
      stop("Error: model must be NULL or one of 'gaus' or 'poiss'")
    }
  }
  
  if(is.null(model)){
    if(!is.integer(x)){
      model = "gaus"
    }else{
      model = "poiss"
    }
  }
  
  if(model == "gaus"){
    return(smash.gaus(x, ...))
  }
  
  if(model == "poiss"){
    return(smash.poiss(x, ...))
  }
}



#' Computes the non-decimated wavelet transform matrix for a given basis
#' @param n: the sample size. Must be a power of 2
#' @param filter.number,family: specifies the type of wavelet basis used
#' @return the NDWT matrix for the specified basis, with the entries squared
ndwt.mat = function(n, filter.number, family) {
  J = log2(n)
  X = diag(rep(1, n))
  W = matrix(0, n * J, n)
  # for(i in 1:n){ W[,i]=wd(X[i,],filter.number,family,type='station')$D }
  W = apply(X, 1, wd.D, filter.number = filter.number, family = family, type = "station")
  # W=Matrix(W,sparse=TRUE) if(post.var==FALSE){
  return(list(W2 = W^2))
  # }else{ Winv=ginv(W) return(list(W=W,Wi=Winv)) #Wsvd=svd(W)
  # #Wmod=Wsvd$u[,1:(n-1)]%*%diag(Wsvd$d[1:(n-1)])%*%t(Wsvd$v[1:(n-1),1:(n-1)]) #return(list(W=W,Wmod=Wmod)) }
}


#' a wrapper function for code in \code{mu.smooth} and \code{var.smooth}
#' @keywords internal 
shrink.wc = function(wc, wc.var.sqrt, prior, pointmass, nullcheck, VB, mixsd, mixcompdist, gridmult, jash, df, SGD) {
    if (jash == FALSE) {
        zdat.ash = suppressWarnings(ash(wc, wc.var.sqrt, prior = prior, outputlevel=2, pointmass = pointmass, 
            VB = VB, mixsd = mixsd, mixcompdist = mixcompdist, gridmult = gridmult,
            df = NULL, control = list(trace = FALSE)))
    } else {
        zdat.ash = jasha(wc, wc.var.sqrt, df = df, SGD = SGD)
    }
    return(zdat.ash)
}


#' @title mu.smooth
#' @return 'mu.est' if posterior variances are not computed, and a list with elements 'mu.est' and 'mu.est.var' otherwise
#' @keywords internal
mu.smooth = function(wc, data.var, basis, tsum, Wl, return.loglr, post.var, prior, pointmass, nullcheck, VB, mixsd, mixcompdist, 
    gridmult, J, n) {
    wmean = matrix(0, J, n)
    wvar = matrix(0, J, n)
    if (basis[[1]] == "haar") {
        y = wc
        vtable = cxxtitable(data.var)$sumtable
        logLR.scale = 0
        for (j in 0:(J - 1)) {
            ind.nnull = (vtable[j + 2, ] != 0)
            zdat.ash = shrink.wc(y[j + 2, ind.nnull], sqrt(vtable[j + 2, ind.nnull]), prior = prior, pointmass = pointmass, 
                nullcheck = nullcheck, VB = VB, mixsd = mixsd, mixcompdist = mixcompdist, gridmult = gridmult, jash = FALSE, 
                df = NULL, SGD = FALSE)
            wmean[j + 1, ind.nnull] = zdat.ash$PosteriorMean/2
            wmean[j + 1, !ind.nnull] = 0
            if (return.loglr == TRUE) {
                spins = 2^(j + 1)
                zdat.ash$model = "EE"
                logLR.temp = ashr:::calc_loglik(zdat.ash, y[j + 2, ind.nnull], sqrt(vtable[j + 2, ind.nnull]), NULL) -  
                    sum(dnorm(y[j + 2, ind.nnull], 0, sqrt(vtable[j + 2, ind.nnull]), log = TRUE))
                logLR.scale[j + 1] = logLR.temp/spins
            }
            if (post.var == TRUE) {
                wvar[j + 1, ind.nnull] = zdat.ash$PosteriorSD^2/4
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
        x.w.v = apply((rep(1, n * J) %o% data.var) * Wl$W2, 1, sum)  #diagonal of W*V*W' 
        x.pm = rep(0, n)
        x.w.v.s = rep(0, n * J)
        logLR.scale = 0
        for (j in 0:(J - 1)) {
            index = (((J - 1) - j) * n + 1):((J - j) * n)
            x.w.j = accessD(x.w, j)
            x.w.v.j = x.w.v[index]
            ind.nnull = (x.w.v.j != 0)
            zdat.ash = shrink.wc(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]), prior = prior, pointmass = pointmass, 
                nullcheck = nullcheck, VB = VB, mixsd = mixsd, mixcompdist = mixcompdist, gridmult = gridmult, jash = FALSE, 
                df = NULL, SGD = FALSE)
            x.pm[ind.nnull] = zdat.ash$PosteriorMean
            x.pm[!ind.nnull] = 0
            x.w = putD(x.w, j, x.pm)
            if (return.loglr == TRUE) {
                spins = 2^(J - j)
                zdat.ash$model = "EE"
                logLR.temp = ashr:::calc_loglik(zdat.ash, x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]), NULL) -  
                    sum(dnorm(x.w.j[ind.nnull], 0, sqrt(x.w.v.j[ind.nnull]), log = TRUE))
                logLR.scale[j + 1] = logLR.temp/spins
            }
            if (post.var == TRUE) {
                x.w.v.s[index[ind.nnull]] = zdat.ash$PosteriorSD^2
                x.w.v.s[index[!ind.nnull]] = 0
            }
        }
        mu.est = AvBasis(convert(x.w))
        if (return.loglr == TRUE) {
            logLR = sum(logLR.scale) 
        }
        if (post.var == TRUE) {
            mv.wd = wd.var(rep(0, n), filter.number = basis$filter.number, family = basis$family, type = "station")
            mv.wd$D = x.w.v.s
            mu.est.var = AvBasis.var(convert.var(mv.wd))
            # mu.est.var=apply((rep(1,n)%o%x.w.v.s)*Wl$Wi^2,1,sum) #diagonal of Winv*D*Winv'
            # vywt=qr.solve(Wl$Wmod,diag(x.w.v.s)) mu.est.var=diag(qr.solve(Wl$Wmod,t(vywt)))
            # mu.est.var=diag(Wl$Wi%*%diag(x.w.v.s)%*%t(Wl$Wi))
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


#' @title var.smooth
#' @return 'var.est' if posterior variances are not computed, and a list with elements 'var.est' and 'var.est.var' otherwise
#' @keywords internal
var.smooth = function(data, data.var, x.var.ini, basis, v.basis, Wl, filter.number, family, post.var, prior, pointmass, 
    nullcheck, VB, mixsd, mixcompdist, gridmult, jash, weight, J, n, SGD) {
    wmean = matrix(0, J, n)
    wvar = matrix(0, J, n)
    if (basis[[1]] == "haar" | v.basis == FALSE) {
        vtable = cxxtitable(data.var)$sumtable
        vdtable = cxxtitable(data)$difftable
        for (j in 0:(J - 1)) {
            ind.nnull = (vtable[j + 2, ] != 0)
            zdat.ash = shrink.wc(vdtable[j + 2, ind.nnull], sqrt(vtable[j + 2, ind.nnull]), prior = prior, pointmass = pointmass, nullcheck = nullcheck, 
                VB = VB, mixsd = mixsd, mixcompdist = mixcompdist, gridmult = gridmult, jash = jash, df = min(50, 2^(j + 
                  1)), SGD = SGD)
            wmean[j + 1, ind.nnull] = zdat.ash$PosteriorMean/2
            wmean[j + 1, !ind.nnull] = 0
            if ((sum(is.na(wmean[j + 1, ])) > 0) & (SGD == TRUE)) {
                zdat.ash = shrink.wc(vdtable[j + 2, ind.nnull], sqrt(vtable[j + 2, ind.nnull]), prior = prior, pointmass = pointmass, 
                  nullcheck = nullcheck, VB = VB, mixsd = mixsd, mixcompdist = mixcompdist, gridmult = gridmult, jash = jash, 
                  df = min(50, 2^(j + 1)), SGD = FALSE)
                wmean[j + 1, ind.nnull] = zdat.ash$PosteriorMean/2
                wmean[j + 1, !ind.nnull] = 0
            }
            if (post.var == TRUE) {
                wvar[j + 1, ind.nnull] = zdat.ash$PosteriorSD^2/4
                wvar[j + 1, !ind.nnull] = 0
            }
        }
        wwmean = -wmean
        var.est = cxxreverse_gwave(weight * sum(data) + (1 - weight) * sum(x.var.ini), wmean, wwmean)
        if (post.var == TRUE) {
            wwvar = wvar
            var.est.var = cxxreverse_gvwave(0, wvar, wwvar)
        }
    } else {
        x.w = wd(data, filter.number = filter.number, family = family, type = "station")
        x.w.v = apply((rep(1, n * J) %o% data.var) * Wl$W2, 1, sum)  #diagonal of W*V*W'
        x.pm = rep(0, n)
        x.w.v.s = rep(0, n * J)
        for (j in 0:(J - 1)) {
            index = (((J - 1) - j) * n + 1):((J - j) * n)
            x.w.j = accessD(x.w, j)
            x.w.v.j = x.w.v[index]
            ind.nnull = (x.w.v.j != 0)
            zdat.ash = shrink.wc(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]), prior = prior, pointmass = pointmass, 
                nullcheck = nullcheck, VB = VB, mixsd = mixsd, mixcompdist = mixcompdist, gridmult = gridmult, jash = jash, 
                df = min(50, 2^(j + 1)), SGD = SGD)
            x.pm[ind.nnull] = zdat.ash$PosteriorMean
            x.pm[!ind.nnull] = 0
            if ((sum(is.na(x.pm)) > 0) & (SGD == TRUE)) {
                zdat.ash = shrink.wc(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]), prior = prior, pointmass = pointmass, 
                  nullcheck = nullcheck, VB = VB, mixsd = mixsd, mixcompdist = mixcompdist, gridmult = gridmult, jash = jash, 
                  df = min(50, 2^(j + 1)), SGD = FALSE)
                x.pm[ind.nnull] = zdat.ash$PosteriorMean
                x.pm[!ind.nnull] = 0
            }
            x.w = putD(x.w, j, x.pm)
            if (post.var == TRUE) {
                x.w.v.s[index[ind.nnull]] = zdat.ash$PosteriorSD^2
                x.w.v.s[index[!ind.nnull]] = 0
            }
        }
        var.est = AvBasis(convert(x.w))
        if (post.var == TRUE) {
            mv.wd = wd.var(rep(0, n), filter.number = basis$filter.number, family = basis$family, type = "station")
            mv.wd$D = x.w.v.s
            var.est.var = AvBasis.var(convert.var(mv.wd))
            # mu.est.var=apply((rep(1,n)%o%x.w.v.s)*Wl$Wi^2,1,sum) #diagonal of Winv*D*Winv'
            # vywt=qr.solve(Wl$Wmod,diag(x.w.v.s)) mu.est.var=diag(qr.solve(Wl$Wmod,t(vywt)))
            # mu.est.var=diag(Wl$Wi%*%diag(x.w.v.s)%*%t(Wl$Wi))
        }
    }
    if (post.var == TRUE) {
        return(list(var.est = var.est, var.est.var = var.est.var))
    } else {
        return(var.est)
    }
}


#' Set default \code{ash} parameters.
#' @keywords internal 
#' @param ashparam: a list of parameters to be passed to ash.
setAshParam.gaus <- function(ashparam) {
    #by default ashparam$df=NULL
    #by default ashparam$mixsd=NULL
    #by default ashparam$g=NULL
    if (!is.list(ashparam))
        stop("Error: invalid parameter 'ashparam'")
    ashparam.default = list(optmethod="mixEM", pointmass=TRUE,
                   prior="nullbiased", gridmult=2, control = list(maxiter=5000,trace=FALSE), 
                   mixcompdist="normal", VB = FALSE, nullweight=10, nonzeromode=FALSE, outputlevel=1, randomstart=FALSE,fixg=FALSE, model="EE")
    ashparam = modifyList(ashparam.default, ashparam)
    if (!is.null(ashparam[["g"]]))
        stop("Error: ash parameter 'g' can only be NULL; if you want to specify ash parameter 'g' use multiseq arguments 'fitted.g' and/or 'fitted.g.intercept'")
    
    if(!((is.null(ashparam[["mixsd"]]))|(is.numeric(ashparam[["mixsd"]]) & (length(ashparam[["mixsd"]])<2)))) stop("Error: invalid parameter 'mixsd', 'mixsd'  must be null or a numeric vector of length >=2")
    if(!((ashparam[["prior"]] == "nullbiased") | (ashparam[["prior"]] == "uniform") | is.numeric(ashparam[["prior"]]))) stop("Error: invalid parameter 'prior', 'prior' can be a number or 'nullbiased' or 'uniform'")
    return(ashparam)
}


#' @title Estimate the underlying mean function from noisy Gaussian data
#'
#' @description This function takes a data vector of length a power of 2 and performs signal denoising using wavelet decomposition and an adaptive shrinkage prior on the wavelet parameters. The data are assumed to be (mostly) independent and Gaussian, but not necessarily identically distributed. 
#' 
#' @details We assume that the data come from the model $Y_t = \mu_t + \epsilon_t$ for t=1,...,T, where $\mu_t$ is an underlying mean, assumed to be spatially
#' structured (or treated as points sampled from a smooth continous function), and $\epsilon_t \sim N(0, \sigma_t)$, and are independent. 
#' Smash provides estimates of $\mu_t$ and $\sigma_t^2$ (and their posterior variances if desired).
#' @param sigma: a vector of standard deviations. Can be provided if known or estimated beforehand.
#' @param v.est: bool, indicating if variance estimation should be performed instead.
#' @param joint: bool, indicating if results of mean and variance estimation should be returned together.
#' @param v.basis: bool, indicating if the same wavelet basis should be used for variance estimation as mean estimation. If false, defaults to Haar basis for variance estimation (this is much faster than other bases).
#' @param post.var: bool, indicating if the posterior variance should be returned for the mean and/or variance estiamtes.
#' @param filter.number, family: wavelet basis to be used, as in \code{wavethresh}
#' @param return.loglr: bool, indicating if a logLR should be returned.
#' @param jash: indicates if the prior from method JASH should be used. This will often provide slightly better variance estimates (especially for nonsmooth variance functions), at the cost of computational efficiency. Defaults to FALSE.
#' @param SGD: bool, indicating if stochastic gradient descent should be used in the EM. Only applicable if jash=TRUE.
#' @param weight: optional parameter used in estimating overall variance. Only works for Haar basis. Defaults to 0.5. Setting this to 1 might improve variance estimation slightly
#' @param min.var: The minimum positive value to be set if the variance estimates are negative.
#' @param ashparam: a list of parameters to be passed to \code{ash}; default values are set by function \code{\link{setAshParam}}.
#' @return \code{smash.gaus} returns the mean estimate by default, or the variance estimate if \code{v.est} is TRUE. However, if \code{joint} is TRUE, then a list with the following is returned:
#' \item{mu.res}{a list with the mean estimate and its posterior variance if \code{post.var} is TRUE, or a vector of mean estimates otherwise}
#' \item{var.res}{a list with the variance estimate and its posterior variance if \code{v.est} is TRUE and \code{post.var} is TRUE, or a vector of variance estimates if only \code{v.est} is TRUE}
#' @examples
#' n=2^10
#' t=1:n/n
#' spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
#' mu.s=spike.f(t)
#' #Gaussian case
#' mu.t=(1+mu.s)/5
#' plot(mu.t,type='l')
#' var.fn=(0.0001+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
#' plot(var.fn,type='l')
#' rsnr=sqrt(5)
#' sigma.t=sqrt(var.fn)/mean(sqrt(var.fn))*sd(mu.t)/rsnr^2
#' X.s=rnorm(n,mu.t,sigma.t)
#' mu.est<-smash.gaus(X.s)
#' plot(mu.t,type='l')
#' lines(mu.est$mu.est,col=2)
#'
#' @export
smash.gaus = function(x, sigma = NULL, v.est = FALSE, joint = FALSE, v.basis = FALSE, post.var = FALSE, filter.number = 1, 
    family = "DaubExPhase", return.loglr = FALSE, jash = FALSE, SGD = TRUE, weight = 0.5, min.var = 1e-08, ashparam = list()) {
    n = length(x)
    J = log2(n)
    if (!isTRUE(all.equal(J, trunc(J)))) {
        stop("Error: number of columns of x must be power of 2")
    }
    if (length(sigma) == 1) {
        sigma = rep(sigma, n)
    }
    ashparam = setAshParam.gaus(ashparam)
    if (v.est == TRUE) {
        weight = 1
    } else {
        weight = 0.5
    }
    
    if (filter.number == 1 & family == "DaubExPhase") {
        basis = "haar"
    } else {
        basis = list(family = family, filter.number = filter.number)
    }
    # if(post.var==TRUE&basis[[1]]!='haar'){stop('Error: posterior variances returned only with Haar basis')}
    if (post.var == TRUE & v.est == TRUE & jash == TRUE) {
        stop("Error: Posterior variances for variance estimate not returned for method JASH")
    }
    if (joint == TRUE) {
        v.est = TRUE
    }
    
    
    tsum = sum(x)
    x.w.d = cxxtitable(x)$difftable
    Wl = NULL
    if (basis[[1]] != "haar") {
        x.w.d = wd(x, filter.number = filter.number, family = family, type = "station")
        Wl = ndwt.mat(n, filter.number = filter.number, family = family)
    }
    
    if (is.null(sigma)) {
        var.est1.ini = (x - lshift(x))^2/2
        var.est2.ini = (rshift(x) - x)^2/2
        var.est.ini = (var.est1.ini + var.est2.ini)/2
        mu.est = mu.smooth(x.w.d, var.est.ini, basis, tsum, Wl, FALSE, FALSE, ashparam$prior, ashparam$pointmass, ashparam$nullcheck, 
            ashparam$VB, ashparam$mixsd, ashparam$mixcompdist, ashparam$gridmult, J, n)
        var.est = (x - mu.est)^2
        var.var.est = 2/3 * var.est^2
        var.est = var.smooth(var.est, var.var.est, var.est.ini, basis, v.basis, Wl, filter.number, family, FALSE, ashparam$prior, 
            ashparam$pointmass, ashparam$nullcheck, ashparam$VB, ashparam$mixsd, ashparam$mixcompdist, ashparam$gridmult, 
            jash, weight, J, n, SGD = SGD)
        var.est[var.est <= 0] = 1e-08
        sigma = sqrt(var.est)
    }
    mu.res = mu.smooth(x.w.d, sigma^2, basis, tsum, Wl, return.loglr, post.var, ashparam$prior, ashparam$pointmass, ashparam$nullcheck, 
        ashparam$VB, ashparam$mixsd, ashparam$mixcompdist, 64, J, n)
    
    if (v.est == FALSE) {
        return(mu.res)
    } else {
        if (return.loglr == FALSE & post.var == FALSE) {
            mu.est = mu.res
        } else {
            mu.est = mu.res$mu.est
        }
        var.est = (x - mu.est)^2
        var.var.est = 2/3 * var.est^2
        var.res = var.smooth(var.est, var.var.est, 0, basis, v.basis, Wl, filter.number, family, post.var, ashparam$prior, 
            ashparam$pointmass, ashparam$nullcheck, ashparam$VB, ashparam$mixsd, ashparam$mixcompdist, ashparam$gridmult, 
            jash, 1, J, n, SGD = SGD)
        if (post.var == FALSE) {
            var.res[var.res <= 0] = min.var
        } else {
            var.res$var.est[var.res$var.est <= 0] = min.var
        }
        if (joint == FALSE) {
            return(var.res = var.res)
        } else {
            return(list(mu.res = mu.res, var.res = var.res))
        }
    }
}





#' This function performs TI thresholding for the case where the errors are heteroskedastic.
#'
#' @param x: the data. Should be a vector of length a power of 2
#' @param sigma: the standard deviation function. Can be provided if known or estimated beforehand.
#' @param method: the method to estimate the variance function. Can be 'rmad' for running MAD as described in Gao (1997), or 'smash'.
#' @param filter.number, family: the wavelet basis to be used.
#' @param min.level: the primary resolution level.
#' @return returns a vector of mean estimates
#' @details The 'rmad' option effectively implements the procedure described in Gao (1997), while the 'smash' option first estimates the variance function using
#' package \code{smash} and then performs thresholding given this variance function.
#' @references Gao, Hong-Ye (1997) Wavelet shrinkage estimates for heteroscedastic regression models. MathSoft, Inc.
#' @examples
#' n=2^10
#' t=1:n/n
#' spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
#' mu.s=spike.f(t)
#' #Gaussian case
#' mu.t=(1+mu.s)/5
#' plot(mu.t,type='l')
#' var.fn=(0.0001+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
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
#' @export
ti.thresh = function(x, sigma = NULL, method = "smash", filter.number = 1, family = "DaubExPhase", min.level = 3, ashparam = list()) {
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
            sigma = sqrt(smash.gaus(x, v.est = TRUE, v.basis = TRUE, filter.number = filter.number, family = family, 
                ashparam = ashparam, weight = 1))
        } else if (method == "rmad") {
            x.w = wd(x, filter.number = filter.number, family = family, type = "station")
            win.size = round(n/10)
            odd.boo = (win.size%%2 == 1)
            win.size = win.size + (1 - odd.boo)
            sigma = runmad(accessD(x.w, J - 1), win.size, endrule = "func")
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
        x.w = wd(x, filter.number = filter.number, family = family, type = "station")
        Wl = ndwt.mat(n, filter.number = filter.number, family = family)
        x.w.v = apply((rep(1, n * J) %o% (sigma^2)) * Wl$W2, 1, sum)  #diagonal of W*V*W'
        x.w.t = threshold.var(x.w, x.w.v, lambda.thresh, levels = (min.level):(J - 1), type = "hard")
        mu.est = AvBasis(convert(x.w.t))
    }
    return(mu.est)
}




######################################################################################################################### 










#' reflects a vector if it has length a power of 2; otherwise extends the vector to have length a power of 2 and then reflects it
#' @param x an n-vector
#' @return an n-vector containing the indices of the original signal x 
reflect <- function(x) {
    n = length(x)
    J = log2(n)
    if ((J%%1) == 0) {
        # if J is an integer, i.e. n is a power of 2
        eval.parent(substitute(x <- c(x, x[n:1])))
        return(1:n)
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
        eval.parent(substitute(x <- c(x.ini, x.mir)))
        return((lnum + 1):(lnum + n))
    }
}





#' Compute posterior mean and var for log(p), log(q), log(p0/p1) and log(q0/q1)
#' This function returns posterior means and variances of log(p), log(q), log(p0/p1) and log(q0/q1) as lp, lq, lpratio and lqratio, respectively, where p
#' is the probability of going left and q=1-p.
#' @param log: indicates if mean estimation is to be performed on the log scale
#' @return a list with elements 'lp.mean', 'lp.var', 'lq.mean', 'lq.var' 
compute.res <- function(alpha, log) {
    if (log == TRUE) {
        lp = ff.moments(alpha$mean, alpha$var)
        lq = ff.moments(-alpha$mean, alpha$var)  #find mean and variance of log(q)
        return(list(lp.mean = lp$mean, lp.var = lp$var, lq.mean = lq$mean, lq.var = lq$var))
    } else {
        lp = ff.moments_exp(alpha$mean, alpha$var)
        lq = ff.moments_exp(-alpha$mean, alpha$var)
        return(list(lp.mean = lp$mean, lp.var = lp$meansq, lq.mean = lq$mean, lq.var = lq$meansq))
    }
}


#' For each resolution, performs shrinkage and returns the posterior means and variances in matrix form
getlist.res = function(res, j, n, zdat, log, shrink, prior, pointmass, nullcheck, gridmult, mixsd, VB, mixcompdist) {
    ind = ((j - 1) * n + 1):(j * n)
    if (shrink == TRUE) {
        # apply ash to vector of intercept estimates and SEs
        zdat.ash = suppressWarnings(ash(zdat[1, ind], zdat[2, ind], prior = prior, outputlevel=2, pointmass = pointmass, 
            gridmult = gridmult, mixsd = mixsd, VB = VB, mixcompdist = mixcompdist, df = NULL, control = list(trace = FALSE)))
        alpha.mv = list(mean = zdat.ash$PosteriorMean, var = zdat.ash$PosteriorSD^2)  #find mean and variance of alpha
    } else {
        alpha.mv = list(mean = fill.nas(zdat[1, ind]), var = fill.nas(zdat[2, ind])^2)  #find mean and variance of alpha   
    }
    res.j = compute.res(alpha.mv, log)
    res = rbindlist(list(res, res.j))
    return(res)
}

#' Reconstructs the signal in data space
#' @param ls: estimated total intensity
#' @param res: the matrix of wavelet proportions/probabilities
#' @param log: bool, indicating if signal reconstruction should be in log space
#' @param n: length of data
#' @param J: =log2(n)
#' 
#' return reconstructed signal in original data space
recons.mv = function(ls, res, log, n, J) {
    if (log == TRUE) {
        # reconstructs estimate from the 'wavelet' space on the log level
        est.mean = cxxreverse_pwave(log(ls), matrix(res$lp.mean, J, n, byrow = TRUE), matrix(res$lq.mean, J, n, byrow = TRUE))
        est.var = cxxreverse_pwave(0, matrix(res$lp.var, J, n, byrow = TRUE), matrix(res$lq.var, J, n, byrow = TRUE))
    } else {
        # reconstruction on non-log level
        est.mean = exp(cxxreverse_pwave(log(ls), log(matrix(res$lp.mean, J, n, byrow = TRUE)), log(matrix(res$lq.mean, 
            J, n, byrow = TRUE))))
        est.ms = exp(cxxreverse_pwave(2 * log(ls), log(matrix(res$lp.var, J, n, byrow = TRUE)), log(matrix(res$lq.var, 
            J, n, byrow = TRUE))))
        est.var = pmax(est.ms - est.mean^2, 0)
    }
    return(list(est.mean = est.mean, est.var = est.var))
}


setAshParam.poiss <- function(ashparam) {
    #by default ashparam$df=NULL
    #by default ashparam$mixsd=NULL
    #by default ashparam$g=NULL
    if (!is.list(ashparam))
        stop("Error: invalid parameter 'ashparam'")
    ashparam.default = list(optmethod="mixEM", pointmass=TRUE,
                   prior="nullbiased", gridmult=2, control = list(maxiter=5000,trace=FALSE), 
                   mixcompdist="normal", VB = FALSE, nullweight=10, nonzeromode=FALSE, outputlevel=1, randomstart=FALSE,fixg=FALSE, model="EE")
    ashparam = modifyList(ashparam.default, ashparam)
    if (!is.null(ashparam[["g"]]))
        stop("Error: ash parameter 'g' can only be NULL; if you want to specify ash parameter 'g' use multiseq arguments 'fitted.g' and/or 'fitted.g.intercept'")
    
    if(!((is.null(ashparam[["mixsd"]]))|(is.numeric(ashparam[["mixsd"]]) & (length(ashparam[["mixsd"]])<2)))) stop("Error: invalid parameter 'mixsd', 'mixsd'  must be null or a numeric vector of length >=2")
    if(!((ashparam[["prior"]] == "nullbiased") | (ashparam[["prior"]] == "uniform") | is.numeric(ashparam[["prior"]]))) stop("Error: invalid parameter 'prior', 'prior' can be a number or 'nullbiased' or 'uniform'")
    return(ashparam)
}



#' @title Estimate the underlying intensity for Poisson counts.
#'
#' @description Main smoothing procedure for Poisson data. Takes a univariate inhomogeneous Poisson process and estimates its mean intensity.
#' 
#' @details We assume that the data come from the model $Y_t \sim Pois(\lambda_t)$ for t=1,...,T, where $\lambda_t$ is the underlying intensity, assumed to be spatially
#' structured (or treated as points sampled from a smooth continous function). The $Y_t$ are assumed to be independent. 
#' Smash provides estimates of $\lambda_t$ (and its posterior variance if desired).
#' 
#' @param x: a vector of Poisson counts of a length a power of 2
#' @param post.var: bool, indicates if the posterior variance should be returned
#' @param reflect: bool, indicates if the signals should be reflected; otherwise periodicity is assumed
#' @param lev: integer from 0 to J-1, indicating primary level of resolution. Should NOT be used (ie shrinkage is performed at all resolutions) unless there is good reason to do so.
#' @param log: bool, determines if smoothed signal is returned on log scale or not
#' @param pseudocounts: bool, a number to be added to counts. Passed to \code{glm.approx}
#' @param all: bool, indicates if pseudocounts should be added too all entries or only cases when either number of successes or number of failures (but not both) is 0. Passed to \code{glm.approx}
#' @param lm.approx: bool, indicates if a WLS shuold be fitted instead of GLM. Passed to \code{glm.approx}
#' @param cxx: bool, indicates if C++ code should be used to create TI tables.
#' @param ashparam: a list of parameters to be passed to \code{ash}; default values are set by function \code{\link{setAshParam}}.
#' @return \code{smash.poiss} returns the mean estimate by default, or the variance estimate if \code{post.var} is TRUE
#' @examples
#' n=2^10
#' t=1:n/n
#' spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
#' mu.s=spike.f(t)
#' mu.t=0.01+mu.s
#' X.s=rpois(n,mu.t)
#' mu.est=smash.poiss(X.s)
#' plot(mu.t,type='l')
#' lines(mu.est,col=2)
#'
#' @export
smash.poiss = function(x, post.var = FALSE, reflect = FALSE, lev = 0, log = FALSE, pseudocounts = 0.5, all = FALSE, 
    lm.approx = FALSE, bound = 0.02, cxx = TRUE, ashparam = list()) {
    if (is.matrix(x)) {
        if (nrow(x) == 1) {
            # change matrix x to vector
            x = as.vector(x)
        } else {
            stop("x cannot have multiple rows")
        }
    }
    
    ashparam = setAshParam.poiss(ashparam)
    
    if (!is.numeric(x)) 
        stop("Error: invalid parameter 'x': 'x' must be numeric")
    if (!is.logical(reflect)) 
        stop("Error: invalid parameter 'reflect', 'reflect' must be bool")
    if (!(is.numeric(pseudocounts) & pseudocounts > 0)) 
        stop("Error: invalid parameter 'pseudocounts', 'pseudocounts' must be a positive number")
    if (!is.logical(all)) 
        stop("Error: invalid parameter 'all', 'all'  must be bool")
    if (!((is.null(ashparam$mixsd)) | (is.numeric(ashparam$mixsd) & (length(ashparam$mixsd) < 2)))) 
        stop("Error: invalid parameter 'mixsd', 'mixsd'  must be null or a numeric vector of length >=2")
    if (!((ashparam$prior == "nullbiased") | (ashparam$prior == "uniform") | is.numeric(ashparam$prior))) 
        stop("Error: invalid parameter 'prior', 'prior' can be a number or 'nullbiased' or 'uniform'")
    
    J = log2(length(x))
    if ((J%%1) != 0) 
        {
            reflect = TRUE
        }  #if ncol(x) is not a power of 2, reflect x
    if (reflect == TRUE) 
        reflect.indices = reflect(x)  #reflect signal; this function is pseudo calling x by reference
    
    n = length(x)
    J = log2(n)
    
    # create the parent TI table for each signal, and put into rows of matrix y
    ls = sum(x)
    
    if (!cxx) {
        y = as.vector(t(ParentTItable(x)$parent))
    } else {
        y = cxxSParentTItable(x)
    }
    
    zdat = glm.approx(y, g = NULL, minobs = 1, pseudocounts = pseudocounts, center = FALSE, all = all, forcebin = TRUE, 
        repara = TRUE, lm.approx = lm.approx, disp = "add", bound = bound)
    # define empty matrices for posterior means and variances of log(p) and log(q)
    
    res = list()
    # loop through resolutions, smoothing each resolution separately
    for (j in 1:(J - lev)) {
        res = getlist.res(res, j, n, zdat, log, TRUE, ashparam$prior, ashparam$pointmass, ashparam$nullcheck, ashparam$gridmult, 
            ashparam$mixsd, ashparam$VB, ashparam$mixcompdist)
    }
    # Do not smooth for coarser levels, with everything the same as above but using the estimate and its variance as
    # the posterior mean and variance ie flat prior
    if (lev != 0) {
        for (j in (J - lev + 1):J) {
            res = getlist.res(res, j, n, zdat, log, FALSE, ashparam$prior, ashparam$pointmass, ashparam$nullcheck, 
                ashparam$gridmult, ashparam$mixsd, ashparam$VB, ashparam$mixcompdist)
        }
    }
    recons = recons.mv(ls, res, log, n, J)
    if (reflect == TRUE) {
        recons$est.mean = recons$est.mean[reflect.indices]
        recons$est.var = recons$est.var[reflect.indices]
    }
    if (post.var == FALSE) {
        # if post.var=TRUE then only estimate (and not variance) is returned
        return(est = recons$est.mean)
    } else {
        return(list(est = recons$est.mean, var = recons$est.var))
    }
}

 
