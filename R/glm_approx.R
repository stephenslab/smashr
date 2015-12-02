## Disclaimer: The functions defined in this file are adapted from the nearly identical functions defined in package
## 'Multiseq' under the GPL license, authored by Ester Pantaleo, Heejung Shim, Matthew Stephens and Zhengrong Xing


#' Compute approximately unbiased variance estimates for the estimators for logit(p) when n is small
#' @param n: number of trials
#' @param s: number of successes
#' @param f: number of failures
v3 = function(n, s, f) {
    return((n + 1)/n * (1/(s + 1) + 1/(f + 1)))
}

#' Compute approximately unbiased variance estimates for the estimators for logit(p) when n is small
#' @param n: number of trials
#' @param s: number of successes
#' @param f: number of failures
vs = function(n, s, f) {
    vv = v3(n, s, f)
    return(vv * (1 - 2/n + vv/2))
}

#' Compute approximately unbiased variance estimates for the estimators for logit(p) when n is small
#' @param n: number of trials
#' @param s: number of successes
#' @param f: number of failures
vss = function(n, s, f) {
    vv = v3(n, s, f)
    return(vs(n, s, f) - 1/2 * vv^2 * (vv - 4/n))
}


#' Modified glm function to return relevant outputs, not allowing for underdispersion
#' @param x: covariante
#' @param y: response
#' @param forcebin: see glm.approx
#' @param repara: see glm.approx
#' @param ...: other inputs to glm.fit
#' @return a vector of intercept and slope estimates and their SEs
safe.quasibinomial.glm.fit = function(x, y, forcebin = FALSE, repara = FALSE, ...) {
    if (forcebin) {
        z = glm.fit(x, y, family = binomial(), ...)
        p1 = 1L:z$rank
        covmat = chol2inv(z$qr[[1]][p1, p1, drop = FALSE])
        se = sqrt(diag(covmat))
        if (repara == TRUE) {
            if (length(covmat) <= 1) {
                covmubeta = NA
            } else {
                covmubeta = covmat[2, 1]
            }
            mbvar = covmubeta/se[2]^2
            z$coef[1] = z$coef[1] - z$coef[2] * mbvar
            se[1] = sqrt(se[1]^2 + mbvar^2 * se[2]^2 - 2 * mbvar * covmubeta)
        }
    } else {
        z = glm.fit(x, y, family = quasibinomial(), ...)
        p1 = 1L:z$rank
        covmat.i = chol2inv(z$qr[[1]][p1, p1, drop = FALSE])
        df = z$df.residual
        if (df == 0) {
            d = 1
        } else {
            d = sum((z$weights * z$residuals^2)[z$weights > 0])/df
            d = d * (d >= 1) + 1 * (d < 1)
        }
        covmat = covmat.i * d
        se = sqrt(diag(covmat))
        if (repara == TRUE) {
            if (length(covmat) <= 1) {
                covmubeta = NA
            } else {
                covmubeta = covmat[2, 1]
            }
            mbvar = covmubeta/se[2]^2
            z$coef[1] = z$coef[1] - z$coef[2] * mbvar
            se[1] = sqrt(se[1]^2 + mbvar^2 * se[2]^2 - 2 * mbvar * covmubeta)
        }
    }
    if (repara == FALSE) {
        return(c(z$coef[1], se[1], z$coef[2], se[2]))
    } else {
        return(c(z$coef[1], se[1], z$coef[2], se[2], mbvar))
    }
}



#' Returns estimates of intercept and slope as well as their SEs, given other input options. Called in glm.approx
#' @param x: a 2n by 1 vector, with first n observations giving number of successes, and next n giving number of failures in a series of binomial experiments.
#' @param g: covariate. Can be null, in which case only the intercept estimate and its SE is returned
#' @param minobs: see glm.approx
#' @param pseudocounts: see glm.approx
#' @param all: see glm.approx
#' @param forcebin: see glm.approx
#' @param repara: see glm.approx
#' @return a vector of intercept and slope estimates and their SEs
bintest = function(x, g, minobs = 1, pseudocounts = 0.5, all = FALSE, forcebin = FALSE, repara = FALSE) {
    xmat = matrix(x, ncol = 2)
    zerosum = (apply(xmat, 1, sum) == 0)
    if (sum(!zerosum) > (minobs - 1)) {
        # check for enough observations
        ind1 = (xmat[, 1] == 0)
        ind2 = (xmat[, 2] == 0)
        if (all == TRUE) {
            xmat = xmat[!zerosum, , drop = F] + pseudocounts
        } else {
            xmat[ind1, 1] = xmat[ind1, 1] + pseudocounts
            xmat[ind1, 2] = xmat[ind1, 2] + pseudocounts
            xmat[ind2, 1] = xmat[ind2, 1] + pseudocounts
            xmat[ind2, 2] = xmat[ind2, 2] + pseudocounts
            xmat = xmat[!zerosum, , drop = F]
        }
        # check if there is enough variance in g among informative individuals
        g = g[!zerosum]
        ng = sum(!zerosum)
        dm = matrix(c(rep(1, ng), g), ncol = 2)
        if (!is.na(sd(g)) & (sd(g) > 0.1)) {
            dm[, 2] = g
            return(safe.quasibinomial.glm.fit(dm, xmat, forcebin, repara = repara))
        } else {
            if (repara == FALSE) {
                return(c(safe.quasibinomial.glm.fit(dm, xmat, forcebin, repara = repara)[1:2], NA, NA))
            } else {
                return(c(safe.quasibinomial.glm.fit(dm, xmat, forcebin, repara = repara)[1:2], NA, NA, NA))
            }
        }
    } else {
        # not enough observations, so just return NAs
        if (repara == FALSE) {
            return(c(NA, NA, NA, NA))
        } else {
            return(c(NA, NA, NA, NA, NA))
        }
    }
}

#' @title extract.sf
#' @return a list with elements 'x.s', 'x.f'
extract.sf = function(x, n) {
    return(list(x.s = as.vector(t(x[, (1:(2 * n))%%2 == 1])), x.f = as.vector(t(x[, (1:(2 * n))%%2 == 0]))))
}

#' @title add.counts
#' @return a list with elements 'x.s', 'x.f'
add.counts = function(x.s, x.f, eps, pseudocounts, all, index1, index2, indexn = NULL) {
    if (pseudocounts == 0) {
        x.s[index1] = x.s[index1] + eps
        x.f[index2] = x.f[index2] + eps
    } else if (pseudocounts != 0 & all == TRUE) {
        x.s = x.s + pseudocounts
        x.f = x.f + pseudocounts
    } else {
        x.s[index1] = x.s[index1] + pseudocounts
        x.f[index1] = x.f[index1] + pseudocounts
        x.s[index2] = x.s[index2] + pseudocounts
        x.f[index2] = x.f[index2] + pseudocounts
    }
    if (!is.null(indexn)) {
        x.s[indexn] = 0
        x.f[indexn] = 0
    }
    return(list(x.s = x.s, x.f = x.f))
}

#' @title compute.approx.z
#' Compute a vector of logit(p) given a vector of successes and failures, as well as its variance estimates (MLE with approximation at endpoints for mean; a mix of Berkso's estimator and Tukey's estimator for variance)
#' @return a list with elements 'mu', 'var' and optionally 'p'
compute.approx.z = function(x.s, x.f, bound, eps, pseudocounts, all, indexn = NULL, return.p = FALSE) {
    # compute mu
    index1 = (x.s/x.f) <= bound  #find indices for which binomial success or failures are too small
    index2 = (x.f/x.s) <= bound
    index1[is.na(index1)] = FALSE
    index2[is.na(index2)] = FALSE  #this is the same as above!!!
    x = add.counts(x.s, x.f, eps, pseudocounts, all, index1, index2, indexn)  #add pseudocounts
    s = x$x.s + x$x.f
    mu = log(x$x.s/x$x.f)  #compute logit(p) to be used as observations
    mu[index1] = mu[index1] - 0.5  #end-point correction
    mu[index2] = mu[index2] + 0.5
    # compute var compute var(logit(p))
    if (all == FALSE) {
        var = vss(s, x$x.s, x$x.f)
        var[index1] = vss(s[index1] - 2 * pseudocounts, x$x.s[index1] - pseudocounts, x$x.f[index1] - pseudocounts)
        var[index2] = vss(s[index2] - 2 * pseudocounts, x$x.s[index2] - pseudocounts, x$x.f[index2] - pseudocounts)
    } else {
        var = vss(s - 2 * pseudocounts, x$x.s - pseudocounts, x$x.f - pseudocounts)
    }
    var[var == Inf] = 1e+20
    if (return.p == TRUE) 
        return(list(mu = mu, var = var, p = x$x.s/s)) else return(list(mu = mu, var = var))
}


#' Compute estimates and standard errors for mu and beta when fitting WLS, as well as the covariance between mu and beta
#' @return a list elements 'coef', 'se', 'mbvar'
wls.coef = function(z, disp, indexnm, n, ng, forcebin, g = NULL, repara = NULL) {
    # compute vector of dfs for all n linear models (disregarding obs with missing data)
    if (is.null(g)) 
        df = pmax(colSums(!indexnm) - 1, 0) else df = pmax(colSums(!indexnm) - 2, 0)
    zm = matrix(z$mu, ncol = n, byrow = T)  #create ng*n matrix of logit(p)
    zm[indexnm] = 0
    vm = matrix(z$var, ncol = n, byrow = T)  #create ng*n matrix of var(logit(p))
    res = wls.mb(zm, vm, disp, indexnm, ng, df, forcebin, g, n)
    if (disp == "add") {
        # return estimates if multiplicative dispersion is assumed
        vm[indexnm] = NA
        vv = pmax((res$wrse2 - 1) * colMeans(vm, na.rm = T), 0)  #computes crude estimate of sigma_u^2 as in documentation
        res = wls.mb(zm, vm, disp, indexnm, ng, df, forcebin, g, n, vv)
    }
    if (is.null(g)) 
        return(list(coef = res$coef, se = res$se, mbvar = NULL)) else {
        coef = c(res$muhat, res$betahat)
        se = c(res$semuhat, res$sebetahat)
        if (repara == TRUE) {
            # return reparametrized muhat and behat as well as their SEs, together with gamma as defined in documentation
            mbvar = res$covmubeta/res$sebetahat^2
            coef[1:n] = res$muhat - res$betahat * mbvar
            se[1:n] = sqrt(res$semuhat^2 + mbvar^2 * res$sebetahat^2 - 2 * mbvar * res$covmubeta)
        } else {
            if (repara == FALSE) 
                mbvar = NULL else stop("Error: invalid argument 'repara'")
        }
        return(list(coef = coef, se = se, mbvar = mbvar))
    }
}

#' wls.mb
#' @return a list with elements 'coef', 'se', 'wrse2' if g is not specified, or a list with elements 'muhat', 'semuhat', 'betahat', 'sebetahat', 'covmubeta', 'wrse2' otherwise
wls.mb = function(z, v, disp, indexnm, ng, df, forcebin, g = NULL, n = NULL, vv = NULL) {
    if (is.null(vv)) {
        # compute weights for each of the n models
        w = 1/v
    } else {
        w = 1/(v + rep(1, ng) %o% vv)
    }
    w[indexnm] = 0
    ws = colSums(w)  #define sum of weights for each of the n models (to be used later)
    if (is.null(g)) {
        muhat = colSums(w * z)/ws  #compute muhat for each of the n models using formula in documentation
        wrse = sqrt(colSums((z - rep(1, ng) %o% muhat)^2 * w)/df)  #compute residual standard error
    } else {
        gwmean = colSums(w * g)/ws  #define weighted center of g for each of the n models (to be used later)
        ggwmeanm = g %o% rep(1, n) - rep(1, ng) %o% gwmean  #define weighted difference between each g and its weighted center for each of the n models (to be used later)
        wgg = colSums(w * ggwmeanm^2)  #compute sum_j w_j^2*(g_j-gwmean)^2 (to be used later)
        wgg.ind = wgg < 1e-06
        wgg[wgg.ind] = 0
        betahat = colSums(w * z * ggwmeanm)/colSums(w * ggwmeanm^2)  #compute betahat using formula in documentation
        g.betahat = g %o% betahat
        # compute betahat*g for each of the n models
        muhat = colSums(w * (z - g.betahat))/ws  #compute muhat using formula in documentation
        wrse = sqrt(colSums((z - rep(1, ng) %o% muhat - g.betahat)^2 * w)/df)  #compute residual standard error
    }
    wrse[is.na(wrse)] = 1
    if (forcebin | (is.null(g) & (ng == 2)) | (!is.null(g) & (ng == length(unique(g))))) {
        # force dispersion to be absent (also in the case with only 1 observation in each group)
        wrse = 1
    } else {
        if (is.null(vv)) {
            wrse[(wrse == Inf) | (wrse < 1)] = 1  #do not allow for 'underdispersion'
        } else {
            wrse[(wrse == Inf) | (vv == 0)] = 1
        }
    }
    wrse2 = wrse^2
    if (is.null(g)) {
        semuhat = sqrt(wrse2/ws)
        return(list(coef = muhat, se = semuhat, wrse2 = wrse2))
    } else {
        sebetahat = sqrt(wrse2/wgg)  #compute se(betahat) using formula in documentation
        sebetahat[wgg.ind] = NA
        semuhat = sqrt((1/ws + gwmean^2/wgg) * wrse2)  #compute se(muhat) using formula in documentation
        semuhat[wgg.ind] = NA
        covmubeta = colSums(w * ggwmeanm)/ws/wgg * wrse2 - gwmean * sebetahat^2  #compute covariance between muhat and betahat
        return(list(muhat = muhat, semuhat = semuhat, betahat = betahat, sebetahat = sebetahat, covmubeta = covmubeta, 
            wrse2 = wrse2))
    }
}

#' Computes the dispersion parameter when fitting glm
#' @return a vector of dispersion parameters for each fitted model, or 1 if dispersion is absent
compute.dispersion = function(p, n, ng, indexnm, forcebin, ind = NULL, ord = NULL, lg = NULL, x = NULL, x.s = NULL, 
    x.f = NULL) {
    if (is.null(lg)) {
        if (forcebin | ng == 1) {
            # force dispersion to be absent
            return(1)
        }
    } else {
        if (forcebin | (ng == lg)) {
            # (or if there is 1 obs in each group)
            return(1)
        }
    }
    
    # find effective number of observations after getting rid of missing data
    ngn = !indexnm
    ngn = colSums(ngn)
    ngn = rep(ngn, times = ng)
    
    if (is.null(lg)) {
        ss = x.s + x.f
        p = rep(p, times = ng)
        d.ini = 1/(ngn - 1) * (x.s - p * ss)^2/(ss * p * (1 - p))  #compute dispersion factor as in McC and Nelder
        d.ini[d.ini == Inf] = 0
        d.ini[is.na(d.ini)] = 0
        d.m = matrix(d.ini, ncol = ng)
        d = rowSums(d.m)
        d[d < 1] = 1  #do not allow underdispersion
    } else {
        x = x[ord, ]
        x.sf = extract.sf(x, n)
        s = x.sf$x.s + x.sf$x.f
        pn = NULL
        for (i in 1:lg) {
            pn = c(pn, rep(p[(n * (i - 1) + 1):(n * i)], times = sum(ind[[i]])))
        }
        d.ini = 1/(ngn - lg) * (x.sf$x.s - pn * s)^2/(s * pn * (1 - pn))  #compute dispersion factor as in McC and Nelder
        d.ini[d.ini == Inf] = 0
        d.ini[is.na(d.ini)] = 0
        d.m = matrix(d.ini, ncol = ng)
        d = rowSums(d.m)
        d[d < 1] = 1  #do not allow underdispersion
        d = rep(d, times = lg)
    }
    return(d)
}




#' Compute estimates and standard errors for mu and beta when fitting WLS, as well as the covariance between mu and beta. 
#' @return a list elements 'coef', 'se' and optionally 'mbvar' if lg=2 and repara=TRUE, or 'covv' if lg=3
glm.coef = function(z, g, n, center, repara) {
    lg = length(levels(g))
    mbvar = NULL
    if (lg == 2) {
        # 2 categories
        covv = NULL
        if (center == TRUE) {
            # considered centered and uncentered covariate separately
            g.num = sort(as.numeric(levels(g))[g])
            g.num = unique(g.num - mean(g.num))
            w1 = g.num[1]  #weights come in because covariate is centered
            w2 = g.num[2]
            coef = w2 * z$mu[1:n] - w1 * z$mu[(n + 1):(2 * n)]  #compute logit(p) for each group
            coef = c(coef, z$mu[(n + 1):(2 * n)] - z$mu[1:n])  #compute intercept and slope
            var = w2^2 * z$var[1:n] + w1^2 * z$var[(n + 1):(2 * n)]  #compute var(logit(p)) for each group
            var = c(var, z$var[(n + 1):(2 * n)] + z$var[1:n])  #compute var for intercept and slope
        } else {
            coef = z$mu - c(rep(0, n), rep(z$mu[1:n], times = (lg - 1)))  #compute intercept and slope
            var = z$var + c(rep(0, n), rep(z$var[1:n], times = (lg - 1)))  #compute var of intercept and slope
        }
        if (repara == TRUE) {
            mbvar = -var[1:n]/var[(n + 1):(2 * n)]  #compute gamma as in documentation if reparametrization is used
            coef[1:n] = coef[1:n] - coef[(n + 1):(2 * n)] * mbvar  #reparametrized estimates
            var[1:n] = var[1:n] - var[1:n]^2/var[(n + 1):(2 * n)]  #reparametrized Ses
        }
    } else if (lg == 3) {
        # 3 groups case as in PoissonBinomial_etc considered centered and uncentered covariate separately
        if (center == TRUE) {
            g.num = sort(as.numeric(levels(g))[g])
            g.num[g.num != 1] = 0
            g.num[g.num == 1] = 1
            g.num = unique(g.num - mean(g.num))
            w1.1 = g.num[1]
            w1.2 = g.num[2]
            g.num = sort(as.numeric(levels(g))[g])
            g.num[g.num != 2] = 0
            g.num[g.num == 2] = 1
            g.num = unique(g.num - mean(g.num))
            w2.1 = g.num[1]
            w2.2 = g.num[2]
            coef = (w1.2 + w2.1) * z$mu[1:n] - w1.1 * z$mu[(n + 1):(2 * n)] - w2.1 * z$mu[(2 * n + 1):(3 * n)]
            coef = c(coef, z$mu[(n + 1):(3 * n)] - rep(z$mu[1:n], times = (lg - 1)))
            var = (w1.2 + w2.1)^2 * z$var[1:n] + (w1.1)^2 * z$var[(n + 1):(2 * n)] + (w2.1)^2 * z$var[(2 * n + 1):(3 * 
                n)]
            var = c(var, z$var[(n + 1):(3 * n)] + rep(z$var[1:n], times = (lg - 1)))
            covv = z$var[1:n]
        } else {
            # 3 groups case
            coef = z$mu - c(rep(0, n), z$mu[1:(2 * n)])
            var = z$var + c(rep(0, n), z$var[1:(2 * n)])
            covv = -z$var[(n + 1):(2 * n)]
        }
    }
    return(list(mu = coef, var = var, mbvar = mbvar, covv = covv))
}

#' @title compute.lm
#' @return a matrix with estimates for mu and beta, as well as their SEs. Optionally returns 'mbvar' if specified
compute.lm = function(g, coef, se, mbvar, n, index, repara) {
    if (is.null(g)) {
        na.ind = is.na(coef[1:n]) | is.na(se[1:n])
        coef[na.ind | index] = NA  #set muhat and se(muhat) to NA for those models with insufficiant data or NAs
        se[na.ind | index] = NA
        return(matrix(c(coef, se), nrow = 2, byrow = T))
    } else {
        index = rep(index, times = 2)
        na.ind = is.na(coef[1:n]) | is.na(se[1:n]) | is.na(coef[(n + 1):(2 * n)]) | is.na(se[(n + 1):(2 * n)])
        na.ind2 = rep(na.ind, times = 2)
        coef[na.ind2 | index] = NA  #set muhat and se(muhat) to NA for those models with insufficiant data or NAs
        se[na.ind2 | index] = NA
        toreturn = array(rbind(coef, se), dim = c(2, n, 2))
        if (repara == TRUE) {
            mbvar[index[1:n]] = NA
            mbvar[na.ind] = NA
            return(matrix(rbind(apply(toreturn, 2, rbind), mbvar), ncol = n))
        } else return(matrix(rbind(apply(toreturn, 2, rbind), mbvar), ncol = n))  #should this be return(apply(toreturn,2,rbind))????
    }
}

#' @title compute.glm
#' @return a matrix with estimates for mu and beta, as well as their SEs. Optionally returns 'mbvar' if specified
compute.glm = function(x, g, d, n, na.index, repara) {
    se = sqrt(x$var * d)  #dispersion
    mu = x$mu
    mu[na.index] = NA  #set muhat and se(muhat) to NA for those models with insufficiant data
    se[na.index] = NA
    
    if (is.null(g)) 
        return(matrix(c(mu, se), nrow = 2, byrow = T)) else {
        if (is.factor(g)) {
            lg = length(levels(g))
            toreturn = array(rbind(mu, se), dim = c(2, n, lg))
            if (lg == 3) {
                if (length(d) == 1) 
                  covv = x$covv * d  #dispersion
 else covv = x$covv * d[1:n]  #dispersion
                covv[na.index[1:n]] = NA  #check that this is correct?????
                return(matrix(rbind(apply(toreturn, 2, rbind), covv), ncol = n))
            } else if (lg == 2) {
                if (repara == FALSE) 
                  return(apply(toreturn, 2, rbind)) else {
                  mbvar = x$mbvar
                  mbvar[na.index[1:n]] = NA
                  return(matrix(rbind(apply(toreturn, 2, rbind), mbvar), ncol = n))
                }
            }
        }
    }
}



#' This function essentially fits the model specified in documentation, using either a weighted least squares approach or a generalized linear model approach, with some modifications. It is optimized for fitting multiple models simultaneously
#' @param x: a matrix of N (# of samples) by 2*T (T: # of WCs or, more precisely, of different scales and locations in multi-scale space); Two consecutive columns correspond to a particular scale and location; The first column (the second column) contains # of successes (# of failures) for each sample at corresponding scale and location.
#' @param g: a vector of covariates. Can be a factor (2 groups only) or quantitative
#' @param pseudocounts: adds pseudocounts in the case of low success of failure counts.
#' @param all: bool,  when success or failure counts are low, controls whether or not to add pseudocounts to all observations or only the low counts.
#' @param eps: adds tiny pseudocounts in the case of 0 successes or 0 failures. Only used when pseoducounts=0
#' @param center: bool, specifies if g should be centered or not
#' @param repara: bool, specifies if model should be reparametrized so that muhat and betahat are uncorrelated
#' @param forcebin: bool, if TRUE don't allow for overdipersion.
#' @param lm.approx: bool, specifies if WLS should be used (TRUE), or GLM (FALSE)
#' @param disp: 'all' or 'mult', indicates which type of overdispersion is assumed when lm.approx=TRUE
#' @param bound: numeric, indicates the threshold of the success vs failure ratio below which pseudocounts will be added
#'
#' @export
#' @return a matrix of 2 (or 5 if g is provided) by T (# of WCs); Each row contains alphahat (1st row), standard error of alphahat (2nd), betahat (3rd), standard error of betahat (4th), covariance between alphahat and betahat (5th) for each WC.
glm.approx = function(x, g = NULL, minobs = 1, pseudocounts = 0.5, all = FALSE, eps = 1e-08, center = FALSE, repara = FALSE, 
    forcebin = FALSE, lm.approx = FALSE, disp = c("add", "mult"), bound = 0.02) {
    disp = match.arg(disp)
    if (is.vector(x)) 
        {
            dim(x) <- c(1, length(x))
        }  #if x is a vector convert to matrix 
    n = ncol(x)/2
    ng = nrow(x)
    
    if (ng == 1) {
        lm.approx = FALSE  #if x has 1 row don't use the lm approximation
        forcebin = TRUE
    } else {
        x.sf = extract.sf(x, n)  #extract success and failure counts separately
        indexn = (x.sf$x.s == 0 & x.sf$x.f == 0)  #find indices for which there is no data
        indexnm = matrix(indexn, nrow = ng, byrow = T)
    }
    if (lm.approx == TRUE) {
        # use WLS approximation
        na.index = colSums(matrix((x.sf$x.s + x.sf$x.f) != 0, ncol = n, byrow = T)) < minobs  #find indices for which there is insufficient data
        z = compute.approx.z(x.sf$x.s, x.sf$x.f, bound, eps, pseudocounts, all, indexn)  #obtain estimates for logit(p) and var(logit(p))
        if (is.null(g)) {
            # smoothing multiple signals without covariate
            res = wls.coef(z, disp, indexnm, n, ng, forcebin)
            return(compute.lm(g, res$coef, res$se, res$mbvar, n, na.index, repara))
        } else {
            if (is.factor(g)) 
                g = as.numeric(levels(g))[g]  #if g is a 2-level factor convert to numeric
            if (center == TRUE) 
                g = g - mean(g)
            res = wls.coef(z, disp, indexnm, n, ng, forcebin, g, repara)
            return(compute.lm(g, res$coef, res$se, res$mbvar, n, na.index, repara))
        }
    } else {
        # use GLM case where g is absent smoothing multiple signals
        if (is.null(g)) {
            if (ng > 1) {
                x = colSums(x)  #pool data together as in GLM
            }
            x = matrix(x, ncol = n)
            na.index = (colSums(x) == 0)
            z = compute.approx.z(x[1, ], x[2, ], bound, eps, pseudocounts, all, NULL, TRUE)  #obtain estimates for logit(p) and var(logit(p)) #indexn=NULL???
            d = compute.dispersion(z$p, n, ng, indexnm, forcebin, x.s = x.sf$x.s, x.f = x.sf$x.f)  #computes dispersion  
            return(compute.glm(z, g, d, n, na.index, repara))
        } else {
            # case where g is present first consider case where g is factor, and hence with closed form solution
            if (is.factor(g)) {
                lg = length(levels(g))
                ord = sort.list(g)
                if (ng > lg) {
                  # pool successes and failures in the same category, depending on if there are more obs than no. of categories or
                  # not
                  x.mer = matrix(0, nrow = lg, ncol = 2 * n)
                  ind = list(0)
                  for (i in 1:lg) {
                    ind[[i]] = (g == levels(g)[i])
                    x.mer[i, ] = colSums(matrix(x[ind[[i]], ], nrow = sum(ind[[i]])))
                  }
                } else {
                  ind = NULL
                  x.mer = x[ord, ]
                }
                x.mer.s = x.mer[, (1:(2 * n))%%2 == 1]  #now consider pooled data as effective raw data and extract successes and failures
                x.mer.f = x.mer[, (1:(2 * n))%%2 == 0]
                na.index = colSums(matrix((x.mer.s + x.mer.f) != 0, ncol = n)) < pmin(minobs, lg)  #find indices with insufficient data #since data are pooled, min no. of obs cannot be less than number of groups in g
                na.index = rep(na.index, times = lg)
                x.s.m = as.vector(t(x.mer.s))  #extract successes and failures from pooled data
                x.f.m = as.vector(t(x.mer.f))
                indexn = (x.s.m == 0 & x.f.m == 0)  #define indices where pooled data still has no data
                z = compute.approx.z(x.s.m, x.f.m, bound, eps, pseudocounts, all, indexn, TRUE)  #obtain estimates for logit(p) and var(logit(p)) #indexn?
                res = glm.coef(z, g, n, center, repara)
                d = compute.dispersion(z$p, n, ng, indexnm, forcebin, ind, ord, lg, x)  #computes dispersion
                return(compute.glm(res, g, d, n, na.index, repara))
            } else {
                # now consider the case when g is quantitative
                x = matrix(x, ncol = n)
                if (center == TRUE) 
                  g = g - mean(g)
                # use the bintest function to fit a GLM separately to each case for quantitative covariate
                return(apply(x, 2, bintest, g = g, minobs = minobs, pseudocounts = pseudocounts, all = all, forcebin = forcebin, 
                  repara = repara))
            }
        }
    }
} 
