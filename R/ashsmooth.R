# Disclaimer: Functions 'interleave', 'lsift', 'rshift', 'simplify', 'setAshParam.gaus' and 'setAshParam.pois' in
# this file are adapted from the same functions defined in package 'Multiseq' under the GPL license, authored by
# Ester Pantaleo, Heejung Shim, Matthew Stephens and Zhengrong Xing. All TI table (eg 'titable', 'ParentTItable')
# and reconstruction (eg 'reverse.pwave', 'reverse.gwave') functions and their respective Rcpp counterparts are
# adapted and modified from functions 'ParentTItable' and 'reverse.pwave' in package 'Multiseq' under the GPL
# license, authored by Ester Pantaleo, Heejung Shim, Matthew Stephens and Zhengrong Xing, which are in turn ported
# into R from Matlab functions 'BMSMShrink' and 'TISumProd' as part of the BMSM project maintained by Eric Kolaczyk
# The 'wd' ('wd.D') and 'threshold' ('threshold.haar', 'threshold.var') class of functions are adapted and modified
# from 'wd' and 'threshold.wd' respectively in package 'wavethresh' under th GPL license, authored buy Guy Nason.


#' Interleave two vectors
#' @param x,y: two vectors of the same length
#' @return a vector of length twice that of x (or y)
interleave = function(x, y) {
    return(as.vector(rbind(x, y)))
}


#' Shift a vector one unit to the right
#' @param x: a vector
#' @return a vector of the same length as that of x 
rshift = function(x) {
    L = length(x)
    return(c(x[L], x[-L]))
}
lshift = function(x) {
    return(c(x[-1], x[1]))
}


#' Produces two TI tables. One table contains the difference between adjacent pairs of data in the same resolution, and the other table contains the sum.
#' @param sig: a signal of length a power of 2
#' @return a list of two TI tables in the form of matrices
titable = function(sig) {
    n = length(sig)
    J = log2(n)
    
    dmat = matrix(0, nrow = J + 1, ncol = n)
    ddmat = matrix(0, nrow = J + 1, ncol = n)
    dmat[1, ] = sig
    ddmat[1, ] = sig
    
    # dmat[1,] = as.matrix(sig) dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table
    
    for (D in 0:(J - 1)) {
        nD = 2^(J - D)
        nDo2 = nD/2
        twonD = 2 * nD
        for (l in 0:(2^D - 1)) {
            ind = (l * nD + 1):((l + 1) * nD)
            # ind2 = (l*twonD+1):((l+1)*twonD)
            x = dmat[D + 1, ind]
            lsumx = x[seq(from = 1, to = nD - 1, by = 2)] + x[seq(from = 2, to = nD, by = 2)]
            ldiffx = x[seq(from = 1, to = nD - 1, by = 2)] - x[seq(from = 2, to = nD, by = 2)]
            rx = rshift(x)
            rsumx = rx[seq(from = 1, to = nD - 1, by = 2)] + rx[seq(from = 2, to = nD, by = 2)]
            rdiffx = rx[seq(from = 1, to = nD - 1, by = 2)] - rx[seq(from = 2, to = nD, by = 2)]
            dmat[D + 2, ind] = c(lsumx, rsumx)
            ddmat[D + 2, ind] = c(ldiffx, rdiffx)
        }
    }
    return(list(sumtable = dmat, difftable = ddmat))
}




#' Produces a TI table containing the log difference between adjacent pairs of data in the same resolution.
#' @param sig: a signal of length a power of 2
#' @return a TI tables in the form of a matrix
tirtable = function(sig) {
    n = length(sig)
    J = log2(n)
    
    dmat = matrix(0, nrow = J + 1, ncol = n)
    ddmat = matrix(0, nrow = J + 1, ncol = n)
    dmat[1, ] = sig
    ddmat[1, ] = sig
    
    # dmat[1,] = as.matrix(sig) dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table
    
    for (D in 0:(J - 1)) {
        nD = 2^(J - D)
        nDo2 = nD/2
        twonD = 2 * nD
        for (l in 0:(2^D - 1)) {
            ind = (l * nD + 1):((l + 1) * nD)
            # ind2 = (l*twonD+1):((l+1)*twonD)
            x = dmat[D + 1, ind]
            lsumx = x[seq(from = 1, to = nD - 1, by = 2)] + x[seq(from = 2, to = nD, by = 2)]
            ldiffx = log(x[seq(from = 1, to = nD - 1, by = 2)]) - log(x[seq(from = 2, to = nD, by = 2)])
            rx = rshift(x)
            rsumx = rx[seq(from = 1, to = nD - 1, by = 2)] + rx[seq(from = 2, to = nD, by = 2)]
            rdiffx = log(rx[seq(from = 1, to = nD - 1, by = 2)]) - log(rx[seq(from = 2, to = nD, by = 2)])
            dmat[D + 2, ind] = c(lsumx, rsumx)
            ddmat[D + 2, ind] = c(ldiffx, rdiffx)
        }
    }
    return(ddmat)
}


#' Turns a TItable into a matrix that one can apply glm to (for Poisson denoising).
#' @param dmat: a k by n matrix
#' return a 2 by nk/2 matrix
simplify = function(dmat) {
    matrix(t(dmat), nrow = 2)
}

#' Reverse wavelet transform a set of wavelet coefficients in TItable format for Gaussian data.
#' @param lp: a J by n matrix of estimated wavelet coefficients. 
#' @param lq: a J by n matrix of complementary wavelet coefficients.
#' @param est: an n-vector. Usually a constant vector with each element equal to the estimated total mean.
#' @return reconstructed signal in the original data space.
reverse.gwave = function(est, lp, lq = NULL) {
    if (is.null(lq)) {
        lq = -lp
    }
    if (length(est) == 1) {
        est = rep(est, ncol(lp))
    }
    
    J = nrow(lp)
    
    for (D in J:1) {
        # print(exp(est)) readline('press a key')
        nD = 2^(J - D + 1)
        nDo2 = nD/2
        for (l in 0:(2^(D - 1) - 1)) {
            
            ind = (l * nD + 1):((l + 1) * nD)
            
            estvec = est[ind]/2
            lpvec = lp[D, ind]
            lqvec = lq[D, ind]
            
            estl = estvec[1:nDo2]
            lpl = lpvec[1:nDo2]
            lql = lqvec[1:nDo2]
            nestl = interleave(estl + lpl, estl + lql)  #interleaves the two
            
            
            estr = estvec[(nDo2 + 1):nD]
            lpr = lpvec[(nDo2 + 1):nD]
            lqr = lqvec[(nDo2 + 1):nD]
            nestr = interleave(estr + lpr, estr + lqr)  #interleaves the two
            nestr = lshift(nestr)
            
            est[ind] = 0.5 * (nestl + nestr)
        }
    }
    return(est)
}




#' Reverse wavelet transform a set of posterior variances for wavelet coefficients in TItable format, for Gaussian data.
#' @param lp: a J by n matrix of estimated variances. 
#' @param lq: a J by n matrix of complementary variances (=lp).
#' @param est: an n-vector. Usually 0.
#' @return reconstructed posterior variance in the original data space.
reverse.gvwave = function(est, lp, lq = NULL) {
    if (is.null(lq)) {
        lq = -lp
    }
    if (length(est) == 1) {
        est = rep(est, ncol(lp))
    }
    
    J = nrow(lp)
    
    for (D in J:1) {
        # print(exp(est)) readline('press a key')
        nD = 2^(J - D + 1)
        nDo2 = nD/2
        for (l in 0:(2^(D - 1) - 1)) {
            
            ind = (l * nD + 1):((l + 1) * nD)
            
            estvec = est[ind]/4
            lpvec = lp[D, ind]
            lqvec = lq[D, ind]
            
            estl = estvec[1:nDo2]
            lpl = lpvec[1:nDo2]
            lql = lqvec[1:nDo2]
            nestl = interleave(estl + lpl, estl + lql)  #interleaves the two
            
            estr = estvec[(nDo2 + 1):nD]
            lpr = lpvec[(nDo2 + 1):nD]
            lqr = lqvec[(nDo2 + 1):nD]
            nestr = interleave(estr + lpr, estr + lqr)  #interleaves the two
            nestr = lshift(nestr)
            
            est[ind] = 0.5 * (nestl + nestr)
        }
    }
    return(est)
}









#' Modified wd() function from package 'wavethresh' to only return the detail coefficients
#' @param data,filter.number,family,type,bc,berbose,min.scale,precond: as in wd()
#' @return the detail coefficients of a standard wavelet decomposition
wd.D = function(data, filter.number = 10, family = "DaubLeAsymm", type = "wavelet", bc = "periodic", verbose = FALSE, 
    min.scale = 0, precond = TRUE) {
    if (verbose == TRUE) 
        cat("wd: Argument checking...")
    if (!is.atomic(data)) 
        stop("Data is not atomic")
    DataLength <- length(data)
    nlevels <- nlevelsWT(data)
    if (is.na(nlevels)) 
        stop("Data length is not power of two")
    if (type != "wavelet" && type != "station") 
        stop("Unknown type of wavelet decomposition")
    if (type == "station" && bc != "periodic") 
        stop("Can only do periodic boundary conditions with station")
    if (verbose == TRUE) 
        cat("...done\nFilter...")
    if (bc != "interval") 
        filter <- filter.select(filter.number = filter.number, family = family)
    if (verbose == TRUE) 
        cat("...selected\nFirst/last database...")
    fl.dbase <- first.last(LengthH = length(filter$H), DataLength = DataLength, type = type, bc = bc)
    if (bc == "interval") {
        ans <- wd.int(data = data, preferred.filter.number = filter.number, min.scale = min.scale, precond = precond)
        fl.dbase <- first.last(LengthH = length(filter$H), DataLength = DataLength, type = type, bc = bc, current.scale = min.scale)
        filter <- list(name = paste("CDV", filter.number, sep = ""), family = "CDV", filter.number = filter.number)
        l <- list(transformed.vector = ans$transformed.vector, current.scale = ans$current.scale, filters.used = ans$filters.used, 
            preconditioned = ans$preconditioned, date = ans$date, nlevels = IsPowerOfTwo(length(ans$transformed.vector)), 
            fl.dbase = fl.dbase, type = type, bc = bc, filter = filter)
        class(l) <- "wd"
        return(l)
    }
    dtsp <- tsp(data)
    C <- rep(0, fl.dbase$ntotal)
    C[1:DataLength] <- data
    if (verbose == TRUE) 
        error <- 1 else error <- 0
    if (verbose == TRUE) 
        cat("built\n")
    if (verbose == TRUE) 
        cat("Decomposing...\n")
    nbc <- switch(bc, periodic = 1, symmetric = 2)
    if (is.null(nbc)) 
        stop("Unknown boundary condition")
    ntype <- switch(type, wavelet = 1, station = 2)
    if (is.null(filter$G)) {
        wavelet.decomposition <- .C("wavedecomp", C = as.double(C), D = as.double(rep(0, fl.dbase$ntotal.d)), H = as.double(filter$H), 
            LengthH = as.integer(length(filter$H)), nlevels = as.integer(nlevels), firstC = as.integer(fl.dbase$first.last.c[, 
                1]), lastC = as.integer(fl.dbase$first.last.c[, 2]), offsetC = as.integer(fl.dbase$first.last.c[, 3]), 
            firstD = as.integer(fl.dbase$first.last.d[, 1]), lastD = as.integer(fl.dbase$first.last.d[, 2]), offsetD = as.integer(fl.dbase$first.last.d[, 
                3]), ntype = as.integer(ntype), nbc = as.integer(nbc), error = as.integer(error), PACKAGE = "wavethresh")
    } else {
        wavelet.decomposition <- .C("comwd", CR = as.double(Re(C)), CI = as.double(Im(C)), LengthC = as.integer(fl.dbase$ntotal), 
            DR = as.double(rep(0, fl.dbase$ntotal.d)), DI = as.double(rep(0, fl.dbase$ntotal.d)), LengthD = as.integer(fl.dbase$ntotal.d), 
            HR = as.double(Re(filter$H)), HI = as.double(-Im(filter$H)), GR = as.double(Re(filter$G)), GI = as.double(-Im(filter$G)), 
            LengthH = as.integer(length(filter$H)), nlevels = as.integer(nlevels), firstC = as.integer(fl.dbase$first.last.c[, 
                1]), lastC = as.integer(fl.dbase$first.last.c[, 2]), offsetC = as.integer(fl.dbase$first.last.c[, 3]), 
            firstD = as.integer(fl.dbase$first.last.d[, 1]), lastD = as.integer(fl.dbase$first.last.d[, 2]), offsetD = as.integer(fl.dbase$first.last.d[, 
                3]), ntype = as.integer(ntype), nbc = as.integer(nbc), error = as.integer(error), PACKAGE = "wavethresh")
    }
    if (verbose == TRUE) 
        cat("done\n")
    error <- wavelet.decomposition$error
    if (error != 0) {
        cat("Error ", error, " occured in wavedecomp\n")
        stop("Error")
    }
    if (is.null(filter$G)) {
        l <- list(C = wavelet.decomposition$C, D = wavelet.decomposition$D, nlevels = nlevelsWT(wavelet.decomposition), 
            fl.dbase = fl.dbase, filter = filter, type = type, bc = bc, date = date())
    } else {
        l <- list(C = complex(real = wavelet.decomposition$CR, imaginary = wavelet.decomposition$CI), D = complex(real = wavelet.decomposition$DR, 
            imaginary = wavelet.decomposition$DI), nlevels = nlevelsWT(wavelet.decomposition), fl.dbase = fl.dbase, 
            filter = filter, type = type, bc = bc, date = date())
    }
    class(l) <- "wd"
    if (!is.null(dtsp)) 
        tsp(l) <- dtsp
    l$D
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
        zdat.ash = suppressWarnings(ash(wc, wc.var.sqrt, prior = prior, multiseqoutput = TRUE, pointmass = pointmass, 
            nullcheck = nullcheck, VB = VB, mixsd = mixsd, mixcompdist = mixcompdist, gridmult = gridmult, lambda1 = 1, 
            lambda2 = 0, df = NULL, control = list(trace = FALSE)))
    } else {
        zdat.ash = jasha(wc, wc.var.sqrt, df = df, SGD = SGD)
    }
    return(zdat.ash)
}


#' @title mu.smooth
#' @return 'mu.est' if posterior variances are not computed, and a list with elements 'mu.est' and 'mu.est.var' otherwise
#' @keywords internal
mu.smooth = function(wc, data.var, basis, tsum, Wl, post.var, prior, pointmass, nullcheck, VB, mixsd, mixcompdist, 
    gridmult, J, n) {
    wmean = matrix(0, J, n)
    wvar = matrix(0, J, n)
    if (basis[[1]] == "haar") {
        y = wc
        vtable = cxxtitable(data.var)$sumtable
        for (j in 0:(J - 1)) {
            ind.nnull = (vtable[j + 2, ] != 0)
            zdat.ash = shrink.wc(y[j + 2, ind.nnull], sqrt(vtable[j + 2, ind.nnull]), prior = prior, pointmass = pointmass, 
                nullcheck = nullcheck, VB = VB, mixsd = mixsd, mixcompdist = mixcompdist, gridmult = gridmult, jash = FALSE, 
                df = NULL, SGD = FALSE)
            wmean[j + 1, ind.nnull] = zdat.ash$PosteriorMean/2
            wmean[j + 1, !ind.nnull] = 0
            if (post.var == TRUE) {
                wvar[j + 1, ind.nnull] = zdat.ash$PosteriorSD^2/4
                wvar[j + 1, !ind.nnull] = 0
            }
        }
        wwmean = -wmean
        mu.est = cxxreverse_gwave(tsum, wmean, wwmean)
        if (post.var == TRUE) {
            wwvar = wvar
            mu.est.var = cxxreverse_gvwave(0, wvar, wwvar)
        }
    } else {
        x.w = wc
        x.w.v = apply((rep(1, n * J) %o% data.var) * Wl$W2, 1, sum)  #diagonal of W*V*W' 
        x.pm = rep(0, n)
        x.w.v.s = rep(0, n * J)
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
            if (post.var == TRUE) {
                x.w.v.s[index[ind.nnull]] = zdat.ash$PosteriorSD^2
                x.w.v.s[index[!ind.nnull]] = 0
            }
        }
        mu.est = AvBasis(convert(x.w))
        if (post.var == TRUE) {
            mv.wd = wd.var(rep(0, n), filter.number = basis$filter.number, family = basis$family, type = "station")
            mv.wd$D = x.w.v.s
            mu.est.var = AvBasis.var(convert.var(mv.wd))
            # mu.est.var=apply((rep(1,n)%o%x.w.v.s)*Wl$Wi^2,1,sum) #diagonal of Winv*D*Winv'
            # vywt=qr.solve(Wl$Wmod,diag(x.w.v.s)) mu.est.var=diag(qr.solve(Wl$Wmod,t(vywt)))
            # mu.est.var=diag(Wl$Wi%*%diag(x.w.v.s)%*%t(Wl$Wi))
        }
    }
    if (post.var == TRUE) {
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
            zdat.ash = shrink.wc(vdtable[j + 2, ], sqrt(vtable[j + 2, ]), prior = prior, pointmass = pointmass, nullcheck = nullcheck, 
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
    # by default ashparam$df=NULL by default ashparam$mixsd=NULL by default ashparam$g=NULL
    if (!is.list(ashparam)) 
        stop("Error: invalid parameter 'ashparam'")
    if (is.null(names(ashparam))) {
        ashparam[["mixsd"]] = NULL
    }
    if (is.null(ashparam[["pointmass"]])) 
        ashparam[["pointmass"]] = TRUE
    if (is.null(ashparam[["prior"]])) 
        ashparam[["prior"]] = "nullbiased"
    if (is.null(ashparam[["gridmult"]])) 
        ashparam[["gridmult"]] = 2
    if (is.null(ashparam[["VB"]])) 
        ashparam[["VB"]] = FALSE
    if (is.null(ashparam[["mixcompdist"]])) 
        ashparam[["mixcompdist"]] = "normal"
    if (is.null(ashparam[["nullcheck"]])) 
        ashparam[["nullcheck"]] = TRUE
    
    if (!((is.null(ashparam[["mixsd"]])) | (is.numeric(ashparam[["mixsd"]]) & (length(ashparam[["mixsd"]]) < 2)))) 
        stop("Error: invalid parameter 'mixsd', 'mixsd'  must be null or a numeric vector of length >=2")
    if (!((ashparam[["prior"]] == "nullbiased") | (ashparam[["prior"]] == "uniform") | is.numeric(ashparam[["prior"]]))) 
        stop("Error: invalid parameter 'prior', 'prior' can be a number or 'nullbiased' or 'uniform'")
    return(ashparam)
}


#' Estimate the underlying mean function from noisy Gaussian data
#'
#' This function takes a data vector of length a power of 2 and performs signal denoising using wavelet decomposition and an adaptive shrinkage prior on the wavelet parameters. The data are assumed to be (mostly) independent and Gaussian, but not necessarily identically distributed. 
#' @param sigma: a vector of standard deviations. Can be provided if known or estimated beforehand.
#' @param v.est: bool, indicating if variance estimation should be performed instead.
#' @param joint: bool, indicating if results of mean and variance estimation should be returned together.
#' @param v.basis: bool, indicating if the same wavelet basis should be used for variance estimation as mean estimation. If false, defaults to Haar basis for variance estimation (this is much faster than other bases).
#' @param post.var: bool, indicating if the posterior variance should be returned for the mean and/or variance estiamtes.
#' @param filter.number, family: wavelet basis to be used, as in \code{wavethresh}
#' @param jash: indicates if the prior from method JASH should be used. This will often provide slightly better variance estimates (especially for nonsmooth variance functions), at the cost of computational efficiency. Defaults to FALSE.
#' @param SGD: bool, indicating if stochastic gradient descent should be used in the EM. Only applicable if jash=TRUE.
#' @param weight: optional parameter used in estimating overall variance. Only works for Haar basis. Defaults to 0.5. Setting this to 1 might improve variance estimation slightly
#' @param min.var: The minimum positive value to be set if the variance estimates are negative.
#' @param ashparam: a list of parameters to be passed to \code{ash}; default values are set by function \code{\link{setAshParam}}.
#' @return \code{ashsmooth.gaus} returns the mean estimate by default, or the variance estimate if \code{v.est} is TRUE. However, if \code{joint} is TRUE, then a list with the following is returned:
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
#' mu.est<-ashsmooth.gaus(X.s)
#' plot(mu.t,type='l')
#' lines(mu.est,col=2)
#'
#' @export
ashsmooth.gaus = function(x, sigma = NULL, v.est = FALSE, joint = FALSE, v.basis = FALSE, post.var = FALSE, filter.number = 1, 
    family = "DaubExPhase", jash = FALSE, SGD = TRUE, weight = 0.5, min.var = 1e-08, ashparam = list()) {
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
        mu.est = mu.smooth(x.w.d, var.est.ini, basis, tsum, Wl, FALSE, ashparam$prior, ashparam$pointmass, ashparam$nullcheck, 
            ashparam$VB, ashparam$mixsd, ashparam$mixcompdist, ashparam$gridmult, J, n)
        var.est = (x - mu.est)^2
        var.var.est = 2/3 * var.est^2
        var.est = var.smooth(var.est, var.var.est, var.est.ini, basis, v.basis, Wl, filter.number, family, FALSE, ashparam$prior, 
            ashparam$pointmass, ashparam$nullcheck, ashparam$VB, ashparam$mixsd, ashparam$mixcompdist, ashparam$gridmult, 
            jash, weight, J, n, SGD = SGD)
        var.est[var.est <= 0] = 1e-08
        sigma = sqrt(var.est)
    }
    mu.res = mu.smooth(x.w.d, sigma^2, basis, tsum, Wl, post.var, ashparam$prior, ashparam$pointmass, ashparam$nullcheck, 
        ashparam$VB, ashparam$mixsd, ashparam$mixcompdist, 64, J, n)
    
    if (v.est == FALSE) {
        return(mu.res)
    } else {
        if (post.var == FALSE) {
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



#' This is a modified threshold.wd function from package 'wavethresh', designed to perform thresholding for heteroskedastic Gaussian errors when using the Haar basis
threshold.haar <- function(vdtable, vtable, lambda.thresh, levels, type = "hard", policy = "universal") {
    wmean <- vdtable[-1, ]/2
    d <- NULL
    n <- dim(vdtable)[2]
    nthresh <- length(levels)
    thresh <- list(0)
    for (i in 1:nthresh) {
        d <- vdtable[levels[i] + 2, ]
        noise.level <- sqrt(vtable[levels[i] + 2, ])
        nd <- length(d)
        if (lambda.thresh == 1) {
            lambda <- 2 * sqrt(log(nd))
        } else {
            lambda <- sqrt(2 * log(nd * log2(nd)))
        }
        thresh[[i]] <- lambda * noise.level
    }
    for (i in 1:nthresh) {
        d <- vdtable[levels[i] + 2, ]
        if (type == "hard") {
            d[abs(d) <= thresh[[i]]] <- 0
        } else if (type == "soft") {
            d <- (d * (abs(d) - thresh[[i]]) * (abs(d) > thresh[[i]]))/abs(d)
            d[is.na(d)] <- 0
        }
        wmean[levels[i] + 1, ] <- d/2
    }
    return(wmean)
}


#' This is a modified threshold.wd function from package 'wavethresh', designed to perform thresholding for heteroskedastic Gaussian errors when using non-Haar basis
threshold.var <- function(x.w, x.w.v, lambda.thresh, levels, type = "hard") {
    d <- NULL
    n <- 2^nlevelsWT(x.w)
    J <- nlevelsWT(x.w)
    nthresh <- length(levels)
    thresh <- list(0)
    for (i in 1:nthresh) {
        d <- accessD(x.w, level = levels[i])
        ind <- (((J - 1) - levels[i]) * n + 1):((J - levels[i]) * n)
        noise.level <- sqrt(x.w.v[ind])
        nd <- length(d)
        if (lambda.thresh == 1) {
            lambda <- 2 * sqrt(log(nd))
        } else {
            lambda <- sqrt(2 * log(nd * log2(nd)))
        }
        thresh[[i]] <- lambda * noise.level
    }
    for (i in 1:nthresh) {
        d <- accessD(x.w, level = levels[i])
        if (type == "hard") {
            d[abs(d) <= thresh[[i]]] <- 0
        } else if (type == "soft") {
            d <- (d * (abs(d) - thresh[[i]]) * (abs(d) > thresh[[i]]))/abs(d)
            d[is.na(d)] <- 0
        }
        x.w = putD(x.w, level = levels[i], v = d)
    }
    return(x.w)
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
ti.thresh = function(x, sigma = NULL, method = "smash", filter.number = 1, family = "DaubExPhase", min.level = 3) {
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
            sigma = sqrt(ashsmooth.gaus(x, v.est = TRUE, v.basis = TRUE, filter.number = filter.number, family = family, 
                weight = 1))
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







#' Produces two TI tables for Poisson counts. One is the standard TI table, and the other are the parent values used to make the TI table
#' @param sig: a signal of length a power of 2
#' @return a list of two TI tables in the form of matrices
ParentTItable = function(sig) {
    n = length(sig)
    J = log2(n)
    
    dmat = matrix(0, nrow = J + 1, ncol = n)
    dmat[1, ] = sig
    # dmat[1,] = as.matrix(sig)
    dmat2 = matrix(0, nrow = J, ncol = 2 * n)  #the parent table
    
    for (D in 0:(J - 1)) {
        nD = 2^(J - D)
        nDo2 = nD/2
        twonD = 2 * nD
        for (l in 0:(2^D - 1)) {
            ind = (l * nD + 1):((l + 1) * nD)
            ind2 = (l * twonD + 1):((l + 1) * twonD)
            x = dmat[D + 1, ind]
            lsumx = x[seq(from = 1, to = nD - 1, by = 2)] + x[seq(from = 2, to = nD, by = 2)]
            rx = rshift(x)
            rsumx = rx[seq(from = 1, to = nD - 1, by = 2)] + rx[seq(from = 2, to = nD, by = 2)]
            dmat[D + 2, ind] = c(lsumx, rsumx)
            dmat2[D + 1, ind2] = c(x, rx)
        }
    }
    return(list(TItable = dmat, parent = dmat2))
}





#' Reverse wavelet transform a set of probabilities in TItable format for Poisson data.
#' @param lp: a J by n matrix of probabilities. 
#' @param lq: a J by n matrix of complementary probabilities.
#' @param est: an n-vector. Usually a constant vector with each element equal to the estimated log(total rate).
#' @return reconstructed signal in the original data space.
reverse.pwave = function(est, lp, lq = NULL) {
    if (is.null(lq)) {
        lq = log(1 - exp(lp))
    }
    if (length(est) == 1) {
        est = rep(est, ncol(lp))
    }
    
    J = nrow(lp)
    
    for (D in J:1) {
        
        nD = 2^(J - D + 1)
        nDo2 = nD/2
        for (l in 0:(2^(D - 1) - 1)) {
            
            ind = (l * nD + 1):((l + 1) * nD)
            
            estvec = est[ind]
            lpvec = lp[D, ind]
            lqvec = lq[D, ind]
            
            estl = estvec[1:nDo2]
            lpl = lpvec[1:nDo2]
            lql = lqvec[1:nDo2]
            nestl = interleave(estl + lpl, estl + lql)  #interleaves the two
            
            estr = estvec[(nDo2 + 1):nD]
            lpr = lpvec[(nDo2 + 1):nD]
            lqr = lqvec[(nDo2 + 1):nD]
            nestr = interleave(estr + lpr, estr + lqr)  #interleaves the two
            nestr = lshift(nestr)
            
            est[ind] = 0.5 * (nestl + nestr)
        }
    }
    return(est)
}






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
        zdat.ash = suppressWarnings(ash(zdat[1, ind], zdat[2, ind], prior = prior, multiseqoutput = TRUE, pointmass = pointmass, 
            nullcheck = nullcheck, gridmult = gridmult, mixsd = mixsd, VB = VB, mixcompdist = mixcompdist, lambda1 = 1, 
            lambda2 = 0, df = NULL, control = list(trace = FALSE)))
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


setAshParam.pois <- function(ashparam) {
    # by default ashparam$df=NULL by default ashparam$mixsd=NULL by default ashparam$g=NULL
    if (!is.list(ashparam)) 
        stop("Error: invalid parameter 'ashparam'")
    if (is.null(names(ashparam))) {
        ashparam[["mixsd"]] = NULL
    }
    if (is.null(ashparam[["pointmass"]])) 
        ashparam[["pointmass"]] = TRUE
    if (is.null(ashparam[["prior"]])) 
        ashparam[["prior"]] = "nullbiased"
    if (is.null(ashparam[["gridmult"]])) 
        ashparam[["gridmult"]] = 0
    if (is.null(ashparam[["VB"]])) 
        ashparam[["VB"]] = FALSE
    if (is.null(ashparam[["mixcompdist"]])) 
        ashparam[["mixcompdist"]] = "normal"
    if (is.null(ashparam[["nullcheck"]])) 
        ashparam[["nullcheck"]] = TRUE
    
    if (!((is.null(ashparam[["mixsd"]])) | (is.numeric(ashparam[["mixsd"]]) & (length(ashparam[["mixsd"]]) < 2)))) 
        stop("Error: invalid parameter 'mixsd', 'mixsd'  must be null or a numeric vector of length >=2")
    if (!((ashparam[["prior"]] == "nullbiased") | (ashparam[["prior"]] == "uniform") | is.numeric(ashparam[["prior"]]))) 
        stop("Error: invalid parameter 'prior', 'prior' can be a number or 'nullbiased' or 'uniform'")
    return(ashparam)
}




#' Main smoothing procedure for Poisson data. Takes a univariate inhomogeneous Poisson process and estimates its mean intensity
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
#' @return \code{ashsmooth.pois} returns the mean estimate by default, or the variance estimate if \code{post.var} is TRUE
#' @examples
#' n=2^10
#' t=1:n/n
#' spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
#' mu.s=spike.f(t)
#' mu.t=0.01+mu.s
#' X.s=rpois(n,mu.t)
#' mu.est=ashsmooth.pois(X.s)
#' plot(mu.t,type='l')
#' lines(mu.est,col=2)
#'
#' @export
ashsmooth.pois = function(x, post.var = FALSE, reflect = FALSE, lev = 0, log = FALSE, pseudocounts = 0.5, all = FALSE, 
    lm.approx = FALSE, bound = 0.02, cxx = TRUE, ashparam = list()) {
    if (is.matrix(x)) {
        if (nrow(x) == 1) {
            # change matrix x to vector
            x = as.vector(x)
        } else {
            stop("x cannot have multiple rows")
        }
    }
    
    ashparam = setAshParam.pois(ashparam)
    
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

 
