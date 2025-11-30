# Disclaimer #1: The functions defined in this file are adapted from
# their counterparts defined in the wavethresh package under the GPL
# license, developed by Guy Nason.

# Disclaimer #2: Some of the functions defined here interface to C
# code from the wavethresh package. This is not Best Practice as the C
# code is not considered stable; it could very easily change in future
# versions of the package.

# Modified wd() function from package 'wavethresh' to only return the
# detail coefficients. It returns the detail coefficients of a
# standard wavelet decomposition.
#
#' @importFrom wavethresh wd
wd.D = function (data, filter.number = 10, family = "DaubLeAsymm",
                 type = "wavelet", bc = "periodic", verbose = FALSE, 
                 min.scale = 0, precond = TRUE) {
  l = wavethresh::wd(data = data, filter.number = filter.number,
                    family = family, type = type, bc = bc,
                    verbose = verbose, min.scale = min.scale,
                    precond = precond)
  return(l$D)
}

# This is a modified threshold.wd function from package 'wavethresh',
# designed to perform thresholding for heteroskedastic Gaussian errors
# when using the Haar basis
threshold.haar <- function (vdtable, vtable, lambda.thresh, levels,
                            type = "hard", policy = "universal") {
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

# This is a modified threshold.wd function from package 'wavethresh',
# designed to perform thresholding for heteroskedastic Gaussian errors
# when using non-Haar basis.
threshold.var <- function (x.w, x.w.v, lambda.thresh, levels, type = "hard") {
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

# This function performs "decomposition" of variances of detail
# coefficients for a given wavelet basis in wavelet transformation.
#
#' @importFrom stats tsp
#' @importFrom stats tsp<-
#' @importFrom wavethresh filter.select
#' @importFrom wavethresh first.last
#' @importFrom wavethresh wd.int
#' @importFrom wavethresh IsPowerOfTwo
wd.var <- function (data, filter.number = 10, family = "DaubLeAsymm",
                    type = "wavelet", bc = "periodic", verbose = FALSE,
                    min.scale = 0, precond = TRUE) {
    if (verbose == TRUE) 
        cat("wd: Argument checking...")
    if (!is.atomic(data)) 
        stop("Data is not atomic")
    DataLength <- length(data)
    nlevels <- nlevelsWT(data)
    filter <- wavethresh::filter.select(filter.number, family)
    filter$H <- wavethresh::filter.select(filter.number,family)$H^2
    filter$G <- wavethresh::filter.select(filter.number,family)$H^2
    if (is.na(nlevels)) 
        stop("Data length is not power of two")
    if (type != "wavelet" && type != "station") 
        stop("Unknown type of wavelet decomposition")
    if (type == "station" && bc != "periodic") 
        stop("Can only do periodic boundary conditions with station")
    if (verbose == TRUE) 
        cat("...done\nFilter...")
    if (verbose == TRUE) 
        cat("...selected\nFirst/last database...")
    fl.dbase <- wavethresh::first.last(LengthH = length(filter$H),
                                       DataLength = DataLength, 
                                       type = type, bc = bc)
    if (bc == "interval") {
        ans <- wavethresh::wd.int(data = data,
                                  preferred.filter.number = filter.number, 
                                  min.scale = min.scale,
                                  precond = precond)
        fl.dbase <- wavethresh::first.last(LengthH = length(filter$H),
                                           DataLength = DataLength, 
            type = type, bc = bc, current.scale = min.scale)
        filter <- list(name = paste("CDV", filter.number, sep = ""), 
            family = "CDV", filter.number = filter.number)
        l <-
          list(transformed.vector = ans$transformed.vector, 
               current.scale = ans$current.scale,
               filters.used = ans$filters.used, 
               preconditioned = ans$preconditioned, date = ans$date, 
               nlevels =
                 wavethresh::IsPowerOfTwo(length(ans$transformed.vector)), 
               fl.dbase = fl.dbase, type = type, bc = bc, filter = filter)
        class(l) <- "wd"
        return(l)
    }
    dtsp <- tsp(data)
    C <- rep(0, fl.dbase$ntotal)
    C[1:DataLength] <- data
    if (verbose == TRUE) 
        error <- 1
    else error <- 0
    if (verbose == TRUE) 
        cat("built\n")
    if (verbose == TRUE) 
        cat("Decomposing...\n")
    nbc <- switch(bc, periodic = 1, symmetric = 2)
    if (is.null(nbc)) 
        stop("Unknown boundary condition")
    ntype <- switch(type, wavelet = 1, station = 2)
    if (is.null(filter$G)) {
        wavelet.decomposition <-
          .C("wavedecomp", C = as.double(C),
             D = as.double(rep(0, fl.dbase$ntotal.d)),
             H = as.double(filter$H), 
             LengthH = as.integer(length(filter$H)),
             nlevels = as.integer(nlevels), 
             firstC = as.integer(fl.dbase$first.last.c[, 1]), 
             lastC = as.integer(fl.dbase$first.last.c[, 2]),
             offsetC = as.integer(fl.dbase$first.last.c[,3]),
             firstD = as.integer(fl.dbase$first.last.d[, 
                1]), lastD = as.integer(fl.dbase$first.last.d[, 
                2]), offsetD = as.integer(fl.dbase$first.last.d[, 
                3]), ntype = as.integer(ntype), nbc = as.integer(nbc), 
            error = as.integer(error))
    }
    else {
        wavelet.decomposition <-
          .C("comwd", CR = as.double(Re(C)), 
             CI = as.double(Im(C)), LengthC = as.integer(fl.dbase$ntotal), 
             DR = as.double(rep(0, fl.dbase$ntotal.d)),
             DI = as.double(rep(0, fl.dbase$ntotal.d)),
             LengthD = as.integer(fl.dbase$ntotal.d), 
             HR = as.double(Re(filter$H)), HI = as.double(-Im(filter$H)), 
             GR = as.double(Re(filter$G)), GI = as.double(-Im(filter$G)), 
             LengthH = as.integer(length(filter$H)),
             nlevels = as.integer(nlevels), 
             firstC = as.integer(fl.dbase$first.last.c[, 1]), 
             lastC = as.integer(fl.dbase$first.last.c[, 2]),
             offsetC = as.integer(fl.dbase$first.last.c[,3]),
             firstD = as.integer(fl.dbase$first.last.d[,1]),
             lastD = as.integer(fl.dbase$first.last.d[,2]),
             offsetD = as.integer(fl.dbase$first.last.d[,3]),
             ntype = as.integer(ntype), nbc = as.integer(nbc), 
             error = as.integer(error))
    }
    if (verbose == TRUE) 
        cat("done\n")
    error <- wavelet.decomposition$error
    if (error != 0) {
        cat("Error ", error, " occured in wavedecomp\n")
        stop("Error")
    }
    if (is.null(filter$G)) {
        l <- list(C = wavelet.decomposition$C, D = wavelet.decomposition$D, 
            nlevels = nlevelsWT(wavelet.decomposition), fl.dbase = fl.dbase, 
            filter = filter, type = type, bc = bc, date = date())
    }
    else {
        l <- list(C = complex(real = wavelet.decomposition$CR, 
                              imaginary = wavelet.decomposition$CI),
                  D = complex(real = wavelet.decomposition$DR, 
                              imaginary = wavelet.decomposition$DI),
                  nlevels = nlevelsWT(wavelet.decomposition), 
                  fl.dbase = fl.dbase, filter = filter, type = type, 
                  bc = bc, date = date())
    }
    class(l) <- "wd"
    if (!is.null(dtsp)) 
        tsp(l) <- dtsp
    return(l)
}

# This function converts variances of detail coefficients for a given
# wavelet basis from wd objects to wst objects.
#
#' @importFrom wavethresh wst
#' @importFrom wavethresh getarrvec
#' @importFrom wavethresh accessD
#' @importFrom wavethresh accessC
#' @importFrom wavethresh putD
#' @importFrom wavethresh putC
convert.var <- function (wd, ...) {
    if (wd$type != "station") 
        stop("Object to convert must be of type \"station\" ")
    n <- 2^nlevelsWT(wd)
    dummy <- rep(0, n)
    tmpwst <- wavethresh::wst(dummy, filter.number = wd$filter$filter.number, 
                              family = wd$filter$family)
    tmpwst$filter <- wd$filter
    tmpwst$date <- wd$date
    arrvec <- wavethresh::getarrvec(nlevelsWT(wd), sort = FALSE)
    for (lev in (nlevelsWT(wd) - 1):1) {
        ds <- wavethresh::accessD.wd(wd, level = lev)
        cs <- wavethresh::accessC.wd(wd, level = lev)
        ds <- ds[arrvec[, nlevelsWT(wd) - lev]]
        cs <- cs[arrvec[, nlevelsWT(wd) - lev]]
        tmpwst <- wavethresh::putD(tmpwst, level = lev, v = ds)
        tmpwst <- wavethresh::putC(tmpwst, level = lev, v = cs)
    }
    tmpwst <- wavethresh::putC(tmpwst, level = nlevelsWT(wd),
                               v = accessC(wd, level = wd$nlevels))
    tmpwst <- wavethresh::putD(tmpwst, level = nlevelsWT(wd),
                               v = accessC(wd, level = wd$nlevels))
    tmpwst <- wavethresh::putC(tmpwst, level = 0,
                               v = wavethresh::accessC(wd, level = 0))
    arrvec <- sort.list(wavethresh::levarr(1:n, levstodo = nlevelsWT(wd)))
    tmpwst <- wavethresh::putD(tmpwst, level = 0,
                               v = wavethresh::accessD(wd, level = 0)[arrvec])
    return(tmpwst)
}

# This function reconstructs variance of the mean estimate for a
# given wavelet basis given a wst object.
#
#' @importFrom wavethresh av.basis
#' @importFrom wavethresh nlevelsWT
AvBasis.var <- function (wst, Ccode = TRUE, ...) {
    nlevels <- nlevelsWT(wst)
    if (is.null(wst$filter$G)) {
        if (Ccode == FALSE) {
            answer <- av.basis(wst, level = nlevels - 1, ix1 = 0, 
                ix2 = 1, filter = wst$filter)
        }
        else {
            error <- 0
            answer <- rep(0, 2^nlevels)
            H <- wst$filter$H
            aobj <-
              .C("av_basisWRAP", wstR = as.double(wst$wp), 
                 wstC = as.double(wst$Carray),
                 LengthData = as.integer(length(answer)), 
                 level = as.integer(nlevels - 1), H = as.double(H), 
                 LengthH = as.integer(length(H)), answer = as.double(answer), 
                 error = as.integer(error))
            if (aobj$error != 0) 
                stop(paste("av_basisWRAP returned error code", 
                  aobj$error))
            answer <- aobj$answer
        }
    }
    else {
        error <- 0
        answerR <- answerI <- rep(0, 2^nlevels)
        H <- wst$filter$H
        G <- wst$filter$G
        aobj <- .C("comAB_WRAP", wstR = as.double(Re(wst$wp)), 
                   wstI = as.double(Im(wst$wp)),
                   wstCR = as.double(Re(wst$Carray)), 
                   wstCI = as.double(Im(wst$Carray)),
                   LengthData = as.integer(length(answerR)), 
                   level = as.integer(nlevels - 1), HR = as.double(Re(H)), 
                   HI = as.double(Im(H)), GR = as.double(Re(G)),
                   GI = as.double(Im(G)), LengthH = as.integer(length(H)),
                   answerR = as.double(answerR), answerI = as.double(answerI),
                   error = as.integer(error))
        if (aobj$error != 0) 
            stop(paste("av_basisWRAP returned error code", aobj$error))
        answer <- aobj$answerR
    }
    answer
}
