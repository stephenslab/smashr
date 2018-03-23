# Disclaimer: All TI table (eg 'ParentTItable')
# and reconstruction (eg 'reverse.pwave') functions and their respective Rcpp counterparts are ported
# into R from Matlab functions 'BMSMShrink' and 'TISumProd' as part of the BMSM project maintained by Eric Kolaczyk

# @description Turns a TItable into a matrix that one can apply glm
#   to (for Poisson denoising).
# 
# @param dmat: A k by n matrix.
# 
# @return A 2 by n*k/2 matrix.
# 
simplify = function (dmat)
  matrix(t(dmat), nrow = 2)

#' @title Create a TI table and a parent table.
#' 
#' @description This function returns both a TItable and a "parent"
#'   table whose pairwise comparisons are used to create a TI table. For
#'   example, in the ith row, elements 1, 2 are the parents of the first
#'   element in the (i+1)th row of the TI table.  This function creates
#'   a decomposition table of signal, using pairwise sums, keeping just
#'   the values that are *not* redundant under the shift-invariant
#'   scheme.
#'
#' @param sig: an n vector of Poisson counts at n locations
#' 
#' @return a list with elements "TItable" and "parent"
#' 
#' @details This is very similar to TI-tables in Donoho and
#'   Coifman's TI-denoising framework.x
#' 
#' @keywords internal
#' 
#' @export
#' 
ParentTItable = function (sig) {
  n = length(sig)
  J = log2(n)
  
  dmat = matrix(0, nrow = J + 1, ncol = n)
  dmat[1, ] = sig

  # Initialize the parent table.
  dmat2 = matrix(0, nrow = J, ncol = 2 * n)  
  
  for (D in 0:(J - 1)) {
    nD = 2^(J - D)
    nDo2 = nD/2
    twonD = 2 * nD
    for (l in 0:(2^D - 1)) {
      ind = (l * nD + 1):((l + 1) * nD)
      ind2 = (l * twonD + 1):((l + 1) * twonD)
      x = dmat[D + 1, ind]
      lsumx = x[seq(from = 1, to = nD - 1, by = 2)] +
              x[seq(from = 2, to = nD, by = 2)]
      rx = rshift(x)
      rsumx = rx[seq(from = 1, to = nD - 1, by = 2)] +
              rx[seq(from = 2, to = nD, by = 2)]
      dmat[D + 2, ind] = c(lsumx, rsumx)
      dmat2[D + 1, ind2] = c(x, rx)
    }
  }
  return(list(TItable = dmat, parent = dmat2))
}

#' @title Reverse wavelet transform a set of probabilities in TItable
#'   format for Poisson data.
#' 
#' @param est An n-vector. Usually a constant vector with each
#'   element equal to the estimated log(total rate).
#' 
#' @param lp A J by n matrix of probabilities.
#' 
#' @param lq A J by n matrix of complementary probabilities.
#' 
#' @return Reconstructed signal in the original data space.
#' 
#' @export
#' 
reverse.pwave = function(est, lp, lq = NULL) {
  if (is.null(lq))
    lq = log(1 - exp(lp))
  if (length(est) == 1)
    est = rep(est, ncol(lp))
  
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
      nestl = interleave(estl + lpl, estl + lql)
      
      estr = estvec[(nDo2 + 1):nD]
      lpr = lpvec[(nDo2 + 1):nD]
      lqr = lqvec[(nDo2 + 1):nD]
      nestr = interleave(estr + lpr, estr + lqr)
      nestr = lshift(nestr)
      
      est[ind] = 0.5 * (nestl + nestr)
    }
  }
  return(est)
}
