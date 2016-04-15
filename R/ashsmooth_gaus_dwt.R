library(smash)


wt.mat = function (n, filter.number, family) 
{
  J = log2(n)
  X = diag(rep(1, n))
  W = matrix(0, n, n)
  W = apply(X, 1, smash:::wd.D, filter.number = filter.number, family = family)
  return(list(W2 = W^2))
}


ashsmooth.gaus.dwt = function(x, sigma, v.est = FALSE, joint = FALSE, v.basis = FALSE, post.var = FALSE, filter.number = 1, 
  family = "DaubExPhase", return.loglr = TRUE, jash = FALSE, SGD = TRUE, weight = 0.5, min.var = 1e-08, ashparam = list()) {
  
  n = length(x)
  J = log2(n)
  if (!isTRUE(all.equal(J, trunc(J)))) {
    stop("Error: number of columns of x must be power of 2")
  }
  
  ashparam = smash:::setAshParam.gaus(ashparam)
  if (v.est == TRUE) {
    weight = 1
  } else {
    weight = 0.5
  }
  
  # if(post.var==TRUE&basis[[1]]!='haar'){stop('Error: posterior variances returned only with Haar basis')}
  if (post.var == TRUE & v.est == TRUE & jash == TRUE) {
    stop("Error: Posterior variances for variance estimate not returned for method JASH")
  }
  if (joint == TRUE) {
    v.est = TRUE
  }
  
  
  tsum = sum(x)
  
  Wl = NULL
  x.w.d = wd(x, filter.number = filter.number, family = family)
  Wl = wt.mat(n, filter.number = filter.number, family = family)
  
  
  wc = x.w.d
  
  data.var = sigma^2
  
  wmean = matrix(0, J, n)
  wvar = matrix(0, J, n)
  
  logLR = 0
  
  x.w = wc
  x.w.v = apply((rep(1, n - 1) %o% data.var) * Wl$W2, 1, sum)  #diagonal of W*V*W' 
  x.w.v.s = rep(0, n * J)
  for (j in 0:(J - 1)) {
    x.pm = rep(0, 2^j)
    if(j!=(J-1)){
      index = (sum(n/(2^(1:(J-j-1))))+1):sum(n/(2^(1:(J-j))))
    }else{
      index = 1:sum(n/2)
    }
    x.w.j = accessD(x.w, j)
    x.w.v.j = x.w.v[index]
    ind.nnull = (x.w.v.j != 0)
    zdat.ash = smash:::shrink.wc(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]), prior = ashparam$prior, pointmass = ashparam$pointmass, 
                                 nullcheck = ashparam$nullcheck, VB = ashparam$VB, mixsd = ashparam$mixsd, mixcompdist = ashparam$mixcompdist, gridmult = ashparam$gridmult, jash = FALSE, 
                                 df = NULL, SGD = FALSE)
    x.pm[ind.nnull] = zdat.ash$PosteriorMean
    x.pm[!ind.nnull] = 0
    x.w = putD(x.w, j, x.pm)
    zdat.ash$model = "EE"
    logLR.temp = ashr:::calc_loglik(zdat.ash,x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]),NULL) -  
      sum(dnorm(x.w.j[ind.nnull], 0, sqrt(x.w.v.j[ind.nnull]), log = TRUE))
    logLR[j+1] = logLR.temp
  }
  
  mu.est = wr(x.w)
  logLR = sum(logLR)

  return(list(mu.est = mu.est, logLR = logLR))
  
}
