# Disclaimer: The functions defined in this file are adapted from the same functions defined in package 'Multiseq'
# under the GPL license, authored by Ester Pantaleo, Heejung Shim, Matthew Stephens and Zhengrong Xing

ff_exp = function(x) {
    return(exp(-log(1 + exp(-x))))
}

ffdash = function(x) {
    return(1/(1 + exp(x)))
}

ffdd = function(x) {
    ex = exp(x)
    return(-ex/(1 + ex)^2)
}

ffdashdd = function(x) {
    ex = exp(x)
    return(-ex/(1 + ex)^3)
}

# first and second derivatives of ff
ff.deriv_exp = function(x) {
    ex = exp(x)
    return(list(first = ex/(1 + 2 * ex + ex^2), second = (ex - ex^3)/(1 + 4 * ex + 6 * ex^2 + 4 * ex^3 + ex^4)))
}

# first and second derivatives of ffdash
ffdash.deriv = function(x) {
    ex = exp(x)
    return(list(first = -ex/(1 + ex)^2, second = ex * (ex - 1) * (1 + ex)^(-3)))
}

ffdd.deriv = function(x) {
    ex = exp(x)
    return(list(first = ex * (ex - 1) * (1 + ex)^(-3), second = ex * (4 * ex - ex^2 - 1) * (1 + ex)^(-4)))
}

ffdashdd.deriv = function(x) {
    ex = exp(x)
    return(list(first = ex * (2 * ex - 1) * (1 + ex)^(-4), second = ex * (7 * ex - 4 * ex^2 - 1) * (1 + ex)^(-5)))
}

# return estimated E(ff(X)) and var(ff(X)) by delta method INPUT: mx =E(X) and vx = Var(X)
ff.moments_exp = function(mx, vx) {
    d = ff.deriv_exp(mx)
    return(list(mean = ff_exp(mx) + d$second * vx/2, var = (d$first^2 * vx), meansq = ff_exp(mx)^2 + (d$second * vx/2)^2 + 
        (d$first^2 + ff_exp(mx) * d$second) * vx))
}

# return estimated E(ff'(X)) and var(ff'(X)) by delta method INPUT: mx =E(X) and vx = Var(X)
ffdash.moments = function(mx, vx) {
    d = ffdash.deriv(mx)
    return(list(mean = ffdash(mx) + d$second * vx/2, var = (d$first^2 * vx), meansq = ffdash(mx)^2 + (d$second * vx/2)^2 + 
        (d$first^2 + ffdash(mx) * d$second) * vx))
}

# return estimated E(ff''(X)) and var(ff''(X)) by delta method INPUT: mx =E(X) and vx = Var(X)
ffdd.moments = function(mx, vx) {
    d = ffdd.deriv(mx)
    return(list(mean = ffdd(mx) + d$second * vx/2, var = (d$first^2 * vx), meansq = ffdd(mx)^2 + (d$second * vx/2)^2 + 
        (d$first^2 + ffdd(mx) * d$second) * vx))
}

# return estimated E(ff'(X)*ff''(X)) and var(ff'(X)*ff''(X)) by delta method INPUT: mx =E(X) and vx = Var(X)
ffdashdd.moments = function(mx, vx) {
    d = ffdashdd.deriv(mx)
    return(list(mean = ffdashdd(mx) + d$second * vx/2, var = (d$first^2 * vx), meansq = ffdashdd(mx)^2 + (d$second * 
        vx/2)^2 + (d$first^2 + ffdashdd(mx) * d$second) * vx))
}


logit = function(x) log(x/(1 - x))


fl = function(x, n) {
    return(logit(x/n))
}


fl.deriv = function(x, n) {
    return(list(first = n/(x * (n - x))))
}

fl.moments = function(mx, vx, n) {
    d = fl.deriv(mx, n)
    return(d$first^2 * vx)
}



ff = function(x) {
    return(-log(1 + exp(-x)))
}

ffdash = function(x) {
    return(1/(1 + exp(x)))
}

# first and second derivatives of ff
ff.deriv = function(x) {
    ex = exp(x)
    return(list(first = 1/(1 + ex), second = -ex/(1 + ex)^2))
}

# first and second derivatives of ffdash
ffdash.deriv = function(x) {
    ex = exp(x)
    return(list(first = -ex/(1 + ex)^2, second = ex * (ex - 1) * (1 + ex)^(-3)))
}

# return estimated E(ff(X)) and var(ff(X)) by delta method INPUT: mx =E(X) and vx = Var(X)
ff.moments = function(mx, vx) {
    d = ff.deriv(mx)
    return(list(mean = ff(mx) + d$second * vx/2, var = (d$first^2 * vx), meansq = ff(mx)^2 + (d$second * vx/2)^2 + 
        (d$first^2 + ff(mx) * d$second) * vx))
}

# return estimated E(ff'(X)) and var(ff'(X)) by delta method INPUT: mx =E(X) and vx = Var(X)
ffdash.moments = function(mx, vx) {
    d = ffdash.deriv(mx)
    return(list(mean = ffdash(mx) + d$second * vx/2, var = (d$first^2 * vx), meansq = ffdash(mx)^2 + (d$second * vx/2)^2 + 
        (d$first^2 + ffdash(mx) * d$second) * vx))
} 
