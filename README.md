# smashr: smoothing using Adaptive Shrinkage in R

This R package implements fast, wavelet-based Empirical Bayes
shrinkage methods for signal denoising. This includes smoothing
Poisson-distributed data and Gaussian-distributed data, with possibly
heteroskedastic error. The algorithms implement the methods described
in [Xing & Stephens (2016)][smash-arxiv].

If you find a bug, please post an [issue][issues].

## License

Copyright (c) 2016-2018, Zhengrong Xing, Peter Carbonetto and Matthew
Stephens.

All source code and software in this repository is free software; you
can redistribute it and/or modify it under the terms of the
[GNU General Public License][gpl] as published by the
[Free Software Foundation][fsf]; either version 3 of the License, or
(at your option) any later version. See the [LICENSE](LICENSE) file
for the full text of the license.

## Citing this work

If you find that this R package useful for your work, please cite our
paper:

> Zhengrong Xing and Matthew Stephens (2016). *Smoothing via Adaptive
> Shrinkage (smash): denoising Poisson and heteroskedastic Gaussian
> signals.* [arXiv:1605.07787](https://arxiv.org/abs/1605.07787).

## Quick Start

Follow these steps to quickly get started using smashr.

1. In R, install the latest version of smashr using
   [devtools](https://github.com/r-lib/devtools):

   ```R
   install.packages("devtools")
   library(devtools)
   install_github("stephenslab/smashr")
   ```

2. Run the smashr demo.


3. To learn more, see the smashr package help:

   ```R
   help(package = smashr)
   ```
   
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
lines(mu.est, col = 2)

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

## Credits

This R package was developed by
[Zhengrong Xing](https://github.com/zrxing) and
[Matthew Stephens](http://stephenslab.uchicago.edu) at the University
of Chicago, with contributions from
[Peter Carbonetto](http://pcarbo.github.io).

[smash-arxiv]: http://arxiv.org/abs/1605.07787
[issues]: https://github.com/stephenslab/smashr/issues
[gpl]: http://www.gnu.org/licenses/gpl.html
[fsf]: https://www.fsf.org
