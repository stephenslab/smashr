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

2. Load the smashr package, and run the smashr demo.

   ```R
   library(smashr)
   demo("smashr")
   ```
   
3. To learn more, see the smashr package help:

   ```R
   help(package = smashr)
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
