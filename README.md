# SimBIID

R package that implements various simulation-based inference routines for infectious disease models.

Package provides some code to run simulations of state-space models, and then use these in the ABC-SMC algorithm of Toni et al. (2009) and the bootstrap particle filter based particle MCMC algorithm (Andrieu et al., 2010). Also provides functions to plot and summarise the outputs.

## Installation

The package depends on the `Rcpp` and `RcppArmadillo` packages, which require the installation of the correct C++ compilers. The guidance below is taken from Sections 2.1.1, 2.1.2 or 2.1.3 here:

https://teuder.github.io/rcpp4everyone_en/020_install.html

### Windows

Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/index.html).

### Mac

Install Xcode command line tools. Execute the command `xcode-select --install` in a Terminal.

You might also need to install the gfortran libraries from:

https://cran.r-project.org/bin/macosx/tools/gfortran-6.1.pkg

### Linux

Install gcc and related packages.

In Ubuntu Linux, execute the command `sudo apt-get install r-base-dev` in a Terminal.

Once the compilers have been installed, then instal the `devtools` package in R and run:

```
library(devtools)
install_github("tjmckinley/SimBIID")
```

