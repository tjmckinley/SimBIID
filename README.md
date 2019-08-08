# SimBIID

R package that implements various simulation-based inference routines for infectious disease models.

Package provides some code to run simulations of state-space models, and then use these in the ABC-SMC algorithm of Toni et al. (2009) and the bootstrap particle filter based particle MCMC algorithm (Andrieu et al., 2010). Also provides functions to plot and summarise the outputs.

## Installation

The package uses `Rcpp`, which requires the installation of the correct C++ compilers. Please see Sections 2.1.1, 2.1.2 or 2.1.3 here:

https://teuder.github.io/rcpp4everyone_en/020_install.html

Once the compilers have been installed, then using the `devtools` package, you can run:

```
library(devtools)
install_github("tjmckinley/SimBIID")
```

