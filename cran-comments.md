## Minor submission

This is a minor submission in which I have added additional 
functionality, and amended the way tolerances are updated.
Specifically I have:

* added bisection method for choosing tolerances when multiple
data points are being used in the ABC-SMC algorithm;

* added option to set minimum tolerances for ABC-SMC.

I've also made some patch changes in the way that random seeds are 
passed to `mclapply()`. This is to aid reproducibility when setting 
seeds and using parallelisation. Runs will only be reproducible if 
using the same number of cores each time (which can be specified 
using `mc.cores` argument to various functions).

## Test environments

* local Ubuntu 20.04 install, R 3.6.3
* R-oldrel (via travis-ci): Ubuntu 16.04.6 LTS, R 3.6.3
* R-release (via travis-ci): Ubuntu 16.04.6 LTS, R 4.4.0
* R-devel (via travis-ci): Ubuntu 16.04.6 LTS
* macOS (via travis-ci), High Sierra 10.13.6, R 4.4.0
* win-builder (devel and release [4.4.0])
* Ubuntu Linux 16.04 LTS, R-release, GCC (via r-hub)
* Fedora Linux, R-devel, clang, gfortran (via r-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (via r-hub)
* Oracle Solaris 10, x86, 32 bit, R-patched (experimental) (via r-hub)

## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:
> checking installed package size ... NOTE
    installed size is  5.3Mb
    sub-directories of 1Mb or more:
      libs   4.9Mb

I don't think this is an issue, since it's just the compiled code.

## Downstream dependencies

There are currently no downstream dependencies for this package.
