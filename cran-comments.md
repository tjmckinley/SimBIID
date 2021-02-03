## Minor submission

This is a minor submission in which I have changed a few tests to use
functions in the third edition of the `testthat` package.

## Test environments

* local Ubuntu 20.10 install, R 4.0.2
* R-oldrel (via travis-ci): Ubuntu 16.04.6 LTS, R 3.6.3
* R-release (via travis-ci): Ubuntu 16.04.6 LTS, R 4.0.2
* R-devel (via travis-ci): Ubuntu 16.04.6 LTS
* macOS (via travis-ci), High Sierra 10.13.6, R 4.0.3
* win-builder (devel, release [4.0.3] and old [3.6.3])
* Ubuntu Linux 20.04 LTS, R-release, GCC (via r-hub)
* Fedora Linux, R-devel, GCC, gfortran (via r-hub)
* Debian Linux, R-release, GCC ASAN/UBSAN (via r-hub)
* Oracle Solaris 10, x86, 32 bit, R-release (via r-hub)

## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:

```
> checking installed package size ... NOTE
    installed size is  5.2Mb
    sub-directories of 1Mb or more:
      libs   4.8Mb
```

I don't think this is an issue, since it's just the compiled code.

## Downstream dependencies

There are currently no downstream dependencies for this package.
