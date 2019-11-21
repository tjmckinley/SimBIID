## Patch submission

This is a patch submission in which I have:

* fixed the order in which `tspan` objects are updated in code 
  produced from `mparseRcpp`.

## Test environments

* local Ubuntu 19.04 install, R 3.5.2
* R-oldrel (via travis-ci): Ubuntu 16.04.6 LTS, R 3.5.3
* R-release (via travis-ci): Ubuntu 16.04.6 LTS, R 3.6.1
* R-devel (via travis-ci): Ubuntu 16.04.6 LTS
* macOS (via travis-ci), High Sierra 10.13.3, R 3.6.1
* win-builder (devel and release)

* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (via r-hub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (via r-hub)
* Fedora Linux, R-devel, clang, gfortran (via r-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (via r-hub)
* Oracle Solaris 10, x86, 32 bit, R-patched (experimental) (via r-hub)

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
