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

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* Maintainer: 'Trevelyan J. McKinley <t.mckinley@exeter.ac.uk>'
   New submission
     Andrieu (10:43)
     al (9:51, 10:54)
     SMC (9:26)
     et (9:48, 10:51)
   Possibly mis-spelled words in DESCRIPTION:
   
This arises since it's a new submission to CRAN, and the highlighted words
are not mispelled.

There was one additional NOTE **only found** in the Windows Server 2008 R2 SP1, 
R-devel, 32/64 bit build submitted via **R-Hub**. This was:

* checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    'examples_x64' 'tests_i386' 'tests_x64'
    'SimBIID-Ex_i386.Rout' 'SimBIID-Ex_x64.Rout' 'examples_i386'
    
These files are not in the source code, and I assume must have been
built during the compilation on the server. This note does not arise
on any other build, and so I assume is an artefact of the R-development
package on the Windows server.

## Downstream dependencies
There are currently no downstream dependencies for this package.