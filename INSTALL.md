# "big" Install
## bigKRLS Installation Guide by Pete Mohanty and Robert Shaffer

`bigKRLS` is an `R` algorithm for Kernel-Regularized Least Squares that uses big data packages for size and C++ for speed. This architecture takes a bit of work to set up because R isn't necessarily properly set to connect C++ to C and FORTRAN (i.e., R throws lengthy error messages about lquadmath, lqfortran, clang++, and/or g++) and older compilers aren't necessarily compatabible with parallel processing paradigms that are increasingly popular in the `R` community, like `OMP`. Once everything is connected under the hood you'll be able to install `R` packages like `rstan` and `biglasso` too. 

## Supported Platforms
`bigKRLS` has been run on multiple platforms including `Mac OS X` `Yosemite 10.10.5` and `Sierra 10.12.6`, `Linux Ubuntu 14.04`, and `Windows 7` and `Windows 8`.


## Pre-Requisites

`bigKRLS` is designed to run on `R version 3.3.0` ("Supposedly Educational" released 2016-05-03) or newer. Older, even fairly recent, versions of `R` will not work with `bigmemory`. 

-- Install the newest `R` at https://cran.r-project.org 

### Windows users must install up-to-date Rtools (3.3 or newer):

[https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/)

RTools [best practices](http://thecoatlessprofessor.com/programming/rcpp/install-rtools-for-rcpp/)

### Mac OSX 
Mac users will need to be sure their compilers are up-to-date. Without fairly current compilers (i.e., compilers that do not come standard with `OS X`), it will not be possible to install [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html). 

For the `g++` family, version `4.6.*` or newer is required. For `g++` and related software, see [The Coatless Professor's OpenMP in R and OSX](http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/). Everything up to the `clang4` instructions on that page are recommended.

Mac users will need `clang4` (the recent version is required for `OMP` multicore processing used by other `R` libraries such as `biglasso`). Though `clang4` can be installed with bash commands, the [installer](https://uofi.box.com/v/r-macos-clang-pkg) developed by the Coatless Professor (@coatless) is highly recommended since it automatically takes care of the configuration files and paths `R` requires. For detail, see [https://github.com/coatless/r-macos-clang](https://github.com/coatless/r-macos-clang). 

If troubles persist, we found the following pages particularly helpful: [A](http://thecoatlessprofessor.com/programming/setting-up-rstudio-to-work-with-rcpparmadillo/), [B](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/), and section 2.16 of: [C](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-FAQ.pdf).


### Current RStudio

To use RStudio, Windows users must use `RStudio` 1.1.129 or newer and Unix-type users (including Mac) must use 1.0.136 or newer. As of 2017-10-09, the current stable build works for both:

https://www.rstudio.com/products/rstudio/download/    


## The Environment
`bigKRLS` has several dependencies, some of which require recent version of their dependencies. To smooth installation, we recommend installing these packages first.

```
install.packages(c("Rcpp", "RcppArmadillo", "bigmemory", "biganalytics", "snow", "shiny", "httpuv", "scales", "lazyeval", "tibble")) 
```
## Install via CRAN
You should now be able to install via `CRAN`.
```
install.packages("bigKRLS")
library(bigKRLS)
vignette("bigKRLS_basics")
```

## Install via GitHub
You should now be also able to install `bigKRLS` via `GitHub`. 
```
install.packages("devtools")   
library(devtools)  
```
Windows users should first run these extra lines:
```{r, eval = F}
find_rtools()
find_rtools(T)  
```
Finally, install the most current version with standard devtools syntax:

```{r, eval = F}
install_github('rdrr1990/bigKRLS')
library(bigKRLS)
vignette("bigKRLS_basics")
```
You should be good to go!


## Memory Limits
Despite improvements, the algorithm is still incredibly memory intensive. We recommend proceeding cautiously and bearing in mind that memory usage is a quadratic function of the number of observations, N (roughly 5N<sup>2</sup>). Users should estimate models N = 1,000, 2,500, and 5,000 to see how their system performs before considering larger models. See https://sites.google.com/site/petemohanty/software for detail.

## License 
Code released under GPL (>= 2).


