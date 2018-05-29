# bigKRLS
[![Travis-CI Build Status](https://travis-ci.org/rdrr1990/bigKRLS.svg?branch=master)](https://travis-ci.org/rdrr1990/bigKRLS)
[![Coverage Status](https://img.shields.io/codecov/c/github/rdrr1990/bigKRLS/master.svg)](https://codecov.io/github/rdrr1990/bigKRLS?branch=master) [![cran checks](https://cranchecks.info/badges/summary/bigKRLS)](https://cranchecks.info/pkgs/bigKRLS) [![Rdoc](http://www.rdocumentation.org/badges/version/bigKRLS)](http://www.rdocumentation.org/packages/bigKRLS) ![](http://cranlogs.r-pkg.org/badges/grand-total/bigKRLS)
 
Kernel Regularized Least Squares (KRLS) is a kernel-based, complexity-penalized method developed by [Hainmueller and Hazlett (2013)](http://pan.oxfordjournals.org/content/22/2/143), and designed to minimize parametric assumptions while maintaining interpretive clarity. Here, we introduce `bigKRLS`, an updated version of the original [KRLS R package](https://CRAN.R-project.org/package=KRLS) with algorithmic and implementation improvements designed to optimize speed and memory usage. These improvements allow users to straightforwardly estimate pairwise regression models with KRLS once *N > 2500*. Since April 15, 2017, `bigKRLS` has been available on  [CRAN](http://www.rdocumentation.org/badges/version/bigKRLS). You may also be interested in our [working paper](https://web.stanford.edu/~pmohanty/mohanty_shaffer_workingpaper.pdf), which has been accepted by *Political Analysis*, and which demonstrates the utility of `bigKRLS` by analyzing the 2016 US presidential election. Our replication materials can be found on [Dataverse](https://doi.org/10.7910/DVN/CYYLOK) and [Github repo](https://github.com/rdrr1990/bigKRLS/tree/master/examples) contains examples too.

# Major Updates found in bigKRLS

1. C++ integration. We re-implement most major computations in the model in `C++` via [Rcpp](https://CRAN.R-project.org/package=Rcpp) and [RcppArmadillo](https://CRAN.R-project.org/package=RcppArmadillo). These changes produce up to a 50% runtime decrease compared to the original `R` implementation.

2. Leaner algorithm. Because of the Tikhonov regularization and parameter tuning strategies used in KRLS, the method of estimation is inherently memory-heavy (O($N^2$)), making memory savings important even in small- and medium-sized applications. We develop and implement a [new marginal effects algorithm](https://github.com/rdrr1990/bigKRLS/blob/master/src/bigderiv_v3.cpp), which reduces peak memory usage by approximately an order of magnitude, and cut the number of computations needed to find regularization parameter [in half](https://github.com/rdrr1990/bigKRLS/blob/master/src/solveforc.cpp).

3. Improved memory management. Most data objects in `R` perform poorly in memory-intensive applications. We use a series of packages in the [bigmemory](https://CRAN.R-project.org/package=bigmemory) environment to ease this constraint, allowing our implementation to handle larger datasets more smoothly.

4. Parallel Processing. In addition to the single-core algorithmic improvements, parallel processing obtains the pointwise marginal effects substantially faster.

5. Interactive data visualization. We've designed an `R` [Shiny](shiny.rstudio.com) app that allows users `bigKRLS` users to easily share results with collaborators or more general audiences. Simply call `shiny.bigKRLS()`. 

6. Honest p values. `bigKRLS` now computes p values that reflect both the regularization process and the number of predictors. For details on how the effective sample size is calculated as well as other options, see `help(summary.bigKRLS)`.

```
out <- bigKRLS(y, X)
out$Neffective
summary(out)
```

7. Cross-validation, including K folds crossvalidation. `crossvalidate.bigKRLS` performs CV, stores a number of in and out of sample statistics, as well as metadata documenting how data the were split and the bigmemory file structure (if applicable). 

```
cv <- crossvalidate.bigKRLS(y, X, seed = 2017, ptesting = 20)
kcv <- crossvalidate.bigKRLS(y, X, seed = 2017, Kfolds = 5)
``` 
See `vignette("bigKRLS_basics")` for details.

8. Eigentruncation. `bigKRLS` now supports two types of eigentruncation to decrease runtime.

```
out <- bigKRLS(y, X, eigtrunc = 0.001)     # defaults to 0.001 if N > 3000 and 0 otherwise
out <- bigKRLS(y, X, Neig = 100)           # only compute 100 vecs and vals (defaults to Neig = nrow(X))
```
 




# Installation

`bigKRLS` requires a series of packages--notably `bigmemory`, `Rcpp`, and `RcppArmadillo`--current versions of which require up-to-date versions of `R` *and* its compilers (`RStudio`, if used, must be current as well). To install the latest stable version from `CRAN`:
```
install.packages("bigKRLS")
```
To install the `GitHub` version, use standard devtools syntax:

```
install.packages("devtools")
library(devtools)
install_github('rdrr1990/bigKRLS')
```
New users may wish to [see our installation notes for specifics](https://github.com/rdrr1990/code/blob/master/bigKRLS_installation.md)



# Getting Going...
For details on syntax, load the library and then open our vignette:
```
library(bigKRLS)
vignette("bigKRLS_basics")
```
Because of the quadratic memory requirement, users working on a typical laptop (8-16 gigabytes of RAM) may wish to start at N = 2,500 or 5,000, particularly if the number of *x* variables is large. When you have a sense of how bigKRLS runs on your system, you may wish to only estimate a subset of the marginal effects at N = 10-15,000 by setting `bigKRLS(..., which.derivatives = c(1, 3, 5))` for the marginal effects of the first, third, and fifth *x* variable. 

# License 
Code released under GPL (>= 2).
