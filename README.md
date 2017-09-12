# bigKRLS
[![Travis-CI Build Status](https://travis-ci.org/rdrr1990/bigKRLS.svg?branch=master)](https://travis-ci.org/rdrr1990/bigKRLS)
[![Coverage Status](https://img.shields.io/codecov/c/github/rdrr1990/bigKRLS/master.svg)](https://codecov.io/github/rdrr1990/bigKRLS?branch=master)

Kernel Regularized Least Squares (KRLS) is a kernel-based, complexity-penalized method developed by [Hainmueller and Hazlett (2013)](http://pan.oxfordjournals.org/content/22/2/143), and designed to minimize parametric assumptions while maintaining interpretive clarity. Here, we introduce *bigKRLS*, an updated version of the original [KRLS R package](https://CRAN.R-project.org/package=KRLS) with algorithmic and implementation improvements designed to optimize speed and memory usage. These improvements allow users to straightforwardly estimate pairwise regression models with KRLS once N > ~2500. *bigKRLS* is now available on CRAN:  
[![Rdoc](http://www.rdocumentation.org/badges/version/bigKRLS)](http://www.rdocumentation.org/packages/bigKRLS) ![](http://cranlogs.r-pkg.org/badges/grand-total/bigKRLS)

# Major Updates

1. C++ integration. We re-implement most major computations in the model in C++ via [Rcpp](https://CRAN.R-project.org/package=Rcpp) and [RcppArmadillo](https://CRAN.R-project.org/package=RcppArmadillo). These changes produce up to a 50% runtime decrease compared to the original R implementation.

2. Leaner algorithm. Because of the Tikhonov regularization and parameter tuning strategies used in KRLS, the method of estimation is inherently memory-heavy (O(N<sup>2</sup>)), making memory savings important even in small- and medium-sized applications. We develop and implement a [new marginal effects algorithm](https://github.com/rdrr1990/code/blob/master/mohanty_shaffer_IMC.pdf), which reduces peak memory usage by approximately an order of magnitude, and cut the number of computations needed to find regularization parameter [in half](https://github.com/rdrr1990/code/blob/master/solveforc.pdf).

3. Improved memory management. Most data objects in R perform poorly in memory-intensive applications. We use a series of packages in the [bigmemory](https://CRAN.R-project.org/package=bigmemory) environment to ease this constraint, allowing our implementation to handle larger datasets more smoothly.

4. Parallel Processing. Parallel processing with [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) makes the algorithm much faster for the marginal effects.

5. Interactive data visualization. We've designed an R [Shiny](shiny.rstudio.com) app that allows users bigKRLS users to easily share results with collaborators or more general audiences. Simply call shiny.bigKRLS() on the outputted regression object. 

For more detail, you may be interested in reading our [working paper](https://people.stanford.edu/pmohanty/sites/default/files/mohanty_shaffer_bigkrls_paper.pdf) or watching our [latest presentation](https://www.youtube.com/watch?v=4WYDIXLUYbc).

# New with bigKRLS 2.0.0 (now on CRAN)

1. Honest p values. `bigKRLS` now computes p values that reflect both the regularization process and the number of predictors. For details and other options, see `help(summary.bigKRLS)`.

```
out <- bigKRLS(y, X)
summary(out)
```

2. Cross-validation, including K folds crossvalidation. `crossvalidate.bigKRLS` performs CV, stores a number of in and out of sample statistics, as well as metadata documenting how the were split, the bigmemory file structure (if appropriate), and so on. See `vignette("bigKRLS_basics")` for syntax.


# Installation
`bigKRLS` is under active development. `bigKRLS` requires `R version 3.3.0` or later. Windows users should use `RTools 3.3` or later. To use `RStudio`, `Windows` must use [RStudio 1.1.129](https://dailies.rstudio.com/) or newer. To install the latest stable version from `CRAN`:
```
install.packages("bigKRLS")
```

To instead install the newest version from GitHub, use standard devtools syntax:

```
install.packages("devtools")
library(devtools)
install_github('rdrr1990/bigKRLS')
```

## Dependencies
`bigKRLS` requirea a series of packages in the `bigmemory` environment as well as `Rcpp` and `RcppArmadillo`, current versions of which require up-to-date compilers. New users may wish to see our [installation notes](https://github.com/rdrr1990/code/blob/master/bigKRLS_installation.md).

# Getting Going...
For details on syntax, load the library and then open our vignette:
```
library(bigKRLS)
vignette("bigKRLS_basics")
```
Because of the quadratic memory requirement, users working on a typical laptop (8-16 gigabytes of RAM) may wish to start at N = 2,500 or 5,000, particularly if the number of *x* variables is large. When you have a sense of how bigKRLS runs on your system, you may wish to only estimate a subset of the marginal effects at N = 10-15,000 by setting bigKRLS(... which.derivatives = c(1, 3, 5)) for the marginal effects of the first, third, and fifth *x* variable. 

Recent slides and other code available at https://github.com/rdrr1990/code/

You may also be interested in our recent presentation to the International Methods Colloquium, viewable at https://youtu.be/4WYDIXLUYbc

# License 
Code released under GPL (>= 2).
