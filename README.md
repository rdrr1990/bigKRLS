# bigKRLS
Kernel Regularized Least Squares (KRLS) is a kernel-based, complexity-penalized method developed by [Hainmueller and Hazlett (2013)](http://pan.oxfordjournals.org/content/22/2/143), and designed to minimize parametric assumptions while maintaining interpretive clarity. Here, we introduce *bigKRLS*, an updated version of the original [KRLS R package](https://CRAN.R-project.org/package=KRLS) with algorithmic and implementation improvements designed to optimize speed and memory usage. These improvements allow users to straightforwardly estimate pairwise regression models with KRLS once N > ~2500. 

# Major Updates

1. C++ integration. We re-implement most major computations in the model in C++ via [Rcpp](https://CRAN.R-project.org/package=Rcpp) and [RcppArmadillo](https://CRAN.R-project.org/package=RcppArmadillo). These changes produce up to a 50% runtime decrease compared to the original R implementation.

2. Leaner algorithm. Because of the Tikhonov regularization and parameter tuning strategies used in KRLS, the method of estimation is inherently memory-heavy (O(N<sup>2</sup>)), making memory savings important even in small- and medium-sized applications. We develop and implement a [new marginal effects algorithm](https://github.com/rdrr1990/code/blob/master/mohanty_shaffer_IMC.pdf), which reduces peak memory usage by approximately an order of magnitude, and cut the number of computations needed to find regularization parameter [in half](https://github.com/rdrr1990/code/blob/master/solveforc.pdf).

3. Improved memory management. Most data objects in R perform poorly in memory-intensive applications. We use a series of packages in the [bigmemory](https://CRAN.R-project.org/package=bigmemory) environment to ease this constraint, allowing our implementation to handle larger datasets more smoothly.

4. Parallel Processing. Parallel processing with [snow](https://CRAN.R-project.org/package=snow) makes the algorithm much faster for the marginal effects.

5. Interactive data visualization. We've designed an R [Shiny](shiny.rstudio.com) app that allows users bigKRLS users to easily share results with collaborators or more general audiences. Simply call shiny.bigKRLS() on the outputted regression object. 

For more detail, you may be interested in reading our [working paper](https://people.stanford.edu/pmohanty/sites/default/files/mohanty_shaffer_bigkrls_paper.pdf) or watching our [latest presentation](https://www.youtube.com/watch?v=4WYDIXLUYbc).

# Installation
bigKRLS is under active development, and currently requires R version 3.3.0 or later. Windows users should use RTools 3.3 or later. To use RStudio, Windows must use RStudio 1.1.129 or newer. To install, use standard devtools syntax:

```
install.packages("devtools")
library(devtools)
install_github('rdrr1990/bigKRLS')
```

## Dependencies
bigKRLS requires Rcpp and RcppArmadillo, as well as a series of packages in the bigmemory environment. Users new to these packages may wish see our [installation notes](https://github.com/rdrr1990/code/blob/master/bigKRLS_installation.md).

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
