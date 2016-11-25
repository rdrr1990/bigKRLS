# bigKRLS
Kernel Regularized Least Squares (KRLS) is a kernel-based, complexity-penalized method developed by [Hainmueller and Hazlett (2013)](http://pan.oxfordjournals.org/content/22/2/143), and designed to minimize parametric assumptions while maintaining interpretive clarity. Here, we introduce *bigKRLS*, an updated version of the original [KRLS R package](https://cran.r-project.org/web/packages/KRLS/index.html) with algorithmic and implementation improvements designed to optimize speed and memory usage. These improvements allow users to straightforwardly fit KRLS models to medium and large datasets (N > ~2500). 

# Major Updates

1. C++ integration. We re-implement most major computations in the model in C++ via [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html). These changes produce up to a 50% runtime decrease compared to the original R implementation.

2. Leaner algorithm. Because of the Tikhonov regularization and parameter tuning strategies used in KRLS, the model is inherently memory-heavy (O(N<sup>2</sup>)), making marginal memory savings important even in small- and medium-sized applications. We develop and implement a new local derivatives algorithm, which reduces peak memory usage by approximately an order of magnitude.

3. Improved memory management. Most data objects in R perform poorly in memory-intensive applications. We use a series of packages in the [bigmemory](https://cran.r-project.org/web/packages/bigmemory/index.html) environment to ease this constraint, allowing our implementation to handle larger datasets more smoothly.

4. Interactive data visualization. We've designed an R [Shiny](shiny.rstudio.com) app that allows users bigKRLS users to easily share results with collaborators or more general audiences. Simply call shiny.bigKRLS() on the outputted regression object. 


# Installation
bigKRLS is under active development, and currently requires R version 3.3.0 or later. Windows users should use RTools 3.3 or later; Windows users should use R (not RStudio). To install, use standard devtools syntax:

```
install.packages("devtools")
library(devtools)
install_github('rdrr1990/bigKRLS')
```

# Dependencies
bigKRLS requires Rcpp and RcppArmadillo, as well as a series of packages in the bigmemory environment. Users new to these packages may wish see our [installation notes](https://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnxwZXRlbW9oYW50eXxneDozYTA1ZGRjZmJkZWY0YWI4).

# Getting Going...
For details on syntax, load the library and then open our vignette:
```
library(bigKRLS)
vignette("bigKRLS_basics")
```
Recent slides and other code available at github.com/rdrr1990/code/

You may also be interested in our recent presentation to the International Methods Colloquium, viewable at https://youtu.be/4WYDIXLUYbc

# License 
Code released under GPL (>= 2).
