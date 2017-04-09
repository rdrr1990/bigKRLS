\name{bigKRLS-package}
\alias{bigKRLS-package}
\alias{bigKRLS}
\docType{package}
\title{
\packageTitle{bigKRLS}
}
\description{
\packageDescription{bigKRLS}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{bigKRLS}
\packageIndices{bigKRLS}

Kernel Regularized Least Squares (KRLS) is a kernel-based, complexity-penalized method developed by Hainmueller and Hazlett (2013) to minimize parametric assumptions while maintaining interpretive clarity. Here, we introduce bigKRLS, an updated version of the original \code{\link[https://cran.r-project.org/web/packages/KRLS/index.html]{KRLS}} R package with algorithmic and implementation improvements designed to optimize speed and memory usage. These improvements allow users to straightforwardly fit KRLS models to medium and large datasets (N > ~2500). 

Major Updates

1. C++ integration. We re-implement most major computations in the model in C++ via Rcpp and RcppArmadillo. These changes produce up to a 50\% runtime decrease compared to the original R implementation.

2. Leaner algorithm. Because of the Tikhonov regularization and parameter tuning strategies used in KRLS, the method of estimation is inherently memory-heavy O(N^2), making memory savings important even in small- and medium-sized applications. We develop and implement a new local derivatives algorithm, which reduces peak memory usage by approximately an order of magnitude, and cut the number of computations needed to find regularization parameter in half.

3. Improved memory management. Most data objects in R perform poorly in memory-intensive applications. We use a series of packages in the bigmemory environment to ease this constraint, allowing our implementation to handle larger datasets more smoothly.

4. Parallel Processing. Parallel processing with snow makes the algorithm much faster for the marginal effects.

5. Interactive data visualization. We've designed an R Shiny app that allows users bigKRLS users to easily share results with collaborators or more general audiences. Simply call shiny.bigKRLS() on the outputted regression object. 

Requirements. bigKRLS is under active development, and currently requires R version 3.3.0 or later. Windows users should use RTools 3.3 or later. RStudio users must have a current version as well.

For details on syntax, load the library and then open our vignette vignette("bigKRLS_basics"). Because of the quadratic memory requirement, users working on a typical laptop (8-16 gigabytes of RAM) may wish to start at N = 2,500 or 5,000, particularly if the number of *x* variables is large. When you have a sense of how bigKRLS runs on your system, you may wish to only estimate a subset of the marginal effects at N = 10-15,000 by setting bigKRLS(... which.derivatives = c(1, 3, 5)) for the marginal effects of the first, third, and fifth *x* variable. 

Recent papers, slides, and other code available at github.com/rdrr1990/code/

You may also be interested in our recent presentation to the International Methods Colloquium, viewable at https://youtu.be/4WYDIXLUYbc

# License 
Code released under GPL (>= 2).

}
\author{
\packageAuthor{bigKRLS}

Maintainer: 
\packageMaintainer{bigKRLS}
}
\references{

Mohanty, Pete and Robert B. Shaffer. 2016. "Messy Data, Robust Inference? Navigating Obstacles to Inference with bigKRLS" Project Presented to the International Methods Colloquium and useR! 2016. Visit https://sites.google.com/site/petemohanty for most recent version.

Hainmueller, Jens and Chad Hazlett. 2014. "Kernel Regularized Least Squares: Reducing Misspecification Bias with a Flexible and Interpretable Machine Learning Approach." Political Analysis. 22:143-68. https://web.stanford.edu/~jhain/Paper/PA2014a.pdf (Accessed May 20th, 2016).
}
\keyword{ package }
\seealso{
https://sites.google.com/site/petemohanty/software
https://github.com/rdrr1990/bigKRLS   
... CRAN Packages:
\code{\link[https://cran.r-project.org/web/packages/KRLS/index.html]{KRLS}}
\code{\link[https://cran.r-project.org/web/packages/np/index.html]{NP}}
}
\examples{
%% ~~ simple examples of the most important functions ~~
}
