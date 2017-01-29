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

bigKRLS estimates a model via Kernel Regularized Least Squares with support for big matrices. 

Installation is more complex than usual, particularly for those who have not worked with Rcpp before; see the read me for details. 

The syntax is simple: bigKRLS(y = y, X=X). However, users should proceed with caution. If you are reading this, this computationally-intensive algorithm has not been extensively tested!! 

The memory usage of this algorithm is proportional to N, the number of observations, squared. Once N crosses 12,000 or so, each of the multiple required N * N matrices will be a gigabyte. Though bigmemory should automatically shift to the hard drive as necessary, the comprimises in terms of speed and stability are unknown. Feedback as to what works and what does not on particular machines and operating systems is very welcome (pete.mohanty@gmail.com). 

It is strongly recommended that users start at a modest N (say 5k for a machine with 8 gigs of RAM) and work up incrementally, bearing the quadratic memory requirement in mind.

Missing data are not allowed; users with missing data may wish to use an R package such as mice (multiple imputation by chained equations). 
}
\author{
\packageAuthor{bigKRLS}

Maintainer: \packageMaintainer{bigKRLS}
}
\references{

Mohanty, Pete and Robert B. Shaffer. 2016. "Optimizing non-parametric regression in R." Presentation Presented to Bay Area Users Group Official May Meetup. Slides available at https://sites.google.com/site/petemohanty/software/bigKRLS_for_BARUG.html

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
