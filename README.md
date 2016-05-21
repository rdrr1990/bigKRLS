# bigKRLS

bigKRLS is an algorithm for Kernel Regularized Least Squares that incorporates big data packages 
for size and C++ for speed. 

## Supported Operating Systems
To date, bigKRLS has only been successfully run on Mac OS X Yosemite 10.10.5 and Linux Ubuntu 14.04. Due to discrepancies between the .dll and .so files that C++ generates, we are not releasing bigKRLS for Windows just yet but hope to soon. Thoughts on how to do this smoothly are most welcome!


## Installation

bigKRLS is designed to run on R version 3.3.0 ("Supposedly Educational" released 2016-05-03) and to be built with R Studio 0.99.48. Older, even fairly recent, versions of R will not work with bigmemory; the newest version of RStudio seems to do better helping Rcpp and RcppArmadillo find the C and Fortran on which the R to C++ interface relies. 

We each had installation issues and recommend doing the following in this sequence. (Thoughts on how to streamline are of course very welcome!)

1. Install R 3.3.0, available at https://cran.r-project.org 

2. Install RStudio 0.99.48, available at https://www.rstudio.com/products/rstudio/download/

3. Run the following commands:

install.packages("devtools")
library(devtools)
install.packages(c("Rcpp", "RcppArmadillo"))

The next step should be taken care of by the build, but can't hurt to run:

install.packages(c("bigmemory", "biganalytics", "bigalgebra"))

4. Many people have trouble running Rcpp and RcppArmadillo because R seems to have trouble finding C and Fortran (i.e., get errors like "cannot find lquadmath" or "lqfortran" or similar troubles with clang++ and g++). The documentation for those packages now deals with that topic extensively and once installation is complete our makevars and namespace files should take care of connecting everything under the hood. 

  -- At this point, we suspect Ubuntu users will be good to go but Apple users will still need to run the following two terminal commands:

curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2

sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /

  -- If troubles persist, we found the following pages particularly helpful:

http://thecoatlessprofessor.com/programming/setting-up-rstudio-to-work-with-rcpparmadillo/

http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/

http://web.mit.edu/insong/www/pdf/rpackage_instructions.pdf

5. Decompress the bigKRLS .tar file. Don't use install.packages; instead open bigKRLS.Rproject from RStudio and build it (command shift B).

6. You should be good to go; as a reminder, despite improvements, the algorithm is still incredibly memory intensive. As the documentation details, we recommend proceeding cautiously (say N = 5k for a laptop with 8 gigs of RAM) and bearing in mind that memory usage is a quadratic function of the number of observations. 

