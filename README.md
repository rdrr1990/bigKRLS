# bigKRLS

bigKRLS is an R algorithm for Kernel Regularized Least Squares that uses big data packages 
for size and C++ for speed. 

## Supported Operating Systems
To date, bigKRLS has been successfully run on Mac OS X Yosemite 10.10.5, Linux Ubuntu 14.04, and Windows 7. Due to discrepancies between the .dll and .so files that C++ generates, we have compiled a separate source package for Windows. 


## Installation

bigKRLS is designed to run on R version 3.3.0 ("Supposedly Educational" released 2016-05-03) and to be built with R Studio 0.99.48. Older, even fairly recent, versions of R will not work with bigmemory; the newest version of RStudio does better connecting R and C++. 

We each had installation issues and recommend doing the following in this sequence. Thoughts on how to streamline are of course very welcome!

1. Install R 3.3.0, available at https://cran.r-project.org 

2. Install RStudio 0.99.48, available at https://www.rstudio.com/products/rstudio/download/

3. Run the following commands:

install.packages("devtools")  
library(devtools)  
install.packages(c("Rcpp", "RcppArmadillo"))  

The next step should be taken care of by the build, but can't hurt to run:

install.packages(c("bigmemory", "biganalytics", "bigalgebra"))

### Missed Connections
Many people have trouble running Rcpp and RcppArmadillo because R isn't set to find C and Fortran (i.e., R throws errors like "cannot find lquadmath" or "lqfortran" or similar troubles with clang++ and g++). Once installation is complete our makevars and namespace files should take care of connecting everything under the hood. 

### Ubuntu 
Ubuntu users should be good to go at this point

### Mac OSX 
Mac OSX users may still need to run the following two terminal commands:

curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2

sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /

If troubles persist, we found the following pages particularly helpful:

http://thecoatlessprofessor.com/programming/setting-up-rstudio-to-work-with-rcpparmadillo/

http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/


### Windows users must install Rtools:

https://cran.r-project.org/bin/windows/Rtools/  

Finally, before installing from the .tar, run the following:

library(devtools)  
find_rtools()  
find_rtools(T)  

http://stackoverflow.com/questions/19885381/rtools-not-being-detected-by-r

## Build or Install
You should now be able to install bigKRLS for your operating system. For example:

install.packages("~/Downloads/bigKRLS_1.1.tar.gz", repos = NULL, type = "source")

Alternatively, you may wish to clone or download the files into a folder "bigKRLS" and then use RStudio to start a new project from those existing files and build the file from there.  
  
You should be good to go!

## Memory Limits
Despite improvements, the algorithm is still incredibly memory intensive. As the documentation details, we recommend proceeding cautiously (say N = 5k for a laptop with 8 gigs of available RAM) and bearing in mind that memory usage is a quadratic function of the number of observations. 



