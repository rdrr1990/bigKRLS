bigKRLS and KRLS convergence
================

This code establishes the numeric convergence of bigKRLS and KRLS estimates.

``` r
library(bigKRLS)

set.seed(2018)
N <- 500  
P <- 6
X <- matrix(rnorm(N*P), ncol=P)
X[,P] <- ifelse(X[,P] > 0.12345, 1, 0)
b <- runif(ncol(X))
y <- X %*% b + rnorm(nrow(X))

KRLS.out <- KRLS::krls(X = X, y = y, print.level = 0, eigtrunc=0.01)
bigKRLS.out <- bigKRLS(y, X, instructions = FALSE, eigtrunc=0.01)
```

    ....................

``` r
cor(KRLS.out$coeffs, bigKRLS.out$coeffs)    # should be 1
```

         [,1]
    [1,]    1

``` r
KRLS.out$avgderivatives
```

                x1        x2          x3         x4        x5        x6
    [1,] 0.2286663 0.1150259 0.006574909 0.09488611 0.3828897 0.7653918

``` r
bigKRLS.out$avgderivatives
```

            x1        x2          x3         x4        x5        x6
     0.2286663 0.1150259 0.006574909 0.09488611 0.3828897 0.7653918

``` r
max(abs(bigKRLS.out$derivatives - KRLS.out$derivatives)) < 0.00000001
```

    [1] TRUE

``` r
sessioninfo::session_info()
```

    ─ Session info ──────────────────────────────────────────────────────────
     setting  value                       
     version  R version 3.4.2 (2017-09-28)
     os       macOS Sierra 10.12.6        
     system   x86_64, darwin15.6.0        
     ui       X11                         
     language (EN)                        
     collate  en_US.UTF-8                 
     tz       America/Los_Angeles         
     date     2018-03-09                  

    ─ Packages ──────────────────────────────────────────────────────────────
     package       * version    date       source                            
     backports       1.1.1      2017-09-25 CRAN (R 3.4.2)                    
     bigalgebra      0.8.4      2014-04-16 CRAN (R 3.4.0)                    
     biganalytics    1.1.14     2016-02-18 CRAN (R 3.4.0)                    
     bigKRLS       * 2.1.0      2018-03-09 CRAN (R 3.4.2)                    
     biglm           0.9-1      2013-05-16 CRAN (R 3.4.0)                    
     bigmemory     * 4.5.19     2016-03-28 CRAN (R 3.4.0)                    
     bigmemory.sri * 0.1.3      2014-08-18 CRAN (R 3.4.0)                    
     clisymbols      1.2.0      2017-05-21 cran (@1.2.0)                     
     codetools       0.2-15     2016-10-05 CRAN (R 3.4.2)                    
     colorspace      1.3-2      2016-12-14 CRAN (R 3.4.0)                    
     DBI             0.7        2017-06-18 CRAN (R 3.4.0)                    
     digest          0.6.13     2017-12-14 cran (@0.6.13)                    
     evaluate        0.10.1     2017-06-24 CRAN (R 3.4.1)                    
     foreach         1.4.3      2015-10-13 CRAN (R 3.4.0)                    
     ggplot2         2.2.1.9000 2017-12-15 Github (tidyverse/ggplot2@bfff1d8)
     gtable          0.2.0      2016-02-26 CRAN (R 3.4.0)                    
     htmltools       0.3.6      2017-04-28 CRAN (R 3.4.0)                    
     httpuv          1.3.5      2017-07-04 CRAN (R 3.4.1)                    
     iterators       1.0.8      2015-10-13 CRAN (R 3.4.0)                    
     knitr           1.17       2017-08-10 CRAN (R 3.4.1)                    
     KRLS            1.0-0      2017-07-10 CRAN (R 3.4.1)                    
     lazyeval        0.2.1      2017-10-29 CRAN (R 3.4.2)                    
     magrittr        1.5        2014-11-22 CRAN (R 3.4.0)                    
     mime            0.5        2016-07-07 CRAN (R 3.4.0)                    
     munsell         0.4.3      2016-02-13 CRAN (R 3.4.0)                    
     pillar          1.2.1      2018-02-27 cran (@1.2.1)                     
     plyr            1.8.4      2016-06-08 CRAN (R 3.4.0)                    
     R6              2.2.2      2017-06-17 CRAN (R 3.4.0)                    
     Rcpp            0.12.14    2017-11-23 cran (@0.12.14)                   
     rlang           0.2.0      2018-02-20 cran (@0.2.0)                     
     rmarkdown       1.6        2017-06-15 CRAN (R 3.4.0)                    
     rprojroot       1.2        2017-01-16 CRAN (R 3.4.0)                    
     scales          0.5.0.9000 2017-12-06 Github (hadley/scales@d767915)    
     sessioninfo     1.0.0      2017-06-21 CRAN (R 3.4.1)                    
     shiny           1.0.5      2017-08-23 cran (@1.0.5)                     
     stringi         1.1.5      2017-04-07 CRAN (R 3.4.0)                    
     stringr         1.2.0      2017-02-18 CRAN (R 3.4.0)                    
     tibble          1.4.2      2018-01-22 cran (@1.4.2)                     
     withr           2.1.0.9000 2017-12-06 Github (jimhester/withr@fe81c00)  
     xtable          1.8-2      2016-02-05 CRAN (R 3.4.0)                    
     yaml            2.1.16     2017-12-12 cran (@2.1.16)
