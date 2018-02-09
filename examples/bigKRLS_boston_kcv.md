K folds crossvalidating with bigKRLS
================

This little demo shows how to use and interpret `crossvalidate.bigKRLS` on the Boston housing data. (`medv` is median house value).

``` r
library(pacman)
p_load(bigKRLS, tidyverse, knitr, MASS)
opts_chunk$set(tidy = TRUE, comment = "")

glimpse(Boston)
```

    ## Observations: 506
    ## Variables: 14
    ## $ crim    <dbl> 0.00632, 0.02731, 0.02729, 0.03237, 0.06905, 0.02985, ...
    ## $ zn      <dbl> 18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.5, 12.5, 12.5, 12.5,...
    ## $ indus   <dbl> 2.31, 7.07, 7.07, 2.18, 2.18, 2.18, 7.87, 7.87, 7.87, ...
    ## $ chas    <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
    ## $ nox     <dbl> 0.538, 0.469, 0.469, 0.458, 0.458, 0.458, 0.524, 0.524...
    ## $ rm      <dbl> 6.575, 6.421, 7.185, 6.998, 7.147, 6.430, 6.012, 6.172...
    ## $ age     <dbl> 65.2, 78.9, 61.1, 45.8, 54.2, 58.7, 66.6, 96.1, 100.0,...
    ## $ dis     <dbl> 4.0900, 4.9671, 4.9671, 6.0622, 6.0622, 6.0622, 5.5605...
    ## $ rad     <int> 1, 2, 2, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, ...
    ## $ tax     <dbl> 296, 242, 242, 222, 222, 222, 311, 311, 311, 311, 311,...
    ## $ ptratio <dbl> 15.3, 17.8, 17.8, 18.7, 18.7, 18.7, 15.2, 15.2, 15.2, ...
    ## $ black   <dbl> 396.90, 396.90, 392.83, 394.63, 396.90, 394.12, 395.60...
    ## $ lstat   <dbl> 4.98, 9.14, 4.03, 2.94, 5.33, 5.21, 12.43, 19.15, 29.9...
    ## $ medv    <dbl> 24.0, 21.6, 34.7, 33.4, 36.2, 28.7, 22.9, 27.1, 16.5, ...

``` r
y <- as.matrix(Boston$medv)
X <- as.matrix(Boston %>% dplyr::select(-medv))

out <- crossvalidate.bigKRLS(y, X, seed = 1234, Kfolds = 5)
```

    ## .........................
    ## .........................
    ## .........................
    ## ........................
    ## .........................

``` r
s <- summary(out)
s[["overview"]]
# kable(s[['overview']] %>% format(digits = 2, scientific = FALSE))
# s[['overview']] %>% format(digits = 2, scientific = FALSE) %>% kable
```

|                         | Fold 1  | Fold 2   | Fold 3   | Fold 4  | Fold 5  |
|:------------------------|:--------|:---------|:---------|:--------|:--------|
| MSE (In Sample)         | 5.547   | 5.339    | 5.332    | 4.303   | 5.504   |
| MSE (Out of Sample)     | 8.464   | 10.896   | 8.970    | 17.715  | 7.753   |
| MSE AME (In Sample)     | 815.101 | 1388.626 | 1080.492 | 489.638 | 801.770 |
| MSE AME (Out of Sample) | 754.644 | 1724.905 | 1328.164 | 392.617 | 753.858 |
| R2 (In Sample)          | 0.937   | 0.938    | 0.940    | 0.944   | 0.935   |
| R2 (Out of Sample)      | 0.883   | 0.868    | 0.871    | 0.857   | 0.910   |
| R2 AME (In Sample)      | 0.068   | 0.048    | 0.070    | 0.046   | 0.068   |
| R2 AME (Out of Sample)  | 0.085   | 0.120    | 0.019    | 0.356   | 0.108   |

The predictiveness of the model is relatively stable across models (85-90% of variance explained). The model fits the best for Fold 5 (lowest Mean Squared Error and highest *R*<sup>2</sup> out of sample). The model is massively non-additive in that giant gaps between how much the model as a whole explains vs. how much the average marginal effects (AMEs) explain. Put differently, most of the variance is explained by interactions. The AMEs are most informative for Fold 4.

``` r
summary(out)
```

    Overview of Model Performance

    N: 506 
    Kfolds: 5 
    Seed: 1234 

                              Fold 1   Fold 2   Fold 3   Fold 4   Fold 5
    MSE (In Sample)           5.5473 5.34e+00 5.33e+00   4.3028   5.5040
    MSE (Out of Sample)       8.4643 1.09e+01 8.97e+00  17.7145   7.7533
    MSE AME (In Sample)     815.1014 1.39e+03 1.08e+03 489.6375 801.7701
    MSE AME (Out of Sample) 754.6435 1.72e+03 1.33e+03 392.6165 753.8584
    R2 (In Sample)            0.9365 9.38e-01 9.40e-01   0.9438   0.9345
    R2 (Out of Sample)        0.8830 8.68e-01 8.71e-01   0.8567   0.9102
    R2 AME (In Sample)        0.0679 4.84e-02 7.03e-02   0.0463   0.0677
    R2 AME (Out of Sample)    0.0853 1.20e-01 1.86e-02   0.3563   0.1085

    MSE denotes Mean Squared Error. AME implies calculations done with Average Marginal Effects only.

    Summary of Training Model1:


    MODEL SUMMARY:

    lambda: 0.2145 
    N: 404 
    N Effective: 310.8598 
    R2: 0.9365 
    R2AME**: 0.0679 

    Average Marginal Effects:

            Estimate Std. Error t value Pr(>|t|)
    crim     -0.0977     0.0324 -3.0207   0.0027
    zn       -0.0118     0.0108 -1.0968   0.2736
    indus     0.0038     0.0332  0.1153   0.9083
    chas*   103.6046    18.9535  5.4662   0.0000
    nox     -10.4226     2.7829 -3.7453   0.0002
    rm        4.4998     0.3261 13.8007   0.0000
    age      -0.0418     0.0082 -5.1280   0.0000
    dis      -1.2196     0.1250 -9.7595   0.0000
    rad       0.1590     0.0404  3.9340   0.0001
    tax      -0.0104     0.0016 -6.2975   0.0000
    ptratio  -0.2368     0.0692 -3.4223   0.0007
    black     0.0013     0.0037  0.3400   0.7341
    lstat    -0.3262     0.0392 -8.3303   0.0000


    Percentiles of Marginal Effects:

                  5%      25%     50%      75%      95%
    crim     -0.3903  -0.1510 -0.0604  -0.0156   0.0713
    zn       -0.0534  -0.0296 -0.0087  -0.0012   0.0459
    indus    -0.1600  -0.0659 -0.0155   0.0550   0.2094
    chas*   -68.0724  29.2456 95.8535 179.0750 297.9587
    nox     -31.2882 -14.0401 -8.9859  -3.9585   6.0274
    rm       -1.8657   0.5196  4.1433   8.1936  11.5775
    age      -0.1160  -0.0725 -0.0414  -0.0115   0.0476
    dis      -2.7053  -1.6986 -1.1764  -0.6760  -0.0046
    rad      -0.0130   0.0940  0.1613   0.2340   0.3409
    tax      -0.0247  -0.0160 -0.0103  -0.0032   0.0007
    ptratio  -1.2797  -0.5753 -0.1221   0.1730   0.3994
    black    -0.0071  -0.0022  0.0011   0.0048   0.0103
    lstat    -1.1181  -0.4552 -0.2188  -0.0766   0.0395

    (*) Reported average and percentiles of dy/dx is for discrete change of the dummy variable from min to max (usually 0 to 1)).


    (**) Pseudo-R^2 computed using only the Average Marginal Effects.


    You may also wish to use predict() for out-of-sample forecasts or shiny.bigKRLS() to interact with results. Type vignette("bigKRLS_basics") for sample syntax. Use save.bigKRLS() to store results and load.bigKRLS() to re-open them.



    Summary of Training Model2:


    MODEL SUMMARY:

    lambda: 0.2145 
    N: 405 
    N Effective: 313.9762 
    R2: 0.9377 
    R2AME**: 0.0484 

    Average Marginal Effects:

            Estimate Std. Error t value Pr(>|t|)
    crim     -0.0517     0.0240 -2.1497   0.0324
    zn       -0.0092     0.0101 -0.9120   0.3625
    indus     0.0020     0.0336  0.0608   0.9515
    chas*   147.9221    24.3181  6.0828   0.0000
    nox     -13.0208     2.7971 -4.6551   0.0000
    rm        4.4851     0.3136 14.3003   0.0000
    age      -0.0310     0.0080 -3.8886   0.0001
    dis      -1.1857     0.1231 -9.6309   0.0000
    rad       0.1869     0.0391  4.7763   0.0000
    tax      -0.0086     0.0016 -5.2814   0.0000
    ptratio  -0.3032     0.0735 -4.1265   0.0000
    black     0.0018     0.0038  0.4702   0.6385
    lstat    -0.3683     0.0392 -9.4026   0.0000


    Percentiles of Marginal Effects:

                   5%      25%      50%      75%      95%
    crim      -0.2288  -0.0955  -0.0174   0.0029   0.0503
    zn        -0.0514  -0.0330  -0.0111  -0.0015   0.0545
    indus     -0.1600  -0.0727  -0.0138   0.0559   0.2338
    chas*   -141.6465  33.1156 111.6110 294.1739 472.7116
    nox      -32.5640 -19.2169 -12.9771  -5.1220   5.2963
    rm        -1.1989   1.0614   3.7261   7.9043  11.2197
    age       -0.0964  -0.0568  -0.0372  -0.0094   0.0687
    dis       -2.5854  -1.6466  -1.0910  -0.6789  -0.1640
    rad        0.0142   0.1224   0.1789   0.2572   0.3693
    tax       -0.0223  -0.0138  -0.0084  -0.0022   0.0022
    ptratio   -1.2227  -0.6010  -0.1910   0.0713   0.2279
    black     -0.0146  -0.0030   0.0041   0.0075   0.0126
    lstat     -1.1521  -0.5104  -0.2757  -0.1081   0.0384

    (*) Reported average and percentiles of dy/dx is for discrete change of the dummy variable from min to max (usually 0 to 1)).


    (**) Pseudo-R^2 computed using only the Average Marginal Effects.


    You may also wish to use predict() for out-of-sample forecasts or shiny.bigKRLS() to interact with results. Type vignette("bigKRLS_basics") for sample syntax. Use save.bigKRLS() to store results and load.bigKRLS() to re-open them.



    Summary of Training Model3:


    MODEL SUMMARY:

    lambda: 0.2145 
    N: 405 
    N Effective: 312.0236 
    R2: 0.9397 
    R2AME**: 0.0703 

    Average Marginal Effects:

            Estimate Std. Error t value Pr(>|t|)
    crim     -0.0679     0.0252 -2.6980   0.0074
    zn       -0.0118     0.0105 -1.1221   0.2627
    indus    -0.0175     0.0330 -0.5306   0.5961
    chas*   130.8481    23.8693  5.4819   0.0000
    nox      -7.3249     2.8442 -2.5754   0.0105
    rm        4.5547     0.3115 14.6239   0.0000
    age      -0.0415     0.0082 -5.0442   0.0000
    dis      -1.0450     0.1243 -8.4035   0.0000
    rad       0.1334     0.0396  3.3670   0.0009
    tax      -0.0103     0.0016 -6.2477   0.0000
    ptratio  -0.3054     0.0717 -4.2578   0.0000
    black     0.0032     0.0037  0.8674   0.3864
    lstat    -0.3375     0.0398 -8.4821   0.0000


    Percentiles of Marginal Effects:

                   5%      25%      50%      75%      95%
    crim      -0.2665  -0.1124  -0.0443  -0.0080   0.0522
    zn        -0.0458  -0.0312  -0.0086  -0.0010   0.0390
    indus     -0.1760  -0.0973  -0.0320   0.0393   0.1923
    chas*   -112.2719  30.8496 118.4596 242.5509 401.1752
    nox      -29.3183 -11.4084  -6.1125   0.0483   9.2596
    rm        -1.1173   1.3039   4.1210   7.9353  11.1596
    age       -0.1098  -0.0684  -0.0463  -0.0203   0.0512
    dis       -2.1330  -1.5585  -1.0529  -0.5551   0.0137
    rad        0.0002   0.0884   0.1395   0.1802   0.2705
    tax       -0.0267  -0.0155  -0.0100  -0.0029   0.0006
    ptratio   -1.3377  -0.5942  -0.1856   0.0899   0.3627
    black     -0.0135  -0.0008   0.0044   0.0083   0.0135
    lstat     -1.0287  -0.4540  -0.2768  -0.1079   0.0490

    (*) Reported average and percentiles of dy/dx is for discrete change of the dummy variable from min to max (usually 0 to 1)).


    (**) Pseudo-R^2 computed using only the Average Marginal Effects.


    You may also wish to use predict() for out-of-sample forecasts or shiny.bigKRLS() to interact with results. Type vignette("bigKRLS_basics") for sample syntax. Use save.bigKRLS() to store results and load.bigKRLS() to re-open them.



    Summary of Training Model4:


    MODEL SUMMARY:

    lambda: 0.2235 
    N: 405 
    N Effective: 311.3317 
    R2: 0.9438 
    R2AME**: 0.0463 

    Average Marginal Effects:

            Estimate Std. Error t value Pr(>|t|)
    crim     -0.0575     0.0236 -2.4341   0.0155
    zn       -0.0146     0.0092 -1.5960   0.1116
    indus    -0.0076     0.0297 -0.2550   0.7989
    chas*    76.2714    19.6619  3.8791   0.0001
    nox      -9.3801     2.4636 -3.8075   0.0002
    rm        4.6342     0.2778 16.6812   0.0000
    age      -0.0391     0.0071 -5.5207   0.0000
    dis      -0.8727     0.1072 -8.1380   0.0000
    rad       0.1596     0.0361  4.4222   0.0000
    tax      -0.0092     0.0015 -6.1669   0.0000
    ptratio  -0.2719     0.0621 -4.3755   0.0000
    black     0.0021     0.0035  0.5899   0.5557
    lstat    -0.3188     0.0350 -9.1168   0.0000


    Percentiles of Marginal Effects:

                   5%      25%     50%      75%      95%
    crim      -0.2153  -0.0979 -0.0333  -0.0015   0.0232
    zn        -0.0552  -0.0336 -0.0154  -0.0039   0.0519
    indus     -0.1482  -0.0723 -0.0242   0.0345   0.1882
    chas*   -164.7485  14.6677 96.0850 153.0838 245.1769
    nox      -25.7832 -15.6798 -9.3561  -2.4803   7.1671
    rm        -1.4179   1.4054  4.1174   7.8686  11.2416
    age       -0.0982  -0.0656 -0.0477  -0.0183   0.0471
    dis       -1.5947  -1.2085 -0.8921  -0.5497  -0.0349
    rad       -0.0099   0.0884  0.1644   0.2346   0.3192
    tax       -0.0214  -0.0141 -0.0094  -0.0035   0.0006
    ptratio   -1.2984  -0.5896 -0.1512   0.0863   0.3100
    black     -0.0103  -0.0016  0.0023   0.0063   0.0129
    lstat     -0.8951  -0.4616 -0.2629  -0.1201   0.0199

    (*) Reported average and percentiles of dy/dx is for discrete change of the dummy variable from min to max (usually 0 to 1)).


    (**) Pseudo-R^2 computed using only the Average Marginal Effects.


    You may also wish to use predict() for out-of-sample forecasts or shiny.bigKRLS() to interact with results. Type vignette("bigKRLS_basics") for sample syntax. Use save.bigKRLS() to store results and load.bigKRLS() to re-open them.



    Summary of Training Model5:


    MODEL SUMMARY:

    lambda: 0.2145 
    N: 405 
    N Effective: 311.4581 
    R2: 0.9345 
    R2AME**: 0.0677 

    Average Marginal Effects:

            Estimate Std. Error t value Pr(>|t|)
    crim     -0.0748     0.0295 -2.5335   0.0118
    zn       -0.0080     0.0102 -0.7845   0.4334
    indus    -0.0035     0.0321 -0.1089   0.9134
    chas*   103.3547    20.1556  5.1279   0.0000
    nox      -9.6089     2.7275 -3.5230   0.0005
    rm        4.4201     0.3310 13.3535   0.0000
    age      -0.0279     0.0086 -3.2618   0.0012
    dis      -1.0859     0.1202 -9.0362   0.0000
    rad       0.1522     0.0413  3.6885   0.0003
    tax      -0.0091     0.0016 -5.5111   0.0000
    ptratio  -0.2777     0.0685 -4.0541   0.0001
    black     0.0012     0.0039  0.3008   0.7638
    lstat    -0.3888     0.0401 -9.7057   0.0000


    Percentiles of Marginal Effects:

                  5%      25%     50%      75%      95%
    crim     -0.3058  -0.1314 -0.0450  -0.0063   0.0627
    zn       -0.0415  -0.0256 -0.0060  -0.0001   0.0434
    indus    -0.1683  -0.0802 -0.0121   0.0516   0.2005
    chas*   -93.7708  29.8238 94.8950 177.4629 309.1480
    nox     -32.0905 -14.5180 -8.7140  -1.4732   7.9573
    rm       -1.4243   1.3069  3.3908   7.7720  11.6479
    age      -0.0923  -0.0594 -0.0337  -0.0008   0.0619
    dis      -2.3266  -1.4901 -1.0199  -0.6114  -0.0825
    rad      -0.0077   0.0856  0.1518   0.2171   0.3327
    tax      -0.0225  -0.0142 -0.0089  -0.0027   0.0012
    ptratio  -1.2556  -0.5653 -0.1724   0.0907   0.3276
    black    -0.0150  -0.0032  0.0024   0.0060   0.0130
    lstat    -1.2151  -0.4925 -0.3027  -0.1570  -0.0414

    (*) Reported average and percentiles of dy/dx is for discrete change of the dummy variable from min to max (usually 0 to 1)).


    (**) Pseudo-R^2 computed using only the Average Marginal Effects.


    You may also wish to use predict() for out-of-sample forecasts or shiny.bigKRLS() to interact with results. Type vignette("bigKRLS_basics") for sample syntax. Use save.bigKRLS() to store results and load.bigKRLS() to re-open them.
