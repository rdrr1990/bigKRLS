library(bigKRLS)

# this script illustrates the cumulative impact of different approaches to eigendecomposition, truncation
# including lambdasearch and variance-covariance calculations 
# (but not kernel calculation, derivative estimation, etc.)

edata <- read.csv("https://raw.githubusercontent.com/rdrr1990/bigKRLS/master/examples/data2016GE.csv")

y <- as.matrix(edata[,1])
X <- as.matrix(edata[,-1])
dim(X)
bX <- as.big.matrix(X)

KRLS.K.time <- system.time(KRLS::gausskernel(X,ncol(X)))[3]
bKRLS.K.time <- system.time(bigKRLS:::bGaussKernel(bX))[3]

bKRLS.time <- system.time(bigKRLS(y,X, eigtrunc=0, derivative = FALSE, noisy=FALSE, instructions = FALSE))[3]
KRLS.time <- system.time(KRLS::krls(X,y, derivative=FALSE))[3]

bKRLS.time.trunc <- system.time(bigKRLS(y,X, derivative = FALSE, noisy=FALSE, instructions = FALSE))[3]
KRLS.time.trunc <- system.time(KRLS::krls(X,y, derivative=FALSE,eigtrunc=0.01))[3]

bKRLS.time.Neig <- system.time(bigKRLS(y,X, Neig=50,eigtrunc=0.01,  
                                       derivative = FALSE, noisy=FALSE, instructions = FALSE))[3]

remove(bX)

wallclock <- matrix(nrow=2, ncol=3)
rownames(wallclock) <- c("KRLS", "bigKRLS")
colnames(wallclock) <- c("Full Decomp", "Eigentruncation", "Estimating Fewer")

wallclock[1,1] <- KRLS.time - KRLS.K.time
wallclock[2,1] <- bKRLS.time - bKRLS.K.time

wallclock[1,2] <- KRLS.time.trunc - KRLS.K.time
wallclock[2,2] <- bKRLS.time.trunc - bKRLS.K.time

wallclock[2,3] <- bKRLS.time.Neig - bKRLS.K.time

wallclock



# wallclock on rice.stanford.edu
#          Full Decomp Eigentruncation Estimating Fewer
# KRLS         56.882          43.859               NA
# bigKRLS     144.045          31.389           18.907
