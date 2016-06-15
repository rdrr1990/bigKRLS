#' Kernel Regularized Least Squares with Big Matrices
#' 
#' @param y A vector of observations on the dependent variable; missing values not allowed 
#' @param X A matrix of observations of the independent variables; missing values and constant vectors (e.g. intercept) not allowed
#' @param lambda Regularization parameter. Default: estimated based (in part) on the eigenvalues of the kernel
#' @param sigma Bandwidth parameter, shorthand for sigma squared. Default: sigma <- ncol(X). Since x variables are standardized, faciltates interprepation of the kernel.
#' @param eigtrunc Proportion of eigen object to be truncated. Must be between 0 (no truncation) and 1.
#' @param noisy Display progress to console (intermediate output, time stamps, etc.)? (logical)
#' @return bigKRLS object containing slope and uncertainty estimates
#' @examples
#'N <- 500  # proceed with caution above N = 5,000
#'k <- 4
#'X <- matrix(rnorm(N*k), ncol=k)
#'X <- cbind(X, sample(0:1, replace = TRUE, size = nrow(X)))
#'b <- runif(ncol(X))
#'y <- X %*% b + rnorm(nrow(X))
#' out <- bigKRLS(X = X, y = y)
#' @useDynLib bigKRLS
#' @importFrom Rcpp evalCpp
#' @importFrom stats pt quantile sd var
#' @importFrom utils timestamp
#' @import bigalgebra biganalytics bigmemory
#' @export
bigKRLS <- function (X = NULL, y = NULL, lambda = NULL, 
                          sigma = NULL, derivative = TRUE, binary = TRUE, vcov = TRUE, 
                          noisy = T, L = NULL, U = NULL, tol = NULL, eigtrunc = NULL) 
{
  # suppressing warnings from bigmatrix
  oldw <- getOption("warn")
  options(warn = -1)
  options(bigmemory.allow.dimnames=TRUE)
  
  if(noisy){print("starting KRLS..."); timestamp()}
  
  if(!is.big.matrix(X)){
    X <- as.big.matrix(X, type='double')
  } 
  if(is.null(colnames(X))){
    colnames(X) <- paste("x", 1:ncol(X), sep="")
  }
  miss.ind <- colna(X)
  if (sum(miss.ind) > 0) { 
    stop(paste("error: the following columns in X contain missing data, which must be removed:", 
               paste((1:length(miss.ind))[miss.ind > 0], collapse = ', '), collapse=''))
  }
  n <- nrow(X)
  d <- ncol(X)
  X.init <- X
  X.init.mean <- colmean(X)
  X.init.range <- t(colrange(X))
  X.init.sd <- colsd(X)
  if (min(X.init.sd) == 0) {
    stop(paste("error: the following columns in X are constant and must be removed:",
               which(X.init.sd == 0)))
  }
  
  if (n != nrow(y)) { stop("error: nrow(X) not equal to number of elements in y.")}
  if(!is.big.matrix(y)){
    y <- as.big.matrix(matrix(y, ncol=1), type='double')
  } 
  if (colna(y) > 0) { stop("error: y contains missing data.") }
  if (colsd(y) == 0) { stop("error: y is a constant.") }
  
  if (!is.null(eigtrunc)) {
    if (!is.numeric(eigtrunc)) 
      stop("error: eigtrunc, if used, must numeric value be between 0 and 1 indicating the proportion of eigen to be truncated.")
    if (eigtrunc > 1 | eigtrunc < 0) 
      stop("error: eigtrunc must be between 0 (no eigen truncation) and 1")
  }
  
  stopifnot(is.logical(derivative), is.logical(vcov), is.logical(binary))
  if (derivative & !vcov) {
      stop("error: vcov is needed to get derivatives (derivative==TRUE requires vcov=TRUE)")
  }
  
  if (is.null(sigma)) { # default bandwidth: sigma = ncol(X) = d
    sigma <- d          # (actually shorthand for "sigma squared") 
  }else{
    if(is.vector(sigma) | length(sigma) == 1 | is.numeric(sigma) | sigma > 0){
      print("error: sigma (bandwidth parameter) must be a positive number.")
      stop()
    }
  }
  
  x.is.binary <- apply(X, 2, function(x){length(unique(x))}) == 2 
  treat.x.as.binary <- matrix((x.is.binary + binary) == 2, nrow=1) # x is binary && user wants first differences
  # modify so that user may estimate derivatives of only some binaries?...
  colnames(treat.x.as.binary) <- colnames(X)
  
  if(noisy & sum(treat.x.as.binary) > 0){
    print(paste("The following x variable(s) is (are) binary: ",
                paste(colnames(X)[treat.x.as.binary], collapse=", "), ". ", 
                ifelse(binary, 
                       "First differences will be computed for those binaries instead of local (pairwise) derivatives (default).",
                       "Even for those binaries, local (pairwise) derivatives will be computed (not default)."),
                sep=""))
  }
  
  y.init <- deepcopy(y)
  y.init.sd <- colsd(y)
  y.init.mean <- colmean(y.init)
  
  for(i in 1:ncol(X)){
    X[,i] <- (X[,i] - mean(X[,i]))/sd(X[,i])
  }
  y[,1] <- (y[,1] - mean(y[,1]))/sd(y[,1])
  
  if(noisy){print("data successfully cleaned... getting Kernel..."); timestamp()}
  
  K <- NULL  # K is the kernel
  K <- bGaussKernel(X, sigma)
  
  if(noisy){print("got Kernel... getting Eigen"); timestamp()}
  
  Eigenobject <- bEigen(K) 
  
  if(noisy){print("got Eigen... getting regularization parameter Lambda"); timestamp()}
  
  if (is.null(lambda)) {

    lambda <- bLambdaSearch(L = L, U = U, y = y, Eigenobject = Eigenobject, 
                           eigtrunc = eigtrunc, noisy = noisy)
    
    if (noisy) {cat("Lambda that minimizes Leave-One-Out-Error (Loo) Loss is:", round(lambda, 5), "\n")}
  }else {
    stopifnot(is.vector(lambda), length(lambda) == 1, is.numeric(lambda), lambda > 0)
  }
  
  out <- bSolveForc(y = y, Eigenobject = Eigenobject, lambda = lambda, 
                   eigtrunc = eigtrunc)
  
  # bSolveForc obtains the vector of weights 
  # that assign importance to the similarity scores (found in K)

  yfitted <- K %*% matrix(out$coeffs, ncol=1)
  
  if(noisy){print("got coefficients, predicted values... getting vcovmatc"); timestamp()}
  
  if (vcov == TRUE) {
    sigmasq <- (1/n) * bCrossProd(y - yfitted)[1,1]
    
    if (is.null(eigtrunc)) {  # default
      m <- bMultDiag(Eigenobject$vectors, 
                     sigmasq * (Eigenobject$values + lambda)^-2)
      vcovmatc <- bTCrossProd(m, Eigenobject$vectors)

    }else{
      
      lastkeeper = max(which(Eigenobject$values >= eigtrunc * Eigenobject$values[1]))
      
      m <- bMultDiag(sub.big.matrix(Eigenobject$vectors, 
                                    firstCol=1, 
                                    lastCol=lastkeeper), 
                     sigmasq * (sub.big.matrix(Eigenobject$vectors, 
                                               firstCol=1, 
                                               lastCol=lastkeeper) + lambda)^-2)
      vcovmatc <- bTCrossProd(m, sub.big.matrix(Eigenobject$vectors, 
                                                firstCol=1, 
                                                lastCol=lastkeeper))
    }
    # to do: add parameter to save eigen to disk 
    remove(Eigenobject)
    gc()

    vcovmatyhat <- bCrossProd(K, vcovmatc %*% K)
  }else {
    vcov.c <- NULL
    vcov.fitted <- NULL
  }
  
  if(noisy){print("got vcovmatc...");timestamp()}  
  
  avgderiv <- varavgderivmat <- derivmat <- NULL
  
  if (derivative == TRUE) {
  
    derivmat <- big.matrix(nrow=n, ncol=d, init=NA)
    varavgderivmat <- big.matrix(nrow=d, ncol=1, init=NA)
    
    if(noisy){print("getting derivatives...."); timestamp()}
    
    for(i in 1:d){
      
      if(treat.x.as.binary[i]){
        
        if(noisy){print(paste("computing first differences for binary variable", colnames(X)[i])); timestamp()}
        X.tmp <- as.matrix(X)
        X.firstdiffs <- rbind(X.tmp, X.tmp, X.tmp)
        remove(X.tmp)
        # putting in standardized 1s and 0s...
        X.firstdiffs[1:n,i] <- (1 - X.init.mean[i])/X.init.sd[i]
        X.firstdiffs[(n + 1):(2*n), i] <- -X.init.mean[i]/X.init.sd[i]
        
        K.firstdiffs <- bTempKernel(as.big.matrix(X.firstdiffs), sigma) 
        remove(X.tmp)
        
        y.fitted.tmp <- as.matrix(K.firstdiffs %*% as.matrix(out$coeffs))
        derivmat[,i] <- (y.fitted.tmp[1:n] - y.fitted.tmp[(n + 1):(2*n)])*X.init.sd[i]
        
        remove(y.fitted.tmp)
        
        vcov.fitted <- bTCrossProd(K.firstdiffs %*% vcovmatc, K.firstdiffs)
        remove(K.firstdiffs)
        
        h <- matrix(rep(c(1/n, -(1/n)), each = n), ncol = 1) 
        # vector of length 2n, first half 1/n, second half -1/n
        h <- as.big.matrix(h)
        varavgderivmat[i] <- 2 * X.init.sd[i]^2 * as.matrix((bCrossProd(h, vcov.fitted) %*% h))
        remove(vcov.fitted)
        remove(h)
        gc()
      }else{
        
        if(noisy){print(paste("computing local derivatives of", colnames(X)[i])); timestamp()}
        x <- as.matrix(X[,i])
        distancek <- apply(x, 1, function(xi, x){x - xi}, x)  # see if dist, lower.tri can expedite (symmetric)
        remove(x)
        L <- bElementwise(distancek, K)
        remove(distancek)
        gc()
        
        derivmat[,i] <- (-2/sigma) * as.matrix(L %*% as.big.matrix(out$coeff))
        
        if(noisy){print(paste("computing variance of local derivatives of", colnames(X)[i])); timestamp()}
        
        varavgderivmat[i] <- (1/n^2) * (-2/sigma)^2 * sum(bCrossProd(L, vcovmatc %*% L))
        remove(L)
        gc()
      }
    }
    
    if(noisy){print("finished major calculations; rescaling, etc..."); timestamp()}
    
    derivmat <- as.matrix(derivmat)   # added june 14
    derivmat <- y.init.sd * derivmat
    for(i in 1:ncol(derivmat)){
      derivmat[,i] <- derivmat[,i]/X.init.sd[i]
    }
    
    attr(derivmat, "scaled:scale") <- NULL
    colnames(derivmat) <- colnames(X)
    # avgderiv <- matrix(colmean(derivmat), nrow=1)
    avgderiv <- colMeans(derivmat)
    attr(avgderiv, "scaled:scale") <- NULL
    varavgderivmat <- matrix((y.init.sd/X.init.sd)^2 * as.matrix(varavgderivmat), nrow=1)
    colnames(varavgderivmat) <- colnames(X)
    attr(varavgderivmat, "scaled:scale") <- NULL
  }
  
  yfitted <- as.matrix(yfitted) * y.init.sd + y.init.mean
  
  if (vcov == TRUE) {
    vcov.c <- (y.init.sd^2) * vcovmatc
    vcov.fitted <- (y.init.sd^2) * vcovmatyhat
  }else {
    vcov.c <- NULL
    vcov.fitted <- NULL
  }
  Looe <- out$Le * y.init.sd
  R2 <- 1 - (var(y.init - yfitted)/(y.init.sd^2))
  q <- list(K = K, coeffs = out$coeffs, Looe = Looe, fitted = yfitted, 
            X = X.init, y = y.init, sigma = sigma, lambda = lambda, 
            R2 = R2, derivatives = as.matrix(derivmat), avgderivatives = avgderiv, 
            var.avgderivatives = varavgderivmat, vcov.c = vcov.c, 
            vcov.fitted = vcov.fitted, binaryindicator = treat.x.as.binary)
  class(q) <- "bigKRLS" 
  
  if (noisy && derivative) {
    cat("Average Marginal Effects: \n")
    print(round(q$avgderivatives, 3))
    cat("\n Percentiles of Local Derivatives: \n")
    print(round(apply(as.matrix(q$derivatives), 2, quantile, 
                probs = c(0.05, 0.25, 0.5, 0.75, 0.95)),3))
    cat("\n For more detail, use summary() on the outputted object.")
  }

  return(q)
  
  options(warn = oldw)
}  

#' @export
bLambdaSearch <- function (L = NULL, U = NULL, y = NULL, Eigenobject = NULL, tol = NULL, 
                          noisy = FALSE, eigtrunc = NULL){
  n <- nrow(y)
  if (is.null(tol)) {
    tol <- 10^-3 * n # tolerance parameter
  } else {
    stopifnot(is.vector(tol), length(tol) == 1, is.numeric(tol), tol > 0)
  }
  if (is.null(U)) {
    U <- n
    while (sum(Eigenobject$values/(Eigenobject$values + U)) < 1) {
      U <- U - 1
    }
  } else {
    stopifnot(is.vector(U), length(U) == 1, is.numeric(U), U > 0)
  }
  if (is.null(L)) {
    q <- which.min(abs((Eigenobject$values - max(Eigenobject$values)/1000)))

    L = .Machine$double.eps
    
    # .Machine is a variable holding information on the numerical characteristics 
    # of the machine R is running on, 
    # such as the largest double or integer and the machine's precision.
    
    # double.eps is the smallest positive floating-point number x 
    # such that 1 + x != 1. Normally 2.220446e-16.
    
    while (sum(Eigenobject$values/(Eigenobject$values + L)) > q) {
      L <- L + 0.05 # L starts near 0 so sum(Eigenobject$values/(Eigenobject$values + L)) 
      # ~ sum(Eigenobject$values/Eigenobject$values) = N
      # as L increases, the sum of the ratios decreases quickly
      # since jumping from L = 0 to L = 0.05 leapfrogs many small eigenvalues...
    } 
  } else {
    stopifnot(is.vector(L), length(L) == 1, is.numeric(L), L >= 0)
  }
  X1 <- L + (0.381966) * (U - L) 
  X2 <- U - (0.381966) * (U - L)

  # bLooLoss is big Leave One Out Error Loss
  
  S1 <- bLooLoss(lambda = X1, y = y, Eigenobject = Eigenobject, 
                eigtrunc = eigtrunc)
  S2 <- bLooLoss(lambda = X2, y = y, Eigenobject = Eigenobject, 
                eigtrunc = eigtrunc)
  if (noisy) {
    cat("L:", L, "X1:", X1, "X2:", X2, "U:", U, "S1:", S1, 
        "S2:", S2, "\n")
  }
  while (abs(S1 - S2) > tol) {
    if (S1 < S2) {
      U <- X2
      X2 <- X1
      X1 <- L + (0.381966) * (U - L)
      S2 <- S1
      S1 <- bLooLoss(lambda = X1, y = y, Eigenobject = Eigenobject, 
                    eigtrunc = eigtrunc)
    }
    else {
      L <- X1
      X1 <- X2
      X2 <- U - (0.381966) * (U - L)
      S1 <- S2
      S2 <- bLooLoss(lambda = X2, y = y, Eigenobject = Eigenobject, 
                    eigtrunc = eigtrunc)
    }
    if (noisy) {
      cat("L:", L, "X1:", X1, "X2:", X2, "U:", U, "S1:", 
          S1, "S2:", S2, "\n")
    }
  }
  out <- ifelse(S1 < S2, X1, X2)
  if (noisy) {
    cat("Lambda:", out, "\n")
  }
  return(invisible(out))
}

#' @export
bSolveForc <- function (y = NULL, Eigenobject = NULL, lambda = NULL, eigtrunc = NULL) {
  nn <- nrow(y)
  if (is.null(eigtrunc)) {
    #split this line into two separate operations
    #the multdiag() call takes a trivial amount of time, so no need to parallelize
    
    m <- bMultDiag(Eigenobject$vectors, 
                   1/(Eigenobject$values + lambda))
    Ginv <- bTCrossProd(m, Eigenobject$vectors)
    
  }else {
    lastkeeper = max(which(Eigenobject$values >= eigtrunc * 
                             Eigenobject$values[1]))
    lastkeeper = max(1, lastkeeper)
    
    m <- bMultDiag(Eigenobject$vectors,
                   1/(Eigenobject$values[1:lastkeeper] + lambda),
                   (1:lastkeeper))
    Ginv <- bTCrossProd(m, Eigenobject$vectors, 1:lastkeeper)
  }
  
  coeffs <- (Ginv %*% y)[,]
  Le <- crossprod(coeffs/bDiag(Ginv))
  return(list(coeffs = coeffs, Le = Le))
}

#' @export
bLooLoss <- function (y = NULL, Eigenobject = NULL, lambda = NULL, eigtrunc = NULL) 
{
  return(bSolveForc(y = y, Eigenobject = Eigenobject, lambda = lambda, 
                   eigtrunc = eigtrunc)$Le)
} # not sure that there's any point to this function
# could just make "bLooLoss" mode a parameter of bSolveForc

#' @export
predict.bigKRLS <- function (object, newdata, se.fit = FALSE, ...) 
{
  if (class(object) != "bigKRLS") {
    warning("Object not of class 'bigKRLS'")
    UseMethod("predict")
    return(invisible(NULL))
  }
  if (se.fit == TRUE) {
    if (is.null(object$vcov.c)) {
      stop("recompute bigKRLS object with bigKRLS(,vcov=TRUE) to compute standard errors")
    }
  }
  newdata <- as.big.matrix(newdata)
  if (ncol(object$X) != ncol(newdata)) {
    stop("ncol(newdata) differs from ncol(X) from fitted krls object")
  }
  Xmeans <- colmean(object$X)
  Xsd <- colsd(object$X)
  
  for(i in 1:ncol(X)){
    X[,i] <- (X[,i] - mean(X[,i]))/sd(X[,i])
  }  
  
  newdata.init <- newdata
  
  for(i in 1:ncol(newdata)){
    newdata[,i] <- (newdata[,i] - mean(newdata[,i]))/sd(newdata[,i])
  }
  
  nn <- nrow(newdata)
  newdata.X <- big.matrix(nrow=(nn + nrow(X)), 
                          ncol=(ncol(X)),
                          init=NA)
  for(i in 1:nn){
    newdata.X[,i] <- newdata[,i]
  }
  for(j in (nn+1):(nn+nrow(X))){
    newdata.X[,j] <- X[,j]
  }
  
  newdataK <- bGaussKernel(newdata.X, object$sigma)
  newdataK <- sub.big.matrix(newdataK, 
                             firstRow = 1, lastRow = nn,
                             firstCol = nn, lastCol=nrow(X))
  
  yfitted <- newdataK %*% as.matrix(object$coeffs, nrow=1)
  if (se.fit) {
    vcov.c.raw <- object$vcov.c * as.vector((1/var(object$y)))
    vcov.fitted <- bTCrossProd(newdataK %*% vcov.c.raw, newdataK)
    vcov.fit <- var(object$y) * vcov.fitted
    se.fit <- matrix(sqrt(diag(vcov.fit)), ncol = 1)
  }
  else {
    vcov.fit <- se.fit <- NULL
  }
  yfitted <- (yfitted * sd(object$y) + mean(object$y))
  return(list(fit = yfitted, se.fit = se.fit, vcov.fit = vcov.fit, 
              newdata = newdata, newdataK = newdataK))
}

#' @export
summary.bigKRLS <- function (object, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), digits=3,...) 
{
  if (class(object) != "bigKRLS") {
    warning("Object not of class 'bigKRLS'")
    UseMethod("summary")
    return(invisible(NULL))
  }
  cat("* *********************** *\n")
  cat("Model Summary:\n\n")
  cat("R2:", object$R2, "\n\n")
  d <- ncol(object$X)
  n <- nrow(object$X)
  coefficients <- matrix(NA, d, 0)
  rownames(coefficients) <- colnames(object$X)
  if (is.null(object$derivatives)) {
    cat("\n")
    cat("recompute krls object with krls(...,derivative = TRUE) to get summary of marginal effects\n")
    return(invisible(NULL))
  }
  est <- t(object$avgderivatives)
  se <- sqrt(object$var.avgderivatives)
  tval <- est/se
  avgcoefficients <- cbind(t(est), t(se), t(tval), t(2 * pt(abs(tval), 
                                                 n - d, lower.tail = FALSE)))
  colnames(avgcoefficients) <- c("Est", "Std. Error", "t value", "Pr(>|t|)")
  if (sum(object$binaryindicator) > 0) {
    rownames(avgcoefficients)[object$binaryindicator] <- paste(rownames(avgcoefficients)[object$binaryindicator], 
                                                               "*", sep = "")
  }
  cat("Average Marginal Effects:\n")
  print(round(avgcoefficients, digits))
  if (sum(object$binaryindicator) > 0) {
    cat("\n(*) average dy/dx is for discrete change of dummy variable from min to max (i.e. usually 0 to 1))\n\n")
  }
  qderiv <- apply(object$derivatives, 2, quantile, probs = probs)
  if (sum(object$binaryindicator) > 0) {
    colnames(qderiv)[object$binaryindicator] <- paste(colnames(qderiv)[object$binaryindicator], 
                                                      "*", sep = "")
  }
  qderiv <- t(qderiv)
  cat("\n")
  cat("Percentiles of Local Derivatives:\n")
  print(round(qderiv, digits))
  if (sum(object$binaryindicator) > 0) {
    cat("\n(*) quantiles of dy/dx is for discrete change of dummy variable from min to max (i.e. usually 0 to 1))\n\n")
  }
  ans <- list(coefficients = avgcoefficients, 
              qcoefficients = qderiv)
  class(ans) <- "summary.bigKRLS"
  return(invisible(ans))
}

##################
# Rcpp Functions #
##################

#' @export
bMultDiag <- function (X, d) {
  #rcpp_multdiag.cpp
  out <- big.matrix(nrow=nrow(X),
                    ncol=ncol(X),
                    init=0,
                    type='double')
  d <- as.big.matrix(matrix(d, nrow=1))
  
  BigMultDiag(X@address, d@address, out@address)
  
  return(out)
}

#' @export
bEigen <- function(X){
  #rcpp_eigen.cpp
  vals <- big.matrix(nrow = 1,
                     ncol = ncol(X),
                     init = 0,
                     type = 'double')
  vecs <- big.matrix(nrow = nrow(X),
                     ncol = ncol(X),
                     init=0,
                     type='double')
  BigEigen(X@address, vals@address, vecs@address)
  return(list('values' = vals[,], 'vectors' = vecs))
}

#' @export
bGaussKernel <- function(X, sigma){
  #rcpp_gauss_kernel.cpp
  out <- big.matrix(nrow=nrow(X), ncol=nrow(X), init=0)
  s <- big.matrix(1,1,type='double', init=sigma)
  
  BigGaussKernel(X@address, out@address, s@address)
  return(out)
}

#' @export
bTempKernel <- function(X, sigma){
  #rcpp_temp_kernel.cpp
  n <- nrow(X)/3
  out <- big.matrix(nrow=(n*2), ncol=n, init=0)
  s <- big.matrix(1,1,type='double', init=sigma)
  
  BigTempKernel(X@address, out@address, s@address)
  return(out)
}

#' @export
bCrossProd <- function(X,Y=NULL){
  if(is.null(Y)){
    Y <- deepcopy(X)
  }
  out <- big.matrix(nrow = ncol(X),
                    ncol = ncol(Y),
                    init = 0,
                    type = 'double')
  BigCrossProd(X@address, Y@address, out@address)
  return(out)
}

#' @export
bTCrossProd <- function(X,Y=NULL){
  if(is.null(Y)){
    Y <- deepcopy(X)
  }
  out <- big.matrix(nrow = nrow(X),
                    ncol = nrow(Y),
                    init = 0,
                    type = 'double')
  BigTCrossProd(X@address, Y@address, out@address)
  return(out)
}

#' @export
bDiag <- function(X){
  # return the diagonal elements of a bigmatrix
  out <- sapply(1:nrow(X), function(i){X[i,i]})
  return(out)
}

#' @export
bElementwise <- function(X,Y=NULL){
  if(!is.big.matrix(X)){
    X <- as.big.matrix(X)
  }
  if(is.null(Y)){
    Y <- deepcopy(X)
  }
  
  out <- big.matrix(nrow=nrow(X), ncol=ncol(X), init=0)
  
  BigElementwise(X@address, Y@address, out@address)
  
  return(out)
}