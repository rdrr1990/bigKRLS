#' Kernel Regularized Least Squares with Big Matrices
#' 
#' @param y A vector of observations on the dependent variable; missing values not allowed. May be base R matrix or library(bigmemory) big.matrix.
#' @param X A matrix of observations of the independent variables; factors, missing values, and constant vectors not allowed. May be base R matrix or library(bigmemory) big.matrix.
#' @param sigma Bandwidth parameter, shorthand for sigma squared. Default: sigma <- ncol(X). Since x variables are standardized, facilitates interprepation of the Gaussian kernel, exp(-dist(X)^2/sigma) a.k.a the similarity score. Of course, if dist between observation i and j is 0, there similarity is 1 since exp(0) = 1. Suppose i and j differ by one standard deviation on each dimension. Then the similarity is exp(-ncol(X)/sigma) = exp(-1) = 0.368.  
#' @param derivative Logical: Estimate derivatives (as opposed to just coefficients)? Recommended for interpretability.
#' @param binary Logical: Estimate first differences for each variable with only two observed outcomes?
#' @param vcov.est Logical: Estimate variance covariance matrix? Required to obtain derivatives and standard errors on predictions (default = TRUE).
#' @param lambda Regularization parameter. Default: estimated based (in part) on the eigenvalues of the kernel via Golden Search with convergence parameter "tolerance." Must be positive, real number. 
#' @param L Lower bound of Golden Search for lambda. 
#' @param U Upper bound of Golden Search for lambda.
#' @param tol tolerance parameter for Golden Search for lambda. Default: N / 1000.
#' @param eigtrunc Number of eigenvectors and values to be use. Must be between 0 (no truncation) and N. Eigentruncation may increase speed but may reduce precision of the regularization parameter and therefore the other estimates.
#' @param noisy Logical: Display progress to console (intermediate output, time stamps, etc.)? (bigKRLS runs a touch faster with noisy = FALSE but it is recommended until you have a sense of how it runs on your system, with your data, etc.)
#' @return bigKRLS Object containing slope and uncertainty estimates; summary and predict defined for class bigKRLS.
#' @examples
#'N <- 500  # proceed with caution above N = 10,000 for system with 8 gigs made avaiable to R
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
#' @import bigalgebra biganalytics bigmemory shiny
#' @export
bigKRLS <- function (y = NULL, X = NULL, sigma = NULL, derivative = TRUE, binary = TRUE, vcov.est = TRUE, 
                     lambda = NULL, L = NULL, U = NULL, tol = NULL, eigtrunc = NULL, noisy = TRUE)
{
  if(noisy){cat("starting KRLS... \n\nvalidating inputs, prepping data, etc... \n")}
  
  # suppressing warnings from bigmatrix
  oldw <- getOption("warn")
  options(warn = -1)
  options(bigmemory.allow.dimnames=TRUE)
  
  if(noisy){
    if(is.big.matrix(X)){
      cat('Input given as big.matrix object; X and derivatives will be returned as bigmatrix objects. Be sure to use save.bigKRLS() to store results!\n')
    }else{
      cat('Input given as base R matrices; X and derivatives will be returned as base R matrices.\n')
    }
  } 
  
  if(!is.big.matrix(X) & nrow(X) < 2500){
    big.matrix.in <- FALSE
    
  } else{
    if(noisy == T){
      cat('Input given as a bigmatrix object or N > 2,500. N x N matrices will be returned as bigmatrices. Be sure to use save.bigKRLS() to store results! \n')
    }
    big.matrix.in <- TRUE
  }
  
  X <- to.big.matrix(X)
  y <- to.big.matrix(y, d=1)
  
  if(is.null(colnames(X))){
    colnames(X) <- paste("x", 1:ncol(X), sep="")
  }
  colnames(X)[which(apply(as.matrix(colnames(X)), 1, nchar) == 0)] <- paste("x", which(apply(as.matrix(colnames(X)), 1, nchar) == 0), sep="")
  miss.ind <- colna(X)
  if (sum(miss.ind) > 0) { 
    stop(paste("the following columns in X contain missing data, which must be removed:", 
               paste((1:length(miss.ind))[miss.ind > 0], collapse = ', '), collapse=''))
  }
  n <- nrow(X)
  d <- ncol(X)
  
  X.init <- deepcopy(X)
  X.init.sd <- colsd(X)
  
  if (min(X.init.sd) == 0) {
    stop(paste("the following columns in X are constant and must be removed:",
               which(X.init.sd == 0)))
  }
  
  if (n != nrow(y)) { stop("nrow(X) not equal to number of elements in y.")}
  if (colna(y) > 0) { stop("y contains missing data.") }
  if (colsd(y) == 0) { stop("y is a constant.") }
  
  if(!is.null(lambda)){stopifnot(is.vector(lambda), length(lambda) == 1, is.numeric(lambda), lambda > 0)}
  
  if(!is.null(sigma)){stopifnot(is.vector(sigma), length(sigma) == 1, is.numeric(sigma), sigma > 0)}
  sigma <- ifelse(is.null(sigma), d, sigma)
  
  if (is.null(tol)) { # tolerance parameter for lambda search
    tol <- n/1000 
  } else {
    stopifnot(is.vector(tol), length(tol) == 1, is.numeric(tol), tol > 0)
  }
  
  if (!is.null(eigtrunc) && (!is.numeric(eigtrunc) | eigtrunc > n | eigtrunc < 0)) {
    stop("eigtrunc, if used, must be a number between 0 and N indicating the number of eigenvalues to be used.")
  }
  
  stopifnot(is.logical(derivative), is.logical(vcov.est), is.logical(binary))
  if (derivative & !vcov.est) { stop("vcov.est is needed to get derivatives (derivative==TRUE requires vcov.est=TRUE)")}

  x.is.binary <- apply(X, 2, function(x){length(unique(x))}) == 2 
  treat.x.as.binary <- matrix((x.is.binary + binary) == 2, nrow=1) # x is binary && user wants first differences
  colnames(treat.x.as.binary) <- colnames(X)
  
  if(noisy & sum(treat.x.as.binary > 0)){
    cat(paste("Since binary == ", binary, ",", 
        ifelse(binary, " first differences will be computed for: ", 
               "first differences will NOT be computed for: "), 
        paste(colnames(X)[treat.x.as.binary], collapse=", "), ".", sep=""))
  }
  
  y.init <- deepcopy(y)
  y.init.sd <- colsd(y.init)
  y.init.mean <- colmean(y.init)
  
  for(i in 1:ncol(X)){
    X[,i] <- (X[,i] - mean(X[,i]))/sd(X[,i])
  }
  y[,1] <- (y[,1] - mean(y[,1]))/sd(y[,1])
  
  if(noisy){cat("\ndata successfully cleaned...\n\nstep 1/5: getting Kernel...\n"); timestamp()}
  
  K <- NULL  # K is the kernel
  K <- bGaussKernel(X, sigma)
  
  if(noisy){cat("\nstep 2/5: getting Eigenvectors and values...\n"); timestamp()}
  
  Eigenobject <- bEigen(K, eigtrunc) 
  
  if(noisy){cat("\nstep 3/5: getting regularization parameter Lambda which minimizes Leave-One-Out-Error Loss via Golden Search...\n"); timestamp(); cat("\n\n")}
  
  if (is.null(lambda)) {
    lambda <- bLambdaSearch(L = L, U = U, y = y, Eigenobject = Eigenobject, eigtrunc = eigtrunc, noisy = noisy)
  }
  
  out <- bSolveForc(y = y, Eigenobject = Eigenobject, lambda = lambda, eigtrunc = eigtrunc)
  
  # bSolveForc obtains the vector of coefficients (weights) 
  # that assign importance to the similarity scores (found in K)

  yfitted <- K %*% matrix(out$coeffs, ncol=1)
  
  if(noisy){cat("\nstep 4/5: getting coefficients & fitted values...\n"); timestamp()}
  
  if (vcov.est == TRUE) {
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
    remove(Eigenobject)
    remove(m)
    gc()

    vcovmatyhat <- bCrossProd(K, vcovmatc %*% K)
  }else {
    vcov.est.c <- NULL
    vcov.est.fitted <- NULL
  }
  
  if(noisy){cat("\nstep 5/5: getting derivatives...\n\t>>> run time proportional to ncol(X) and nrow(X)^2...\n\n");timestamp()}  
  
  if (derivative == TRUE) {
  
    deriv_out <- bDerivatives(X,sigma,K,out$coeffs,vcovmatc,X.init.sd)
    
    derivmat <- deriv_out$derivatives
    varavgderivmat <- deriv_out$varavgderiv
    
    if(noisy){cat("finished major calculations; rescaling, etc...\n"); timestamp()}
    
    derivmat <- y.init.sd * derivmat
    for(i in 1:ncol(derivmat)){
      derivmat[,i] <- derivmat[,i]/X.init.sd[i]
    }
    
    attr(derivmat, "scaled:scale") <- NULL
    avgderiv <- matrix(colmean(derivmat), nrow=1)
    attr(avgderiv, "scaled:scale") <- NULL
    
    varavgderivmat <- matrix((y.init.sd/X.init.sd)^2 * as.matrix(varavgderivmat), nrow=1)
    attr(varavgderivmat, "scaled:scale") <- NULL
  }
  
  yfitted <- as.matrix(yfitted) * y.init.sd + y.init.mean
  
  if (vcov.est == TRUE) {
    vcov.est.c <- (y.init.sd^2) * vcovmatc
    vcov.est.fitted <- (y.init.sd^2) * vcovmatyhat
  }else {
    vcov.est.c <- NULL
    vcov.est.fitted <- NULL
  }
  Looe <- out$Le * y.init.sd
  R2 <- 1 - (var(y.init - yfitted)/(y.init.sd^2))
  R2AME <- cor(y.init[,], (X %*% matrix(avgderiv, ncol=1))[,])^2
  # Pseudo R2 using only Average Marginal Effects
  
  # return estimates as base R matrices if inputted as base R (regular) matrices
  if(!big.matrix.in){
    X.init <- X.init[]
    derivmat <- derivmat[]
    if(n < 2500){ # but only for small N for the N*N matrices 
      K <- K[]
      vcov.est.c <- vcov.est.c[]
      vcov.est.fitted <- vcov.est.fitted[]
    }
  }
  
  # y.init is just a vector so can be safely returned as a base R object
  y.init <- y.init[]
  
  w <- list(K = K, coeffs = out$coeffs, Looe = Looe, fitted = yfitted, 
            X = X.init, y = y.init, sigma = sigma, lambda = lambda, 
            R2 = R2, derivatives = derivmat, avgderivatives = avgderiv, 
            var.avgderivatives = varavgderivmat, vcov.est.c = vcov.est.c,
            vcov.est.fitted = vcov.est.fitted, binaryindicator = treat.x.as.binary, R2AME=R2AME)
  
  colnames(w$derivatives) <- colnames(w$avgderivatives) <- colnames(X.init)
  class(w) <- "bigKRLS" 
  
  if (noisy && derivative) {
    cat("Average Marginal Effects: \n")
    print(round(w$avgderivatives, 3))
    cat("\n Percentiles of Local Derivatives: \n")
    print(round(apply(as.matrix(w$derivatives), 2, quantile, 
                probs = c(0.25, 0.5, 0.75)),3))
    cat("\n For more detail, use summary() on the outputted object. Use save.bigKRLS() to store results.")
  }

  return(w)
  
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
  
  if (noisy) {cat("Lambda:", round(out, 5), "\n")}
  
  return(invisible(out))
}

#' @export
bSolveForc <- function (y = NULL, Eigenobject = NULL, lambda = NULL, eigtrunc = NULL) {
  nn <- nrow(y)
  if (is.null(eigtrunc)) {
    m <- bMultDiag(Eigenobject$vectors, 
                   1/(Eigenobject$values + lambda))
    Ginv <- bTCrossProd(m, Eigenobject$vectors)
    
    rm(m)
    gc()
    
  }else {
    lastkeeper = max(which(Eigenobject$values >= eigtrunc * 
                             Eigenobject$values[1]))
    lastkeeper = max(1, lastkeeper)
    
    m <- bMultDiag(Eigenobject$vectors,
                   1/(Eigenobject$values[1:lastkeeper] + lambda),
                   (1:lastkeeper))
    Ginv <- bTCrossProd(m, Eigenobject$vectors, 1:lastkeeper)
    
    rm(m)
    gc()
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
  if(se.fit == TRUE) {
    if (is.null(object$vcov.est.c)) {
      stop("recompute bigKRLS object with bigKRLS(,vcov.est=TRUE) to compute standard errors")
    }
  }
  
  # convert everything to a bigmatrix for internal usage
  object$X <- to.big.matrix(object$X)
  object$K <- to.big.matrix(object$K)
  object$derivatives <- to.big.matrix(object$derivatives)
  object$vcov.est.c <- to.big.matrix(object$vcov.est.c)
  object$vcov.est.fitted <- to.big.matrix(object$vcov.est.fitted)
  
  # set bigmatrix flag for input data for later
  if(!is.big.matrix(newdata)){
    bigmatrix.in <- FALSE
  } else{
    bigmatrix.in <- TRUE
  }
  
  newdata <- to.big.matrix(newdata)
  
  if (ncol(object$X) != ncol(newdata)) {
    stop("ncol(newdata) differs from ncol(X) from fitted bigKRLS object")
  }
  Xmeans <- colmean(object$X)
  Xsd <- colsd(object$X)
  
  for(i in 1:ncol(object$X)){
    object$X[,i] <- (object$X[,i] - Xmeans[i])/Xsd[i]
  }  
  
  newdata.init <- newdata
  
  for(i in 1:ncol(newdata)){
    newdata[,i] <- (newdata[,i] - Xmeans[i])/Xsd[i]
  }
  
  newdataK <- bTempKernel(newdata, object$X, object$sigma)
  
  # convert to regular matrix
  yfitted <- (newdataK %*% as.matrix(object$coeffs, ncol=1))[]
  
  if (se.fit) {
    vcov.est.c.raw <- object$vcov.est.c * (1/var(object$y))
    vcov.est.fitted <- bTCrossProd(newdataK %*% vcov.est.c.raw, newdataK)
    vcov.est.fit <- var(object$y) * vcov.est.fitted
    se.fit <- matrix(sqrt(diag(vcov.est.fit[])), ncol = 1)
  }
  else {
    vcov.est.fit <- se.fit <- NULL
  }
  
  yfitted <- (yfitted * sd(object$y) + mean(object$y))
  
  
  
  if(!bigmatrix.in){
    newdata <- newdata[]
    vcov.est.fit <- vcov.est.fit[]
    newdataK <- newdataK[]
  }
  
  return(list(fit = yfitted, se.fit = se.fit, vcov.est.fit = vcov.est.fit, 
              newdata = newdata, newdataK = newdataK))
}

#' @export
summary.bigKRLS <- function (object, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), digits=4,...) 
{
  if (class(object) != "bigKRLS") {
    warning("Object not of class 'bigKRLS'")
    UseMethod("summary")
    return(invisible(NULL))
  }
  cat("* *********************** *\n")
  cat("Model Summary:\n\n")
  cat("R2:", round(object$R2, digits), "\n")
  cat("R2AME**:", round(object$R2AME, digits), "\n\n")
  d <- ncol(object$X)
  n <- nrow(object$X)
  coefficients <- matrix(NA, d)
  rownames(coefficients) <- colnames(object$X)
  if (is.null(object$derivatives)) {
    cat("\n")
    cat("recompute with bigKRLS(..., derivative = TRUE) for summary of marginal effects\n")
    return(invisible(NULL))
  }
  est <- object$avgderivatives
  se <- sqrt(object$var.avgderivatives)
  tval <- est/se
  pval <- 2 * pt(abs(tval), n - d, lower.tail = FALSE)
  avgcoefficients <- t(rbind(est, se, tval, pval))
  colnames(avgcoefficients) <- c("Est", "Std. Error", "t value", "Pr(>|t|)")
  rownames(avgcoefficients) <- colnames(object$X)
  if (sum(object$binaryindicator) > 0) {
    rownames(avgcoefficients)[object$binaryindicator] <- paste(rownames(avgcoefficients)[object$binaryindicator], 
                                                               "*", sep = "")
  }
  cat("Average Marginal Effects:\n")
  print(round(avgcoefficients, digits))
  qderiv <- apply(object$derivatives, 2, quantile, probs = probs)
  colnames(qderiv) <- colnames(object$X)
  if (sum(object$binaryindicator) > 0) {
    colnames(qderiv)[object$binaryindicator] <- paste(colnames(qderiv)[object$binaryindicator], 
                                                      "*", sep = "")
  }
  qderiv <- t(qderiv)
  cat("\n")
  cat("Percentiles of Local Derivatives:\n")
  print(round(qderiv, digits))
  if (sum(object$binaryindicator) > 0) {
    cat("\n(*) Reported average and quantiles of dy/dx is for discrete change of the dummy variable from min to max (usually 0 to 1)).\n\n")
  }
  cat("\n(**) Pseudo-R^2 computed only using the Average Marginal Effects.\n\n")
  ans <- list(coefficients = avgcoefficients, 
              qcoefficients = qderiv)
  class(ans) <- "summary.bigKRLS"
  return(invisible(ans))
}

#' @export
save.bigKRLS <- function (object, model_subfolder_name) 
{
  if (class(object) != "bigKRLS") {
    warning("Object not of class 'bigKRLS'")
    UseMethod("save")
    return(invisible(NULL))
  }
  stopifnot(is.character(model_subfolder_name))
  
  dir.create(model_subfolder_name)
  wd.original <- getwd()
  setwd(paste(c(wd.original, .Platform$file.sep, model_subfolder_name), collapse=""))
  cat("Saving model estimates to:\n\n", getwd(), "\n\n")
  
  bigKRLS_out <- object  
  for(i in 1:length(object)){
    if(is.big.matrix(object[[i]])){
      cat("\twriting", paste(c(names(object)[i], ".txt"), collapse = ""), "...\n")
      write.big.matrix(x = object[[i]], 
                       filename = paste(c(names(object)[i], ".txt"), collapse = ""))
      bigKRLS_out[[i]] <- NULL
    }
  }
  
  Nbm <- sum(lapply(object, class) == "big.matrix")
  cat("\n", Nbm, " matrices saved as big matrices", 
      ifelse(Nbm == 0, " (base R save() may be used safely in this case too).\n",
             ", which should be loaded back into R with bigmemory::read.big.matrix()\nSmaller outputted objects saved in estimates.rdata. \n"), 
      "Total file size approximately ", round(sum(file.info(list.files())$size)/1024^2), " megabytes.",
      sep="")
    
  save(bigKRLS_out, file="estimates.rdata")
  setwd(wd.original) 
  
}


#' @export
to.big.matrix <- function(obj, d=NULL){
  if(is.null(d)){
    d <- ifelse(!is.null(ncol(obj)), ncol(obj), 1)
  }
  
  if(!is.big.matrix(obj)){
    obj <- as.big.matrix(matrix(obj, ncol=d))
  }
  return(obj)
}

#' @export
shiny.bigKRLS <- function(out, export=F, main.label = NULL, plot.main.label = NULL, labs = NULL,
                          shiny.palette = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                                            "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")){
  
  if(!export){cat("export set to false; set export to true to prepare files for server or other machine.")}
  if(!is.null(labs)){
    colnames(out$X) <- colnames(out$derivatives) <- names(out$avgderivatives) <- names(out$var.avgderivatives) <- labs
  }
  
  palette(shiny.palette)
  
  bigKRLS_server <- shinyServer(function(input, output, session) {
    
    selectedData <- reactive({
      return(list(cbind(out$derivatives[, input$dydxp], 
                        out$X[, c(input$xp)]), input$type))
    })
    
    output$graph <- renderPlot({
      
      if(selectedData()[[2]] == "Smooth"){
        
        L <- loess.smooth(x=selectedData()[[1]][,2], 
                          y=selectedData()[[1]][,1])
        
        plot(y=L$y, x=L$x, ylab=paste("Marginal Effect of",input$dydxp), pch = 19, bty = "n",
             main=plot.main.label, 
             xlab=paste("Observed Level of", input$xp), cex=2, cex.axis=1.5,  cex.lab=1.4,
             col = colorRampPalette(c("blue", "red"))(length(L$y))[rank(L$y)])
        
      }else{
        plot(x=(selectedData()[[1]][,2]), xlab = paste("Observed Level of", input$xp),
             y=(selectedData()[[1]][,1]), ylab = paste("Marginal Effect of",input$dydxp), 
             pch = 4, bty = "n", cex=2, cex.axis=1.5,  cex.lab=1.4,
             main=plot.main.label,
             col = colorRampPalette(c("green", "purple"))(nrow(out$X))[rank(out$coeffs^2)], 
             ylim = range(selectedData()[[1]][,1])*c(.8, 1.25), 
             xlim = range(selectedData()[[1]][,2])*c(.8, 1.25))
        
        fields::image.plot(legend.only = T, zlim=c(1/nrow(out$X), 1), 
                           legend.cex = 0.75,legend.shrink = .4,   
                           col = colorRampPalette(c("purple", "green"))(nrow(out$X)))
        text(x = 1.2*range(selectedData()[[1]][,2])[2], 
             y = .75*range(selectedData()[[1]][,1])[2], 
             "Relative Fit \nIn Full Model") 
      }
    })})
  
  bigKRLS_ui <- shinyUI(fluidPage(
    
    titlePanel(main.label),
    
    sidebarPanel(
      selectInput('dydxp', 'Local Derivatives (dy/dx)', colnames(out$derivatives)),
      selectInput('xp', 'x', colnames(out$X)), 
      selectInput('type', 'Plot Type', c("Smooth", "Scatter"))
    ),
    
    mainPanel(plotOutput('graph'))
    
  ))
  
  if(export){
    
    out <- out
    out$K <- tmp$vcov.c <- tmp$vcov.fitted <- NULL
    for(i in which(unlist(lapply(out, is.big.matrix)))){
      out[[i]] <- as.matrix(out[[i]])
    }
    
    save(out, file="shiny_out.rdata")
    
    cat("A re-formatted version of your output has been saved with file name \"shiny_out.rdata\" in your current working directory:\n", getwd(),
        "\nFor a few technical reasons, the big N * N matrices have been removed and the smaller ones converted back to base R;\nthis should make your output small enough for the free version of Shiny's server.\nTo access the Shiny app later or on a different machine, simply execute this script with the following commands:\n",
        "\nload(\"shiny_out.rdata\")\nNext, execute this script to make sure Shiny is initialized with current values. \nshiny_bigKRLS(out)")
  }else{
    shinyApp(ui = bigKRLS_ui, server = bigKRLS_server)
  }
}


##################
# Rcpp Functions #
##################

#' @export
bMultDiag <- function (X, v) {
  #rcpp_multdiag.cpp
  out <- big.matrix(nrow=nrow(X),
                    ncol=ncol(X),
                    init=0,
                    type='double')
  v <- to.big.matrix(v, d=1)
  
  BigMultDiag(X@address, v@address, out@address)
  
  return(out)
}

#' @export
bEigen <- function(X, eigtrunc){
  #rcpp_eigen.cpp
  vals <- big.matrix(nrow = 1,
                     ncol = ncol(X),
                     init = 0,
                     type = 'double')
  vecs <- big.matrix(nrow = nrow(X),
                     ncol = ncol(X),
                     init=0,
                     type='double')
  if(is.null(eigtrunc)){
    eigtrunc <- ncol(X)
  }
  eigtrunc <- to.big.matrix(eigtrunc)
  BigEigen(X@address, eigtrunc@address, vals@address, vecs@address)
  return(list('values' = vals[,], 'vectors' = vecs*-1))
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
bTempKernel <- function(X_new, X_old, sigma){
  #rcpp_temp_kernel.cpp
  out <- big.matrix(nrow=nrow(X_new), ncol=nrow(X_old), init=0)
  s <- big.matrix(1,1,type='double', init=sigma)
  
  BigTempKernel(X_new@address, X_old@address, out@address, s@address)
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
  X <- to.big.matrix(X)
  
  if(is.null(Y)){
    Y <- deepcopy(X)
  }
  
  out <- big.matrix(nrow=nrow(X), ncol=ncol(X), init=0)
  
  BigElementwise(X@address, Y@address, out@address)
  
  return(out)
}

#' @export
bDerivatives <- function(X,sigma,K,coeffs,vcovmatc, X.sd){  
  sigma <- to.big.matrix(sigma)
  coeffs <- to.big.matrix(coeffs, d=1)
  derivatives <- big.matrix(nrow=nrow(X), ncol=ncol(X), init=-1)
  varavgderiv <- big.matrix(nrow=1, ncol=ncol(X), init=-1)
  X.sd <- to.big.matrix(X.sd, d=1)
  
  BigDerivMat(X@address, sigma@address, K@address, coeffs@address, vcovmatc@address, X.sd@address,
              derivatives@address, varavgderiv@address)
  
  return(list('derivatives'=derivatives, 'varavgderiv'=varavgderiv[]))
}
