#' bigKRLS
#' 
#' Runtime and Memory Optimized Kernel Regularized Least Squares
#' Pete Mohanty (Stanford University) and Robert Shaffer (University of Texas at Austin)
#' 
#' Kernel Regularized Least Squares (KRLS) is a kernel-based, complexity-penalized method developed by Hainmueller and Hazlett (2014) to minimize parametric assumptions while maintaining interpretive clarity. Here, we introduce bigKRLS, an updated version of the original KRLS R package with algorithmic and implementation improvements designed to optimize speed and memory usage. These improvements allow users to straightforwardly fit KRLS models to medium and large datasets (N > ~2500). 
#'
#' Major Updates
#'
#' 1. C++ integration. We re-implement most major computations in the model in C++ via Rcpp and RcppArmadillo. These changes produce up to a 50\% runtime decrease compared to the original R implementation even on a single core.
#'
#' 2. Leaner algorithm. Because of the Tikhonov regularization and parameter tuning strategies used in KRLS, the method of estimation is inherently memory-heavy O(N^2), making memory savings important even in small- and medium-sized applications. We develop and implement a new local derivatives algorithm, which reduces peak memory usage by approximately an order of magnitude, and cut the number of computations needed to find regularization parameter in half.
#'
#' 3. Improved memory management. Most data objects in R perform poorly in memory-intensive applications. We use a series of packages in the bigmemory environment to ease this constraint, allowing our implementation to handle larger datasets more smoothly.
#'
#' 4. Parallel Processing. Parallel processing with parallel makes the algorithm much faster for the marginal effects.
#'
#' 5. Interactive data visualization. We've designed an R Shiny app that allows users bigKRLS users to easily share results with collaborators or more general audiences. Simply call shiny.bigKRLS() on the outputted regression object. 
#'
#' 6. Honest p values. bigKRLS now computes p values that reflect both the regularization process and the number of predictors. For details and other options, see help(summary.bigKRLS).
#' 
#' 7. Cross-validation, including K folds crossvalidation. crossvalidate.bigKRLS performs CV, stores a number of in and out of sample statistics, as well as metadata documenting how the were split, the bigmemory file structure (if appropriate), and so on. See vignette("bigKRLS_basics") or help("crossvalidate.bigKRLS") for syntax.
#'
#' Requirements. bigKRLS is under active development. bigKRLS, as well its dependencies, require current versions of R and its compilers (and RStudio if used). For details, see \url{https://github.com/rdrr1990/code/blob/master/bigKRLS_installation.md}.
#'
#' For details on syntax, load the library and then open our vignette vignette("bigKRLS_basics"). Because of the quadratic memory requirement, users working on a typical laptop (8-16 gigabytes of RAM) may wish to start at N = 2,500 or 5,000, particularly if the number of *x* variables is large. When you have a sense of how bigKRLS runs on your system, you may wish to only estimate a subset of the marginal effects at N = 10-15,000 by setting bigKRLS(... which.derivatives = c(1, 3, 5)) for the marginal effects of the first, third, and fifth x variable. 
#' 
#' Mohanty, Pete and Robert B. Shaffer. 2016. "Messy Data, Robust Inference? Navigating Obstacles to Inference with bigKRLS" Project Presented to the International Methods Colloquium and useR! 2016. Visit \url{https://sites.google.com/site/petemohanty} for most recent version.
#' 
#' Hainmueller, Jens and Chad Hazlett. 2014. "Kernel Regularized Least Squares: Reducing Misspecification Bias with a Flexible and Interpretable Machine Learning Approach." Political Analysis. 22:143-68. \url{https://web.stanford.edu/~jhain/Paper/PA2014a.pdf} (Accessed May 20th, 2016).
#' 
#' Recent papers, presentations, and other code available at \url{github.com/rdrr1990/code/}
#' 
#' License 
#' Code released under GPL (>= 2).
#' @useDynLib bigKRLS, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats pt quantile cor sd var
#' @importFrom parallel detectCores
#' @importFrom grDevices palette
#' @importFrom ggplot2 aes element_blank geom_hline geom_point geom_smooth ggplot labs theme theme_minimal xlab ylab
#' @import bigalgebra biganalytics bigmemory shiny parallel
#' @docType package
#' @name bigKRLS
"_PACKAGE"

#' bigKRLS
#' @param y A vector of numeric observations on the dependent variable; missing values not allowed. May be base R matrix or library(bigmemory) big.matrix.
#' @param X A matrix of numeric observations of the independent variables; factors, missing values, and constant vectors not allowed. May be base R matrix or library(bigmemory) big.matrix.
#' @param sigma Bandwidth parameter, shorthand for sigma squared. Default: sigma <- ncol(X). Since x variables are standardized, facilitates interprepation of the Gaussian kernel, exp(-dist(X)^2/sigma) a.k.a the similarity score. Of course, if dist between observation i and j is 0, there similarity is 1 since exp(0) = 1. Suppose i and j differ by one standard deviation on each dimension. Then the similarity is exp(-ncol(X)/sigma) = exp(-1) = 0.368.  
#' @param derivative Logical: Estimate derivatives (as opposed to just coefficients)? Recommended for interpretability.
#' @param which.derivatives Optional. For which columns of X should marginal effects be estimated ("variables of interest"). If derivative=TRUE and which.derivative=NULL, all will marginal effects estimated (default settings). Example: out = bigKRLS(..., which.derivatives = c(1, 3, 5))
#' @param vcov.est Logical: Estimate variance covariance matrix? Required to obtain derivatives and standard errors on predictions. Default is TRUE.
#' @param Neig Number of eigenvectors and eigenvalues to calculate. The default is to calculate all N and only use those where eigval >= 0.01 max(eigval) (see eigtrunc). See out$eigenvalues and out$lastkeeper. When out$lastkeeper is much smaller than N, set Neig to be approximately out$lastkeeper for similar models, data, etc. to decrease runtime.
#' @param eigtrunc Eigentruncation, default 0.001. eigtrunc = 0.25 keeps only those eigenvectors/values such that the eigenvalue is at least 25\% of the max. If eigtrunc == 0 (default), all Neig are used to select lambda and to estimate variances.
#' @param lambda Regularization parameter. Default: estimated based (in part) on the eigenvalues of the kernel via Golden Search with convergence parameter "tolerance." Must be positive, real number. 
#' @param L Lower bound of Golden Search for lambda. 
#' @param U Upper bound of Golden Search for lambda.
#' @param tol tolerance parameter for Golden Search for lambda. Default: N / 1000.
#' @param model_subfolder_name If not null, will save estimates to this subfolder of your current working directory. Alternatively, use save.bigKRLS() on the outputted object.
#' @param overwrite.existing Logical: overwrite contents in folder 'model_subfolder_name'? If FALSE, appends lowest possible number to model_subfolder_name name (e.g., ../myresults3/). 
#' @param Ncores Number of processor cores to use. Default = ncol(X) or N - 2 (whichever is smaller). More than N - 2 NOT recommended. Uses library(parallel) unless Ncores = 1.
#' @param acf Logical. Experimental; default == FALSE. Calculate Neffective as function of mean absolute auto-correlation in X to correct p-values? Requires ncol(X) > 2. Intended for data that may violate i.i.d. To correct P values with this effective sample size, call summary(out, pval_type = "acf").
#' @param noisy Logical: Display detailed version of progress to console (intermediate output, time stamps, etc.) as opposed to minimal display? Default: if(N > 2000) TRUE else FALSE. SSH users should use X11 forwarding to see Rcpp progress display.  
#' @param instructions Display syntax after estimation with other library(bigKRLS) functions that can be used on output? Logical. (This parameter is different from noisy for the sake of crossvalidation.bigKRLS().)
#' @return bigKRLS Object containing slope and uncertainty estimates; summary() and predict() defined for class bigKRLS, as is shiny.bigKRLS().
#' @examples
#'# weight of chickens toy dataset
#'# y <- as.matrix(ChickWeight$weight) 
#'# X <- matrix(cbind(ChickWeight$Time, ChickWeight$Diet == 1), ncol = 2)
#'# out <- bigKRLS(y, X)
#'# out$R2                                     # 0.7635361
#'# summary(out, labs = c("Time", "Diet")) 
#'# save.bigKRLS(out, "exciting_results") 
#'# don't use save() unless out$has.big.matrices == FALSE
#'# out2 <- bigKRLS(y, X, which.derivatives = 2) 
#'# if x2 is variable of interest 
#' @export
bigKRLS <- function (y = NULL, X = NULL, sigma = NULL, 
                     derivative = TRUE, which.derivatives = NULL, vcov.est = TRUE,
                     Neig = NULL, eigtrunc = 0.001,
                     lambda = NULL, L = NULL, U = NULL, tol = NULL,
                     model_subfolder_name = NULL, 
                     overwrite.existing = FALSE, Ncores = NULL, 
                     acf = FALSE, noisy = NULL, instructions = TRUE)
{
  
  # Ensure RStudio is new enough for dependencies, see init.R
  check_platform()
  # Note Windows requires newer RStudio ( >= 1.1.129) than Mac or Linux ( >= 1.0.136) 
  # due to BH compatility issues. 

  if(!is.null(model_subfolder_name)){
    stopifnot(is.character(model_subfolder_name))
    
    if(!overwrite.existing & (model_subfolder_name %in% dir())){
      
      i <- 1
      tmp.name <- paste(model_subfolder_name, i, sep="")
      while(tmp.name %in% dir()){
        tmp.name <- paste(model_subfolder_name, i, sep="")
        i <- i + 1
      }
      if(model_subfolder_name %in% dir()){
        warning(cat("\na subfolder named", model_subfolder_name, 
                    "exists in your current working directory.\nYour output will be saved to", tmp.name, 
                    "instead.\nTo disable this safeguard, set bigKRLS(..., overwrite.existing=TRUE) next time.\n"))
      }
      model_subfolder_name <- tmp.name
    }
    
    dir.create(model_subfolder_name, showWarnings=FALSE)
    cat("\nmodel estimates will be saved to:\n\n", model_subfolder_name, "\n\n")
    
  }
  
  # create a folder for file backings (metadata for each big matrix etc)
  # must call to.big.matrix() with
  # to.big.matrix(... path = big.meta)
  
  big.meta <- create.metadata.dir()
  
  # suppress warnings from bigmatrix
  oldw <- getOption("warn")
  options(warn = -1)
  options(bigmemory.allow.dimnames = TRUE)
  
  stopifnot(is.matrix(X) | is.big.matrix(X))
  
  w <- list()                                             # w will become bigKRLS object
  return.big.rectangles <- is.big.matrix(X)               # X matrix, derivatives -- how to return?
  return.big.squares <- is.big.matrix(X) | nrow(X) > 2500 # Kernel, variance matrices -- how to return?
  w[["has.big.matrices"]] <- return.big.squares | return.big.rectangles
  
  noisy <- if(is.null(noisy)) nrow(X) > 2000 else noisy
  stopifnot(is.logical(noisy))
  
  if(noisy){
    
    msg <- if(w[["has.big.matrices"]]) 
      "the output will contain big.matrix objects. To avoid crashing R, use save.bigKRLS() on the output, not save().\n\n" 
    else 
      "the output will consist entirely of base R objects.\n\n"
    cat("Based on sample size (whether or not N > 2,500) and input type (base R vs. 'big' matrices),", msg)
    
  }
  
  # all X columns must have labels to prevent various post-estimation nuissance errors
  colnames(X) <- if(is.null(colnames(X))) paste0("x", 1:ncol(X)) else colnames(X)
  for(i in 1:ncol(X)){
    if(nchar(colnames(X)[i]) == 0){
      colnames(X)[i] <- paste0("x", i)
    }
  }
  xlabs <- colnames(X)
  
  w[["X"]] <- if(return.big.rectangles) to.big.matrix(X, deepcopy = TRUE, path = big.meta) else X 
  # deepcopy(X) prevents pointer to X from being inadvertently standardized 
  # in AND outside of bigKRLS()
  X <- to.big.matrix(X, path = big.meta, name = "X")
  X.init.sd <- colsd(X)
  y <- to.big.matrix(y, p = 1, path = big.meta, name = "y")
  y.init <- deepcopy(y)
  
  miss.ind <- colna(X)
  if (sum(miss.ind) > 0) { 
    stop(paste("the following columns in X contain missing data, which must be removed:", 
               paste((1:length(miss.ind))[miss.ind > 0], collapse = ', '), collapse=''))
  }
  n <- nrow(X)
  p <- ncol(X)
  # correcting p values as f(pairwise correlation of rows of X) 
  # only possible + nontrivial when ncol(X) > 2 
  acf <- acf & p > 2
  
  Neig <- if(is.numeric(Neig)) min(n, as.integer(Neig)) else n
  if(!is.numeric(eigtrunc) | eigtrunc < 0 | eigtrunc > 1)
    stop("eigtrunc must be between 0 (no truncation) and 1 (keep largest only).")
  
  if(!is.null(which.derivatives)){
    
    if(!derivative)
      stop("which.derivative requires derivative = TRUE\n\nDerivative is a logical indicating whether derivatives should be estimated (as opposed to just coefficients); which.derivatives is a vector indicating which one (with NULL meaning all).")
    stopifnot(sum(which.derivatives %in% 1:p) == length(which.derivatives))
    if(noisy){
      cat("\nMarginal effects will be calculated for the following x variables:\n")
      cat(which.derivatives, sep=", ")
    }
  }
  
  if(min(X.init.sd) == 0) 
    stop("The following columns in X are constant and must be removed: ", which(X.init.sd == 0))
  if(n != nrow(y)) 
    stop("nrow(X) not equal to number of elements in y.")
  if(colna(y) > 0) 
    stop("y contains missing data.") 
  if(colsd(y) == 0) 
    stop("y is a constant.") 
  if(!is.null(lambda))
    stopifnot(is.vector(lambda), length(lambda) == 1, is.numeric(lambda), lambda > 0)
  if(!is.null(sigma))
    stopifnot(is.vector(sigma), length(sigma) == 1, is.numeric(sigma), sigma > 0)
  
  sigma <- if(is.null(sigma)) p else sigma
  
  if (is.null(tol)) { # tolerance parameter for lambda search
    tol <- n/1000
  } else {
    stopifnot(is.vector(tol), length(tol) == 1, is.numeric(tol), tol > 0)
  }
  
  stopifnot(is.logical(derivative), is.logical(vcov.est))
  if (derivative & !vcov.est) 
    stop("vcov.est is needed to get derivatives (derivative==TRUE requires vcov.est=TRUE).")
  
  x.is.binary <- apply(X, 2, function(x){length(unique(x))}) == 2 
  if(noisy & sum(x.is.binary) > 0){
    cat(paste("\nFirst differences will be computed for the following (binary) columns of X: ", 
              toString((1:p)[x.is.binary], sep=', '), sep=""), '\n\n')
  }
  
  y.init.sd <- colsd(y.init)
  y.init.mean <- colmean(y.init)
  
  for(i in 1:ncol(X)){
    X[,i] <- (X[,i] - mean(X[,i]))/sd(X[,i])
  }
  y[,1] <- (y[,1] - mean(y[,1]))/sd(y[,1])
  
  # by default uses the same number of cores as X variables or N available - 2, whichever is smaller
  Ncores <- if(is.null(Ncores)) min(c(parallel::detectCores() - 2, ncol(X))) else Ncores
  if(noisy) cat(Ncores, "cores will be used.\n")
  
  if(noisy) cat("\nStep 1/5: Kernel (started at ", Time(), ").", sep="")
  
  K <- bGaussKernel(X, sigma)

  if(noisy) cat("\nStep 2/5: Spectral Decomposition (started at ", Time(), ").", sep="")
  
  Eigenobject <- bEigen(K, Neig, eigtrunc)

  w[["K.eigenvalues"]] <- Eigenobject$values
  w[["lastkeeper"]] <- Eigenobject$lastkeeper
  
  if (is.null(lambda)) {
    if(noisy){cat("\nStep 3/5: Golden Search for regularization parameter lambda (started at ", 
                  Time(), ").", sep="")}
    lambda <- bLambdaSearch(L = L, U = U, y = y, 
                            Eigenobject = Eigenobject, noisy = noisy)
  }else{
    if(noisy){cat("\nSkipping step 3/5, proceeding with user-inputted lambda.\n")}
  }
  
  w[["Neffective"]] <- n - sum(w[["K.eigenvalues"]]/(w[["K.eigenvalues"]] + lambda))
  if(noisy){cat("Effective Sample Size: ", w[["Neffective"]], '.', sep='')}
  
  if(noisy){cat("\n\nStep 4/5: Calculate coefficients & related estimates (started at ", 
                Time(), ").", sep="")}
  
  out <- bSolveForc(y = y, Eigenobject = Eigenobject, lambda = lambda)
  
  # bSolveForc obtains the vector of coefficients (weights) 
  # that assign importance to the similarity scores (found in K)
  if(noisy){cat("\n\nFitting values.")}
  yfitted <- K %*% to.big.matrix(out$coeffs, p = 1, path = big.meta)

  if (vcov.est == TRUE) {
    sigmasq <- bCrossProd(y - yfitted)[]/n
    if(noisy){cat("\nIn standardized units, sigmasq = ", round(sigmasq, 5), ".", sep='')}
    if(noisy){cat("\nCalculating variance-covariance of the coefficients.")}
    
    # subsetting from Neig < N and/or eigtrunc now handled by bEigen()
    m <- bMultDiag(Eigenobject$vectors, 
                   sigmasq * (Eigenobject$values + lambda)^-2)
    vcovmatc <- bTCrossProd(m, Eigenobject$vectors)
    
    remove(Eigenobject)
    remove(m)
    gc()
    if(noisy) cat("\nEstimating variance covariance of the fitted values.")
    vcovmatyhat <- bCrossProd(K, vcovmatc %*% K)
    if(!is.null(model_subfolder_name) & return.big.squares){
      vcovmatyhat <- (y.init.sd^2) * vcovmatyhat
      cat("\nsaving vcovmatyhat to", getwd())
      write.big.matrix(x = vcovmatyhat, 
                       filename = file.path(model_subfolder_name, "vcovmatyhat.txt"))
      remove(vcovmatyhat)
      cat("\nvcovmatyhat successfully saved to disk (and removed from memory for speed).\n")
    }
    
  }else {
    vcov.est.c <- NULL
    vcov.est.fitted <- NULL
  }
  if (derivative == TRUE) {
    
    if(noisy){cat("\n\nStep 5/5: Estimate marginal effects and their uncertainty (started at ", 
                  Time(), ").\n\n", sep="")}
    
    X_estimate <- if(!is.null(which.derivatives)) deepcopy(X, cols = which.derivatives) else X
    
    if(Ncores == 1){
        deriv_out <- bDerivatives(X_estimate, sigma, K, out$coeffs, vcovmatc)
    }else{
      
      X_index <- if(is.null(which.derivatives)) 1:p else which.derivatives

      # each core will need to know how to find the big matrices
      # writing their description to disk will allow each core to do that...
      
      K <- to.big.matrix(K, name = "K", path = big.meta)
      vcovmatc <- to.big.matrix(vcovmatc, name = "V", path = big.meta)
      
      if(!("cl" %in% ls())){
        cl <- if(noisy) makeCluster(Ncores, outfile='') else makeCluster(Ncores)
        clusterEvalQ(cl, suppressPackageStartupMessages(library(bigKRLS)))
      } 
      
      tmp <- parLapply(cl, X_index, function(i, sigma, coefficients, X.init.sd, path){

        # each core finds the big matrices like so...
        X <- attach.resource(dget(file.path(path, "X.desc")), path = path)
        K <- attach.resource(dget(file.path(path, "K.desc")), path = path)
        V <- attach.resource(dget(file.path(path, "V.desc")), path = path)
        
        # the description that dget obtains doesn't contain a path variable
        # which prevents attach.big.matrix() from working. but it's just a wrapper
        # for attach.resource()...
          
        x <- deepcopy(X, cols = i)
        
        output <- bDerivatives(x, sigma, K, coefficients, V)
        # can't return pointers
        list(output[[1]][], output[[2]])
        # output is small but could also use to.big.matrix for reverse direction (should if N * N needed)
      }, sigma, out$coeffs, X.init.sd, big.meta)
      stopCluster(cl) 
      remove(cl)
      
      derivs <- matrix(nrow = n, ncol = ncol(X_estimate))
      varavgderiv <- c()
      for(i in 1:ncol(X_estimate)){
        derivs[,i] <- tmp[[i]][[1]]
        varavgderiv[i] <- tmp[[i]][[2]]
      }
      deriv_out <- list()
      deriv_out[["derivatives"]] <- to.big.matrix(derivs, path = big.meta, name = "derivs") 
      deriv_out[["varavgderiv"]] <- varavgderiv
      remove(tmp, derivs, varavgderiv)
    }
    
    
    if(noisy){
      cat('\n\nFinshed (', Time(), ').', sep="") 
      cat('\n\nPrepping bigKRLS output object...\n')
    }
    
    derivmat <- deriv_out$derivatives
    varavgderivmat <- deriv_out$varavgderiv
    remove(deriv_out)
    
    # Pseudo R2 using only Average Marginal Effects
    yhat_ame <- (X_estimate[] %*% colMeans(derivmat[]))^2

    w[["R2AME"]] <- cor(y.init[], yhat_ame)^2
    
    derivmat <- y.init.sd * derivmat
    for(i in 1:ncol(derivmat)){
      derivmat[,i] <- derivmat[,i]/X.init.sd[i]
    }
    
    attr(derivmat, "scaled:scale") <- NULL
    avgderiv <- matrix(colmean(derivmat), nrow=1)
    attr(avgderiv, "scaled:scale") <- NULL
    
    if(is.null(which.derivatives)){
      varavgderivmat <- matrix((y.init.sd/X.init.sd)^2 * as.matrix(varavgderivmat), nrow = 1)
    }else{
      varavgderivmat <- matrix((y.init.sd/X.init.sd[which.derivatives])^2 * as.matrix(varavgderivmat), nrow=1)
    }
    
    attr(varavgderivmat, "scaled:scale") <- NULL
  }
  
  if(acf){
    if(noisy){cat('Accumulating absolute pairwise correlations within X to correct p-values; see help(bigKRLS).')}
    Neffective.acf <- bNeffective(X)
    if(noisy){cat("\nEffective Sample Size as f(absolute correlation of X): ", Neffective.acf, '.', sep='')}
  }
  
  # w is the output object
  
  w[["coeffs"]] <- out$coeffs
  w[["y"]] <- y.init[]
  w[["sigma"]] <- sigma
  w[["lambda"]] <- lambda 
  w[["binaryindicator"]] <- x.is.binary
  w[["which.derivatives"]] <- which.derivatives
  w[["xlabs"]] <- xlabs
  
  w[["yfitted"]] <- yfitted <- as.matrix(yfitted) * y.init.sd + y.init.mean
  w[["R2"]] <- 1 - (var(y.init - yfitted)/(y.init.sd^2))
  w[["Looe"]] <- out$Le * y.init.sd
  w[["Neffective.acf"]] <- if(exists("Neffective.acf")) Neffective.acf else NULL
  
  # returning base R matrices when sensible...
  w[["K"]] <- if(return.big.squares) K else K[] 
  
  if (vcov.est) {
    
    vcovmatc <- (y.init.sd^2) * vcovmatc
    
    if(return.big.squares){
      
      w[["vcov.est.c"]] <- vcovmatc
      
      if(is.null(model_subfolder_name)){
        vcovmatyhat <- (y.init.sd^2) * vcovmatyhat
        w[["vcov.est.fitted"]] <- vcovmatyhat
      } # vcovmatyhat already saved otherwise
      
    }else{
      w[["vcov.est.c"]] <- vcovmatc[]
      w[["vcov.est.fitted"]] <- vcovmatyhat[]
    }
  }
  
  w[["derivative.call"]] <- derivative
  
if(derivative){
    rownames(avgderiv) <- rownames(varavgderivmat) <- ""
    
    w[["avgderivatives"]] <- avgderiv
    w[["var.avgderivatives"]] <- varavgderivmat
    w[["derivatives"]] <- if(return.big.rectangles) derivmat else derivmat[]
    
    if(p == 1 & !return.big.rectangles){
      w$derivatives <- matrix(w$derivatives)
      w$avgderivatives <- matrix(w$avgderivatives)
    }
    colnames(w$derivatives) <- colnames(w$avgderivatives) <- if(is.null(which.derivatives)) xlabs else xlabs[which.derivatives]
  }
  
  if(!is.null(model_subfolder_name)){
    
    cat("\nsaving output to", model_subfolder_name, "\n")
    w[["path"]] <- normalizePath(model_subfolder_name)
    
    for(i in which(unlist(lapply(w, is.big.matrix)))){
      output_file = file.path(model_subfolder_name, paste0(names(w)[i], ".txt"))
      cat("\twriting", output_file, "...\n")
      write.big.matrix(x = w[[i]], filename = output_file) 
                       # col.names = !is.null(colnames(w[[i]])) 
                       # deprecating, handling with object$xlabs
    }
    
    Nbm <- sum(unlist(lapply(w, is.big.matrix))) + return.big.squares
    cat("\n\n", Nbm, "matrices saved as big matrices.\n") 
    if(Nbm == 0){
      cat(" (base R save() may be used safely in this case too).\n")
    }else{
      cat("\nto reload, use syntax like:\n\nload.bigKRLS(\"", w$path, "\")\n or\n",
          "load.bigKRLS(\"", w$path, "\", newname=\"my_estimates\")\n", sep="")}
    if(Nbm > 0){
      bigKRLS_out <- w[-which(unlist(lapply(w, is.big.matrix)))]
      class(bigKRLS_out) <- "bigKRLS"
    }else{
      bigKRLS_out <- w
    }
    stopifnot(sum(unlist(lapply(bigKRLS_out, is.big.matrix))) == 0)
    save(bigKRLS_out, file=file.path(model_subfolder_name, "estimates.RData"))
    cat("\nbase R elements of the output saved to estimates.RData.\n")
    cat("Total file size approximately", 
        round(sum(file.info(list.files(path = model_subfolder_name, full.names = TRUE))$size)/1024^2), 
        "megabytes.\n\n")
    model_subfolder_name
  }
  
  unlink(big.meta, recursive = TRUE)
  # file.remove(dir(path = big.meta, full.names = TRUE))
  # description are pointers that will crash R outside of current R session; 
  # removing their footprint
  
  if(instructions) cat("\nAll done. You may wish to use summary() for more detail, predict() for out-of-sample forecasts, or shiny.bigKRLS() to interact with results. For an alternative approach, see help(crossvalidate.bigKRLS). Type vignette(\"bigKRLS_basics\") for sample syntax. Use save.bigKRLS() to store results and load.bigKRLS() to re-open them.\n\n")
  class(w) <- "bigKRLS" 
  return(w)
  
  options(warn = oldw)
}  

###################
# bigKRLS exports # 
###################

#' predict.bigKRLS
#' 
#' Predict function for bigKRLS object. crossvalidate.bigKRLS() provides additional functionality.
#' 
#' @param object bigKRLS output
#' @param newdata new data. ncol(X) == ncol(newdata) but nrow(X) need not be the same as nrow(newdata).
#' @param se.pred get standard errors on predictions?
#' @param ytest Provide testing data to have it returned with the object. Optional. To automatically generate out-of-sample test statistics, use crossvalidate.bigKRLS() instead.
#' @param ... ignore
#' @return Returns bigKRLS_predicted list object.
#' @examples  
#'# y <- as.matrix(ChickWeight$weight)
#'# X <- matrix(cbind(ChickWeight$Time, ChickWeight$Diet == 1), ncol = 2)
#'# N <- length(y)
#'# set.seed(123)
#'# train <- sample(N, 100, replace = FALSE)
#'# outOfSample <- c(1:N)[-train]
#'# test <- sample(outOfSample, 10, replace = FALSE)
#'# fit <- bigKRLS(y[train], X[train,], instructions = FALSE) 
#'# p <- predict(fit, X[test,])
#'# range(p$predicted) # 44.04614 257.76520
#' @method predict bigKRLS
#' @export
predict.bigKRLS <- function (object, newdata, se.pred = FALSE, ytest = NULL, ...) 
{
  if (class(object) != "bigKRLS") {
    warning("Object not of class 'bigKRLS'")
    UseMethod("predict")
    return(invisible(NULL))
  }
  if(se.pred == TRUE) {
    if (is.null(object$vcov.est.c)) {
      stop("recompute bigKRLS object with bigKRLS(,vcov.est=TRUE) to compute standard errors")
    }
  }
  
  # convert everything to a bigmatrix for internal usage
  
  big.meta <- create.metadata.dir()
  
  object$X <- to.big.matrix(object$X, deepcopy = TRUE, path = big.meta)
  object$K <- to.big.matrix(object$K, path = big.meta) #, description = "K")
  object$vcov.est.c <- to.big.matrix(object$vcov.est.c, path = big.meta)

  if(!is.null(object$vcov.est.fitted)){
    object$vcov.est.fitted <- to.big.matrix(object$vcov.est.fitted, path = big.meta)  
  }else{
    cat("vcov.est.fitted not found in bigKRLS object, attempting to load from object's path,\n ", object$path)
    if("vcov.est.fitted.txt" %in% dir(object$path)){
      object$vcov.est.fitted <- read.big.matrix(filename = file.path(object$path, "vcov.est.fitted.txt"),
                                                type = 'double')
      cat("\nVariance(y hat) matrix loaded successfully\n")
    }else{
      cat("\nFile not found.\n")
    }
  }
  
  # flag: return big matrices? (new kernel, etc...)
  bigmatrix.in <- is.big.matrix(newdata) | object$has.big.matrices
  
  newdata.init <- newdata
  newdata <- to.big.matrix(newdata, deepcopy = TRUE, path = big.meta)
  
  if (ncol(object$X) != ncol(newdata)) {
    stop("ncol(newdata) differs from ncol(X) from fitted bigKRLS object")
  }
  Xmeans <- colmean(object$X, na.rm = TRUE)
  Xsd <- colsd(object$X, na.rm = TRUE)
  
  for(i in 1:ncol(object$X)){
    object$X[,i] <- (object$X[,i] - Xmeans[i])/Xsd[i]
  }  
  
  for(i in 1:ncol(newdata)){
    newdata[,i] <- (newdata[,i] - Xmeans[i])/Xsd[i]
  }
  
  newdataK <- bTempKernel(newdata, object$X, object$sigma)
  
  ypred <- (newdataK %*% to.big.matrix(object$coeffs, path = big.meta))[]
  
  if (se.pred) {
    
    # vcov.est.c.raw <- object$vcov.est.c * (1/var(object$y))
    # vcov.est.pred <- var(object$y) * bTCrossProd(newdataK %*% vcov.est.c.raw, newdataK)
    # remove(vcov.est.c.raw)
    vcov.est.pred <- var(object$y) * bTCrossProd(newdataK %*% (object$vcov.est.c * (1/var(object$y))), newdataK)
    se.pred <- sqrt(bDiag(vcov.est.pred))
    # se.pred <- matrix(sqrt(diag(vcov.est.pred[])), ncol = 1) 
    
  }
  else {
    vcov.est.pred <- se.pred <- NULL
  }
  
  ypred <- ypred * sd(object$y) + mean(object$y)
  
  if(!bigmatrix.in){
    vcov.est.pred <- vcov.est.pred[]
    newdataK <- newdataK[]
  }
  
  out <- list(predicted = ypred, se.pred = se.pred, vcov.est.pred = vcov.est.pred, 
              newdata = newdata.init, newdataK = newdataK, 
              has.big.matrices = bigmatrix.in, # TRUE if bigKRLS returned big OR user inputted to predict()
              ytest = ytest)
  
  class(out) <- "bigKRLS_predicted"
  unlink(big.meta, recursive = TRUE)
  return(out)
  
}

#' summary.bigKRLS
#' 
#' Summary function for bigKRLS output. 
#' 
#' @param object bigKRLS output. If you saved with save.bigKRLS(), only the .RData file is needed for this function.
#' @param degrees "Neffective" (default) or "N". What value should be used as the sample size for the t-tests of the the AMEs (average marginal effects)? If 'Neffective' (default), degrees of freedom for t tests reflects degrees of freedom used to obtain regularization parameter, lambda. Neffective = N - sum(eigenvalues/(eigenvalues + lambda)); see e.g. Hastie et al. (2015, 61-68). 'N' is simply the observed sample size (note this is the default for library(KRLS)). Degrees of freedom for t-tests is either Neffective - P or N - P.
#' @param probs For quantiles of the marginal effects of each x variable.
#' @param digits Number of signficant digits.
#' @param labs Optional vector of x labels.
#' @param ... ignore
#' @return Returns list with "ttests" (Average Marginal Effect estimates, standard errors, t-values, and p values) and "percentiles" (of the marginal effects).
#' @method summary bigKRLS
#' @examples
#'# y <- as.matrix(ChickWeight$weight)
#'# X <- matrix(cbind(ChickWeight$Time, ChickWeight$Diet == 1), ncol = 2)
#'# out <- bigKRLS(y, X, Ncores=1)
#'# summary(out)
#'# s = summary(out, digits = 2, labs = c("Time", "ChickWeightDiet"))
#'# knitr::kable(s[["ttests"]])     # to format with RNotebook or RMarkdown
#'# knitr::kable(s[["percentiles"]])
#'# summary(out, degrees = "N")     # don't adjust p value for cost of lambda. See above.
#'@export
summary.bigKRLS <- function (object, degrees = "Neffective", probs = c(0.05, 0.25, 0.5, 0.75, 0.95), 
                             digits = 4, labs = NULL, ...) 
{
  if (class(object) != "bigKRLS") {
    warning("Object not of class 'bigKRLS'")
    UseMethod("summary")
    return(invisible(NULL))
  }
  
  N <- n <- nrow(object$X)
  
  stopifnot(degrees %in% c("acf", "Neffective", "N"))
  
  if(degrees == "Neffective"){
    n <- object$Neffective
  }
  if(degrees == "acf"){
    if(is.null(object$Neffective.acf)){
      big.meta <- create.metadata.dir()
      n <- bNeffective(to.big.matrix(scale(object$X[]), path = big.meta))
      unlink(big.meta)
      cat("\n\n\n")
    }else{
      n <- object$Neffective.acf
    } 
  } 
  
  cat("\n\nMODEL SUMMARY:\n\n")
  cat("lambda:", round(object$lambda, digits), "\n")
  cat("N:", N, "\n")
  if(n != N) cat("N Effective:", n, "\n")
  
  p <- ncol(object$X)
  cat("R2:", round(object$R2, digits), "\n")
  
  if (is.null(object$derivatives)) {
    cat("\nrecompute with bigKRLS(..., derivative = TRUE) for estimates of marginal effects\n")
    return(invisible(NULL))
  }
  
  if(!is.null(object$R2AME))
    cat("R2AME**:", round(object$R2AME, digits), "\n\n")
  
  if(!is.null(labs)){
    stopifnot(length(labs) == p)
    colnames(object$X) <- labs
  }else{
    if(is.big.matrix(object$X)) options(bigmemory.allow.dimnames=TRUE)
    colnames(object$X) <- object$xlabs
  }
  
  if(is.null(object$which.derivatives)){
    object$which.derivatives <- 1:p
  }
  
  est <- object$avgderivatives
  se <- sqrt(object$var.avgderivatives)
  if(degrees != "Neffective"){
    se <- se*nrow(object$X)/n # correcting variance estimate
  }
  tval <- est/se
  pval <- 2 * pt(abs(tval), n - p, lower.tail = FALSE)
  AME <- t(rbind(est, se, tval, pval))
  colnames(AME) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(AME) <- colnames(object$X)[object$which.derivatives]
  if (sum(object$binaryindicator[object$which.derivatives]) > 0) {
    tmp <- rownames(AME)[object$binaryindicator[object$which.derivatives]]
    rownames(AME)[object$binaryindicator[object$which.derivatives]] <- paste(tmp, "*", sep="")
  }
  cat("Average Marginal Effects:\n\n")
  print(round(AME, digits))
  
  cat("\n\nPercentiles of Marginal Effects:\n\n")
  
  deriv <- matrix(object$derivatives[], ncol = length(object$which.derivatives))
  qderiv <- t(apply(deriv, 2, quantile, probs = probs, na.rm = TRUE))
  rownames(qderiv) <- rownames(AME)
  print(round(qderiv, digits))
  
  if (sum(object$binaryindicator) > 0) {
    cat("\n(*) Reported average and percentiles of dy/dx is for discrete change of the dummy variable from min to max (usually 0 to 1)).\n\n")
  }
  cat("\n(**) Pseudo-R^2 computed using only the Average Marginal Effects.") 
  if(length(object$which.derivatives) != ncol(object$X)) cat(" NOTE: If only a subset of marginal effects were estimated, Pseudo-R^2 calculated with that subset.")
  cat("\n\n")
  cat("\nYou may also wish to use predict() for out-of-sample forecasts or shiny.bigKRLS() to interact with results. Type vignette(\"bigKRLS_basics\") for sample syntax. Use save.bigKRLS() to store results and load.bigKRLS() to re-open them.\n\n")
  ans <- list(ttests = AME, 
              percentiles = qderiv)
  class(ans) <- "summary.bigKRLS"
  return(invisible(ans))
  
}

#' summary.bigKRLS_CV
#' 
#' Summary function for bigKRLS crossvalidated output.
#' 
#' @param object bigKRLS_CV output. If you saved with save.bigKRLS(), only the .RData file is needed for this function (for K folds CV, that means only the .RData in the top level folder).
#' @param ... Additional parameters to be passed to summary() for the training model(s). For example, summary(cv, digits = 3). See ?bigKRLS.summary for details.
#' @method summary bigKRLS_CV
#' @examples
#'# y <- as.matrix(ChickWeight$weight)
#'# X <- matrix(cbind(ChickWeight$Time, ChickWeight$Diet == 1), ncol = 2)
#'# cv.out <- crossvalidate.bigKRLS(y, X, seed = 123, ptesting = 20)
#'# summary(cv.out)
#'# cv <- summary(cv.out, labs = c("Alpha", "Beta", "Gamma", "Delta", "Epsilon"))
#'# cv$training.ttests            # hypothesis tests, training model
#'# kcv.out <- crossvalidate.bigKRLS(y, X, seed = 123, Kfolds = 3)
#'# summary(kcv.out)           
#'# kcv <- summary(kcv.out) 
#'# kcv$overview                 # test stats, in + out of sample, all folds
#'# kcv$training2.ttests         # hypothesis tests, fold 2 
#' @export
summary.bigKRLS_CV <- function (object, ...) 
{
  if (class(object) != "bigKRLS_CV") {
    warning("Object not of class 'bigKRLS_CV'")
    UseMethod("summary")
    return(invisible(NULL))
  }
  
  arguments = list(...)
  digits <- if("digits" %in% names(arguments)) digits else 3
  
  if(object$type == "crossvalidated"){
    
    cat("Overview of Model Performance\n\n")
    
    cat("N:", length(unlist(object$indices)), "\n")
    cat("Seed:", object$seed, "\n\n")
    
    overview <- matrix(ncol = 2, nrow = 6)
    colnames(overview) <- c("In Sample", "Out of Sample")
    rownames(overview) <- c("Mean Squared Error (Full Model)", 
                            "Mean Squared Error (Average Marginal Effects Only)",
                            "Pseudo-R^2 (Full Model)",
                            "Pseudo-R^2 (Average Marginal Effects Only)",
                            "",
                            "N")
    
    overview[1, 1] <- object$MSE_is
    overview[1, 2] <- object$MSE_oos
    
    overview[2, 1] <- object$MSE_AME_is
    overview[2, 2] <- object$MSE_AME_oos
    
    overview[3, 1] <- object$pseudoR2_is
    overview[3, 2] <- object$pseudoR2_oos
    
    overview[4, 1] <- object$pseudoR2AME_is
    overview[4, 2] <- object$pseudoR2AME_oos
    
    overview[6, 1] <- length(object$indices$train.set)
    overview[6, 2] <- length(object$indices$test.set)
    
    print(round(overview, digits = digits), na.print = "")
    
    cat("\n\nSummary of Training Model:\n")
    z <- summary(object$trained, ...)
    
    ans <- list(overview = overview,
                training.ttests = z$ttests, 
                training.percentiles = z$percentiles)
    class(ans) <- "summary.bigKRLS_CV"
    
  }else{
    
    stopifnot(object$type == "KfoldsCV")
    # object[lapply(object, class) == "numeric"]
    
    cat("Overview of Model Performance\n\n")
    
    cat("N:", length(unlist(object$folds)), "\n")
    cat("Kfolds:", object$Kfolds, "\n")
    cat("Seed:", object$seed, "\n\n")
    
    overview <- matrix(unlist(object[lapply(object, class) == "numeric"][-c(1:2)]), 
                       ncol=object$Kfolds, byrow = TRUE)
    colnames(overview) <- paste("Fold", 1:object$Kfolds)
    
    # somewhat cumbersome but the test stats will differ depending on whether 
    # user computes with bigKRLS(... derivative = TRUE)
    labs <- unlist(strsplit(names(unlist(object[lapply(object, class) == "numeric"][-c(1:2)])), "\\."))
    labs <- unique(labs[-grep("fold", labs)])
    labs <- gsub("_is", " (In Sample)", labs)
    labs <- gsub("_oos", " (Out of Sample)", labs)
    labs <- gsub("_", " ", labs)
    labs <- gsub("R2AME", "R2 AME", labs)
    
    rownames(overview) <- labs
    overview <- overview[match(sort(labs), labs), ]
    
    print(overview, digits = digits)
    class(overview) <- "summary.bigKRLS_CV"
    ans <- list(overview = overview)
    
    cat("\nMSE denotes Mean Squared Error. AME implies calculations done with Average Marginal Effects only.")
    
    for(k in 1:object$Kfolds){
      cat("\n\nSummary of Training Model", k , ":\n", sep="")
      z <- summary(object[[paste0("fold_", k)]][["trained"]], ...)
      ans[[paste0("training", k, ".ttests")]] <- z$ttests
      ans[[paste0("training", k, ".percentiles")]] <- z$percentiles
    }
    
  }
  
  return(invisible(ans))
  
}

#' save.bigKRLS
#' 
#' save function, recommended when bigKRLS output contains big matrices (once N > 2,500 the kernel is stored this way).
#' Base R data will be stored in a list in an .RData file, big matrices will be stored in .txt files. 
#' Call load.bigKRLS() to retrieve. 
#' 
#' @param object bigKRLS output (regression, prediction, and crossvalidation). Use load.bigKRLS(model_subfolder_name), not load().
#' @param model_subfolder_name A name of a folder where the file(s) will be written. 
#' @param overwrite.existing Logical -- write over folders with the same name? Default == FALSE.
#' @param noisy Logical -- display progress, additional instructions? Default == TRUE.
#' @examples
#'# y <- as.matrix(ChickWeight$weight)
#'# X <- matrix(cbind(ChickWeight$Time, ChickWeight$Diet == 1), ncol = 2)
#'# out <- bigKRLS(y, X, Ncores=1)
#'# save.bigKRLS(out, "exciting_results") 
#'# don't use save() unless out$has.big.matrices == FALSE
#'# load.bigKRLS("/path/to/exciting_results") 
#'# path not necessary if in current working directory
#' @export
save.bigKRLS <- function (object, model_subfolder_name, 
                          overwrite.existing = FALSE, noisy = TRUE) 
{
  
  bClasses <- c("bigKRLS", "bigKRLS_predicted", "bigKRLS_CV")
  if (!(class(object) %in% bClasses)) {
    warning("Object not a bigKRLS class.")
    UseMethod("save")
    return(invisible(NULL))
  }
  stopifnot(is.character(model_subfolder_name))
  
  object <- make_path(object, model_subfolder_name, overwrite.existing)
  
  if(class(object) == "bigKRLS" | class(object) == "bigKRLS_predicted"){
    bSave(object, noisy)
  }else{
    
    for(k in grep("fold_", names(object))){
      
      object[[k]][["trained"]] <- make_path(object[[k]][["trained"]], 
                                            file.path(object[["model_subfolder_name"]], names(object)[k], "trained"), 
                                            overwrite.existing = TRUE) # overwrite TRUE here in keeping with above choice
      bSave(object[[k]][["trained"]], noisy)
      
      object[[k]][["tested"]] <- make_path(object[[k]][["tested"]], 
                                           file.path(object[["model_subfolder_name"]], names(object)[k], "tested"), 
                                           overwrite.existing = TRUE)
      bSave(object[[k]][["tested"]], noisy)
      
    }
    object[["dir"]] <- dir(object[["model_subfolder_name"]], recursive = TRUE)
    tmp <- class(object)
    object <- object[-grep("fold_", names(object))]
    class(object) <- tmp
    save(object, file = file.path(object[["model_subfolder_name"]], "estimates.RData"))
    if(noisy) cat("\nBase R elements, summary stats, and metadata of the entire cross-validated outputted object saved in:", 
                  file.path(object[["model_subfolder_name"]], "estimates.RData\n"))
  }
  
  if(noisy) cat("\nTotal file size approximately", 
                round(sum(file.info(list.files(path = object[["model_subfolder_name"]], full.names = TRUE, recursive = TRUE))$size)/1024^2), 
                "megabytes.\n")
  
}

#' load.bigKRLS
#' 
#' Reconstructs bigKRLS output object as list.
#' 
#' @param path Path to folder where bigKRLS object was saved. 
#' @param newname If NULL (default), bigKRLS regression and prediction output will appear as "bigKRLS_out", while crossvalidation results will appear as "object".
#' @param pos position. Default == 1 (global environment). NULL means don't assign (return only).
#' @param return_object Logical: return library(bigKRLS) object? Default == FALSE. 
#' @param noisy Logical: display updates?
#' @examples
#'# save.bigKRLS(out, "exciting_results") # don't use save()
#'# load.bigKRLS("exciting_results") # don't use load()
#' @export
load.bigKRLS <- function(path, newname = NULL, pos = 1, noisy = TRUE, return_object = FALSE){
  
  stopifnot(is.null(newname) | is.character(newname))
  
  files <- dir(path = path)
  if(!(tolower("estimates.RData") %in% tolower(files))){
    stop("estimates.RData not found. Check the path to the output folder.\n\nNote: for any files saved manually, note that load.bigKRLS() anticipates the convention used by save.bigKRLS: estimates.RData stores the base R objects in a list called bigKRLS_out, big matrices stored as text files named like they are in bigKRLS objects (object$K becomes K.txt, etc.).\n\n")
  }
  
  # bigKRLS version > 1.5 uses more standard convention .RData but loads .rdata too...
  # name = load(file.path(path, "estimates.RData"))
  name <- load(file.path(path, files[match(tolower("estimates.RData"), tolower(files))]))
  # name will be 'bigKRLS_out' for predict and regression objects
  # and 'object' for CV
  
  # equivalent to class(get(name)) == "bigKRLS" | class(get(name)) == "bigKRLS_predicted"
  if(name == "bigKRLS_out"){
    
    bigKRLS_out <- bLoad(bigKRLS_out, path, noisy) # load big.matrix objects, if any
    
  }else{
    
    stopifnot(name == "object" & class(get(name)) == "bigKRLS_CV")
    
    est <- file.path(path, object$dir[grep("estimates.RData", object$dir)])
    
    for(i in 1:length(est)){
      
      load(est[i])  # loads next estimates.RData, necessarily named bigKRLS_out
      kind <- if(i %% 2 == 0) "trained" else "tested"
      fold <- paste0("fold_", ceiling(i/2))
      object[[fold]][[kind]] <- bLoad(bigKRLS_out, path, noisy)
      
    }
    
  } 
  
  
  if(!is.null(pos)) {
    
    if(is.null(newname)){
      newname <- name
    }
    
    if(class(get(name)) == "bigKRLS" | class(get(name)) == "bigKRLS_predicted"){
      assign(newname, bigKRLS_out, envir = as.environment(pos))
    }else{
      assign(newname, object, envir = as.environment(pos))
    }
    
    if(noisy) cat("\nNew object created named", 
                  newname, ".\n\nRun vignette(\"bigKRLS_basics\") for examples for a", 
                  class(bigKRLS_out), "object.")
  }
  
  if(return_object) return(bigKRLS_out)
  
}


#' shiny.bigKRLS
#' 
#' @param out bigKRLS output. Does not require any N * N matrices.
#' @param export Logical -- instead of running Shiny, prepare the key estimates as a separate file? (The N * Ns are too large for most Shiny servers but the key estimates are only N * P).
#' @param main.label Main label (upper left of app)
#' @param plot.label Optional character
#' @param xlabs Optional character vector for the x variables.
#' @param font_size Font size. Default == 20. calls "ggplot2::theme_minimal(base_size = font_size)"
#' @param shiny.palette color scheme for main app. 9 colors.
#' @param hline horizontal line. Default == 0 (x axis)
#' @examples
#'# N <- 500  # proceed with caution above N = 5,000 for system with 8 gigs made avaiable to R
#'# P <- 4
#'# X <- matrix(rnorm(N*P), ncol=P)
#'# b <- 1:P 
#'# y <- sin(X[,1]) + X %*% b + rnorm(N)
#'# out <- bigKRLS(y, X, Ncores=1)
#'# shiny.bigKRLS(out, "exciting_results", "The Results", c("Frequency", "xA", "xB", "xC")) # not run
#' @export
shiny.bigKRLS <- function(out, export=FALSE, main.label = "bigKRLS estimates", 
                          plot.label = NULL, xlabs = NULL, font_size = 20, hline = 0,
                          shiny.palette = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                                            "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")){
  
  if(!export){cat("export set to FALSE; set export to TRUE to prepare files for another machine.")}
  
  if(is.null(xlabs)) xlabs = out$xlabs
  
  colnames(out$X) <- xlabs
  dydxlabs <- if(is.null(out$which.derivatives)) xlabs else xlabs[out$which.derivatives]
  colnames(out$derivatives) <- names(out$avgderivatives) <- names(out$var.avgderivatives) <- out$dydxlabs <- dydxlabs
  
  palette(shiny.palette)
  
  bigKRLS_server <- shinyServer(function(input, output, session) {
    
    selectedData <- reactive({
      
      return(list(x = as.numeric(out$X[, input$xp]),
                  derivatives = as.numeric(out$derivatives[, input$dydxp])))
    })
    
    output$graph <- renderPlot({
      
      P = ggplot(NULL) 
      P = P + geom_point(aes(x = selectedData()[["x"]], y = selectedData()[["derivatives"]]), 
                         alpha = 1, size=.1, color='grey') 
      P = P +  geom_smooth(aes(x = selectedData()[["x"]], y = selectedData()[["derivatives"]]),
                           method='loess') + xlab(input$xp) 
      P = P +  ylab(paste('Marginal Effects of ', input$dydxp)) 
      P = P +  geom_hline(aes(yintercept=hline))
      P = P +  theme_minimal(base_size = font_size)
      P = P +  theme(panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.background=element_blank()) 
      P = P + labs(title = plot.label)
      P
      
    })
  })
  
  bigKRLS_ui <- shinyUI(fluidPage(
    
    titlePanel(main.label),
    
    sidebarPanel(
      selectInput('dydxp', 'Marginal Effects (dy/dx)', colnames(out$derivatives)),
      selectInput('xp', 'x', colnames(out$X))
    ),
    
    mainPanel(plotOutput('graph'))
    
  ))
  
  if(export){
    
    output_baseR <- out
    output_baseR$K <- output_baseR$vcov.c <- output_baseR$vcov.fitted <- NULL
    for(i in which(unlist(lapply(output_baseR, is.big.matrix)))){
      output_baseR[[i]] <- as.matrix(output_baseR[[i]])
    }
    
    save(output_baseR, file="shiny_out.RData")
    
    cat("A re-formatted version of your output has been saved with file name \"shiny_out.RData\" in your current working directory:\n", getwd(),
        "\nFor a few technical reasons, the big N * N matrices have been removed and the smaller ones converted back to base R;\nthis should make your output small enough for the free version of Shiny's server.\nTo access the Shiny app later or on a different machine, simply execute this script with the following commands:\n",
        "\nload(\"shiny_out.RData\")\nNext, execute this code:\n\nshiny.bigKRLS(output_baseR)")
  }else{
    shinyApp(ui = bigKRLS_ui, server = bigKRLS_server)
  }
}


#' crossvalidate.bigKRLS
#' 
#' @param y A vector of numeric observations on the dependent variable; missing values not allowed. May be base R matrix or library(bigmemory) big.matrix.
#' @param X A matrix of numeric observations of the independent variables; factors, missing values, and constant vectors not allowed. May be base R matrix or library(bigmemory) big.matrix.
#' @param seed Seed to be used when partitioning data. For example, crossvalidate.bigKRLS(..., seed = 123). ?set.seed for details.
#' @param Kfolds Number of folds for cross validation. Requires ptesting == NULL. Note KRLS assumes variation in each column; rare events or rarely observed factor levels may violate this assumption if Kfolds is too large given the data.
#' @param ptesting Percentage of data to be used for testing (e.g., ptesting = 20 means 80\% training, 20\% testing). Requires Kfolds == NULL. Note KRLS assumes variation in each column; rare events or rarely observed factor levels may violate this assumptions if ptesting is too small given the data.
#' @param estimates_subfolder If non-null, saves all model estimates in current working directory.
#' @param ... Additional arguments to be passed to bigKRLS() or predict(). E.g., crossvalidate.bigKRLS(y, X, derivative = FALSE) will run faster but compute fewer test stats comparing in and out of sample performance (because the marginal effects will not be estimated).
#' @return bigKRLS_CV (list) Object of estimates and summary stats; summary() is defined. For train/test, contains a bigKRLS regression object and a predict object. For Kfolds, contains a nested series of training and testing models. 
#' @examples
#'# y <- as.matrix(ChickWeight$weight)
#'# X <- matrix(cbind(ChickWeight$Time, ChickWeight$Diet == 1), ncol = 2)
#'# cv.out <- crossvalidate.bigKRLS(y, X, seed = 123, ptesting = 20)
#'# cv.out$pseudoR2_oos # cor(y[test], cv.out$tested$predicted)^2 == 0.7488783
#'# cv <- summary(cv.out)
#'# cv$training.ttests                      # hypothesis tests, training model
#'# kcv.out <- crossvalidate.bigKRLS(y, X, seed = 123, Kfolds = 3)
#'# kcv <- summary(kcv.out, digits = 3) 
#'# kcv$overview                   # test stats, in + out of sample, all folds
#'# kcv$training2.ttests                            # hypothesis tests, fold 2 
#'# save.bigKRLS(kcv.out, "myKfolds")
#'# load.bigKRLS("/path/to/myKfolds")     
#' @export 
crossvalidate.bigKRLS <- function(y, X, seed, Kfolds = NULL, ptesting = NULL, estimates_subfolder = NULL, ...){
  
  if(is.null(Kfolds) + is.null(ptesting) != 1) stop("Specify either Kfolds or ptesting but not both.")
  
  # suppressing warnings from bigmatrix
  oldw <- getOption("waxfrn")
  options(warn = -1)
  options(bigmemory.allow.dimnames=TRUE)
  
  stopifnot(is.big.matrix(X) | is.matrix(X))
  
  arguments <- list(...)
  marginals <- TRUE
  if("derivative" %in% names(arguments)){
    marginals <- arguments[["derivative"]]
  } # flag: compute test stats that require derivatives?
  
  Noisy <- nrow(X) > 2000 
  if("noisy" %in% names(arguments)){
    Noisy <- arguments[["noisy"]]
  } # flag: make CV output in line with user wishes, bigKRLS() defaults
  
  set.seed(seed)
  N <- nrow(X)
  
  big.meta <- create.metadata.dir()
  
  if(!is.null(ptesting)){
    
    if(ptesting < 0 | ptesting > 100) stop("ptesting, the percentage of data to be used for validation, must be between 0 and 100.")
    
    Ntesting <- round(N * ptesting/100, 0)
    Ntraining <- N - Ntesting
    train.set <- sample(N, Ntraining, replace = FALSE)
    test.set <- matrix(1:N)[which(!(1:N %in% train.set))]
    
    Xtrain <- submatrix(X, train.set)
    Xtest <- submatrix(X, test.set)
    ytrain <- submatrix(y, train.set)    
    ytest <- submatrix(y, test.set)
    
    trained <- bigKRLS(ytrain, Xtrain, instructions = FALSE, ...)
    tested <- predict.bigKRLS(trained, Xtest)
    tested[["ytest"]] <- ytest
    
    cv_out <- list(trained = trained, tested = tested, type = "crossvalidated")
    cv_out[["seed"]] <- seed
    cv_out[["indices"]] <- list(train.set = train.set, test.set = test.set)
    cv_out[["pseudoR2_is"]] <- trained$R2
    cv_out[["pseudoR2_oos"]] <- cor(tested$predicted, ytest[])^2
    cv_out[["MSE_oos"]] <- mean((tested$predicted - ytest[])^2)
    cv_out[["MSE_is"]] <- mean((trained$yfitted - trained$y[])^2)
    
    if(marginals){
      
      cv_out[["pseudoR2AME_is"]] <- trained$R2AME
      
      delta <- if(is.big.matrix(trained$X)) 
        to.big.matrix(matrix(trained$avgderivatives), p = 1, path = big.meta) else
          t(trained$avgderivatives)
      cv_out[["MSE_AME_is"]] <- mean((trained[["y"]] - (trained[["X"]] %*% delta)[])^2)
      
      delta <- if(is.big.matrix(Xtest)) 
        to.big.matrix(matrix(trained$avgderivatives), p = 1, path = big.meta) else
          t(trained$avgderivatives)
      yhat_ame <- (Xtest %*% delta)[]
      cv_out[["pseudoR2AME_oos"]] <- cor(tested[["ytest"]][], yhat_ame)^2
      cv_out[["MSE_AME_oos"]] <- mean((ytest - (Xtest %*% delta)[])^2)
      
    }
    
    cv_out[["ptesting"]] <- ptesting
    class(cv_out) <- "bigKRLS_CV" 
    # one bigKRLS object, one bigKRLS.predict object, type either "crossvalidated" or "kfolds"
    if("big.matrix" %in% lapply(trained, class) & is.null("estimates_subfolder")) 
      cat("NOTE: Outputted object contains big.matrix objects. To avoid crashing R, use save.bigKRLS(), not base R save() to store results.")
    if(Noisy) cat("You may wish to use summary() on the outputted object.")
    return(cv_out)
    
  }
  
  if(!is.null(Kfolds)){
    
    stopifnot(is.numeric(Kfolds) & Kfolds > 0 & Kfolds %% 1 == 0)
    
    # randomly places observations into (approximately) equal folds
    folds <- as.integer(cut(sample(N), breaks = Kfolds))
    
    for(k in 1:Kfolds){
      
      if(Noisy) cat("\n\n Performing pre-check of data for fold ", k, ".\n\n", sep="")
      
      Xtrain <- submatrix(X, folds != k)
      ytrain <- submatrix(y, folds != k)
      
      check_data(ytrain, Xtrain, instructions = FALSE)
      
    }
    
    out <- list(type = "KfoldsCV") # object to be returned
    class(out) <- "bigKRLS_CV"
    # out contains measures of fit, meta data, and (nested within each nested fold) bigKRLS, predict objects
    out[["Kfolds"]] <- Kfolds
    out[["seed"]] <- seed
    out[["folds"]] <- folds
    warn.big <- FALSE # dummy flag variable: warn re: big.matrix objects?
    
    # K measures of fit for each fold...
    out[["R2_is"]] <- c() # in sample R2, based on yhat = kernel[train, ] %*% coefs.hat
    out[["R2_oos"]] <- c() # out of sample R2, based on kernel(X[cbind(test, train])
    out[["MSE_is"]] <- c() # in sample mean squared error
    out[["MSE_oos"]] <- c() # out of sample mean squared error
    
    if(marginals){
      
      out[["R2AME_is"]] <- c() # in sample R2 average marginal effects, yhat = X[train, ] %*% colMeans(delta)
      out[["R2AME_oos"]] <- c() # oos R2, yhat = X[test, ] %*% colMeans(delta)
      out[["MSE_AME_is"]] <- c() # is for MSE for AMEs
      out[["MSE_AME_oos"]] <- c() # oos for MSE for AMEs
      
    }
    
    
    for(k in 1:Kfolds){
      
      if(Noisy) cat("\n\n Starting fold ", k, ".\n\n", sep="")
      
      Xtrain <- submatrix(X, folds != k)
      Xtest <- submatrix(X, folds == k)
      ytrain <- submatrix(y, folds != k)    
      ytest <- submatrix(y, folds == k)
      
      trained <- bigKRLS(ytrain, Xtrain, instructions = FALSE, ...)
      if(Noisy) summary(trained)
      tested <- predict.bigKRLS(trained, Xtest)
      tested[["ytest"]] <- ytest
      
      cv_out <- list(trained = trained, tested = tested)
      class(cv_out) <- "bigKRLS_CV"
      # cv_out contains one bigKRLS object (trained), one bigKRLS.predict object (tested)
      
      out[[paste0("fold_", k)]] <- cv_out
      
      ytest <- as.matrix(ytest) # in case of big.matrix objects...
      ytrain <- as.matrix(ytrain) # (loading vectors into R generally harmless)
      
      # measures of fit...
      out[["R2_is"]][k] <- trained$R2
      out[["R2_oos"]][k] <- cv_out[["tested"]][["pseudoR2"]] <- cor(ytest, as.matrix(tested$predicted))^2
      out[["MSE_is"]][k] <- cv_out[["trained"]][["MSE"]] <- mean((ytrain - as.matrix(trained$yfitted))^2)
      out[["MSE_oos"]][k] <- cv_out[["tested"]][["MSE"]] <- mean((ytest - tested$predicted)^2)
      
      if(marginals){
        
        out[["R2AME_is"]][k] <- trained$R2AME
        
        delta <- if(is.big.matrix(trained$X)) 
          to.big.matrix(matrix(trained$avgderivatives), p = 1, path = big.meta) else
            t(trained$avgderivatives)
        out[["MSE_AME_is"]][k] <- cv_out[["trained"]][["MSE_AME"]] <- mean((trained[["y"]] - (trained[["X"]] %*% delta)[])^2)
        
        delta <- if(is.big.matrix(Xtest)) 
          to.big.matrix(matrix(trained$avgderivatives), p = 1, path = big.meta) else
            t(trained$avgderivatives)
        yhat_ame <- (Xtest %*% delta)[]
        out[["R2AME_oos"]][k] <- cor(ytest, yhat_ame)^2
        out[["MSE_AME_oos"]][k] <- cv_out[["tested"]][["MSE_AME"]] <- mean((ytest - yhat_ame)^2)
      }
      
      warn.big <- warn.big | ("big.matrix" %in% lapply(trained, class) & is.null("estimates_subfolder")) 
      cat("\n")  
    }
    
    names(out[["R2_is"]]) <- names(out[["R2_oos"]]) <- names(out[["MSE_is"]]) <- 
      names(out[["MSE_oos"]]) <- paste0("fold", 1:Kfolds)
    
    if(marginals){
      names(out[["R2AME_is"]]) <- names(out[["R2AME_oos"]]) <- 
        names(out[["MSE_AME_is"]]) <- names(out[["MSE_AME_oos"]]) <- paste0("fold", 1:Kfolds)
    }
    
    if(warn.big) cat("NOTE: Outputted object contains big.matrix objects. To avoid crashing R, use save.bigKRLS(), not base R save() to store results.")
    
    if(!is.null(estimates_subfolder)) save.bigKRLS(out)
    unlink(big.meta, recursive = TRUE)
    
    return(out)
    
  }
  
}
