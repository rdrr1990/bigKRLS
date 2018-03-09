####################################
# functions that support bigKRLS() #
####################################

bLambdaSearch <- function (L = NULL, U = NULL, y = NULL, Eigenobject = NULL, tol = NULL, 
                           noisy = FALSE){
  
  if(sum(is.na(Eigenobject$values)) > 0) 
    stop("Missing eigenvalues prevent bigKRLS from obtaining the regularization parameter lambda.\n\tCheck for repeated observations (or other perfect linear combinations in X).")
  n <- nrow(y)
  if (is.null(tol)) {
    tol <- 10^-3 * n # tolerance parameter
  } else {
    stopifnot(is.vector(tol), length(tol) == 1, is.numeric(tol), tol > 0)
  }
  if (is.null(U)) {
    U <- n
    #if(Eigenobject$lastkeeper == n){
      while (sum(Eigenobject$values/(Eigenobject$values + U)) < 1) {
        U <- U - 1
      }
    #}
    #  } else {
    stopifnot(is.vector(U), length(U) == 1, is.numeric(U), U > 0)
  }
  if (is.null(L)) {
    
    L <- .Machine$double.eps # smallest double such that 1 + x != 1. Normally 2.220446e-16.
    #if(Eigenobject$lastkeeper == n){
      q <- which.min(abs((Eigenobject$values - max(Eigenobject$values)/1000)))
      while (sum(Eigenobject$values/(Eigenobject$values + L)) > q) {
        L <- L + 0.05 
      } 
    #}
  } #else {
    stopifnot(is.vector(L), length(L) == 1, is.numeric(L), L >= 0)
  #}
  X1 <- L + (0.381966) * (U - L) 
  X2 <- U - (0.381966) * (U - L)
  
  # bLooLoss is big Leave One Out Error Loss
  
  if(noisy) cat("\n\nGetting S1.")
  S1 <- bLooLoss(lambda = X1, y = y, Eigenobject = Eigenobject)
  if(noisy){cat("\nGetting S2.")}
  S2 <- bLooLoss(lambda = X2, y = y, Eigenobject = Eigenobject)
  f3 <- function(x){format(round(x, digits=3), nsmall=3)}
  if (noisy) {
    cat("\n\nL: ", f3(L), 
        " X1: ", f3(X1), " X2: ", f3(X2), 
        " U: ", f3(U), " S1: ", f3(S1), " S2: ", f3(S2), 
        "\n(", Time(),  
        ").", sep = "")
  }
  while (abs(S1 - S2) > tol) {
    if (S1 < S2) {
      U <- X2
      X2 <- X1
      X1 <- L + (0.381966) * (U - L)
      S2 <- S1
      S1 <- bLooLoss(lambda = X1, y = y, Eigenobject = Eigenobject)
    }
    else {
      L <- X1
      X1 <- X2
      X2 <- U - (0.381966) * (U - L)
      S1 <- S2
      S2 <- bLooLoss(lambda = X2, y = y, Eigenobject = Eigenobject)
    }
    if (noisy) {
      cat("\nL: ", f3(L), 
          " X1: ", f3(X1), " X2: ", f3(X2), 
          " U: ", f3(U), " S1: ", f3(S1), " S2: ", f3(S2), 
          "\n(", Time(),  
          ").", sep = "")    }
  }
  out <- ifelse(S1 < S2, X1, X2)
  
  if (noisy) {cat("\n\nlambda = ", round(out, 5), ".\n", sep='')}
  
  return(invisible(out))
}

bSolveForc <- function (y = NULL, Eigenobject = NULL, lambda = NULL) {
  
  out <- BigSolveForc(Eigenobject$vectors@address, 
                        Eigenobject$values, y[], lambda) #, Eigenobject$lastkeeper)
  return(list(Le = out[[1]], coeffs = out[[2]]))
  
}

bLooLoss <- function (y = NULL, Eigenobject = NULL, lambda = NULL) 
{
  return(bSolveForc(y = y, Eigenobject = Eigenobject, lambda = lambda)$Le)
} 


#######################################
# Rcpp and bigmemory Helper Functions #
#######################################

noise <- function() sample(as.numeric(gsub("\\.", "", as.character(as.numeric(format(Sys.time(), "%OS"))))), 1)
Time <- function() format(Sys.time(), format = "%H:%M:%S")

create.metadata.dir <- function(){
  
  tmp <- tempdir()
  Nextant <- length(dir(path = tmp, pattern = "filebacks"))
  big.meta <- file.path(tmp, paste0("filebacks", Nextant + 1, "_", noise()))   
  dir.create(big.meta)
  return(big.meta)
  
}


to.big.matrix <- function(object, p = NULL, deepcopy = FALSE, name = NULL, path = NULL){
  
  # returns file-backed big.matrix 
  # (takes `object` as matrix, matrix-like object, or big.matrix)
  # ensures vectors are column matrices
  # coerces integers to numeric so that they are properly cast as doubles on Rcpp side
  # name is optional. 
  # name = "X" stores info on the bigmatrix being created as "X.desc".
  # name = NULL leads to "tmp.desc" (or tmp1, tmp2, etc. as necessary)
  #   
  # optionally returns deep copy  
  
  # p <- if(is.null(p) && is.null(ncol(object))) 1 else ncol(object)
  
  if(is.null(p)){
    p <- ifelse(is.null(ncol(object)), 1, ncol(object))
  } 
  
  if(is.null(path)){
    path <- create.metadata.dir()
  }
  
  if(is.null(name)) {
    #name <- paste0("bigmatrix", length(dir(path = path, pattern = "desc")) + 1, "_", atomic())
    name <- basename(tempfile(tmpdir = path))
  }
  
  if(!is.big.matrix(object)){
    object <- as.big.matrix(matrix(as.numeric(object), ncol = p), 
                            backingfile = name, 
                            backingpath = path,
                            descriptorfile = paste0(name, ".desc"))
  }else{
    
    # use filebacked function instead ?
    dput(describe(object), file.path(path, paste0(name, ".desc")))
    
  }
  if(deepcopy) return(deepcopy(object)) else return(object)
  
}


bMultDiag <- function (X, v, check_platform = FALSE) {
  
  if(check_platform) check_platform()
  # multdiag.cpp
  out <- big.matrix(nrow = nrow(X),
                    ncol = ncol(X),
                    init = 0,
                    type = 'double')
  
  BigMultDiag(X@address, v, out@address)
  
  return(out)
}

bEigen <- function(A, Neig = NULL, eigtrunc = 0, check_platform = FALSE){
  # A, typically the kernel, is assumed symmetric by underlying functions
  
  if(check_platform) check_platform()
  
  Neig <- if(is.null(Neig)) nrow(A) else Neig
  # corresponds to arma::eig_sym and arma::eigs_sym, respectively
  
  vals <- big.matrix(nrow = 1, ncol = Neig, type = 'double')
  vecs <- to.big.matrix(big.matrix(nrow = nrow(A), ncol = Neig, type = 'double'), name = "vecs")
  
  BigEigen(A@address, Neig, vals@address, vecs@address)
  
  vecs <- -1*vecs
  
  out <- list('values' = vals[])
  
  out[["lastkeeper"]] <- max(which(out$values >= eigtrunc*out$values[1]))
  
  if(out$lastkeeper == ncol(vecs)){
    out[["vectors"]] <- vecs
  }else{
    out[["vectors"]] <- deepcopy(vecs, cols = 1:out[["lastkeeper"]])
    remove(vecs)
  } 
  return(out)
}

bGaussKernel <- function(X, bandwidth = NULL, check_platform = FALSE){ 
  
  if(check_platform) check_platform()
  
  bandwidth <- if(is.null(bandwidth)) ncol(X) else bandwidth
  out <- big.matrix(nrow=nrow(X), ncol=nrow(X), init=0)
  BigGaussKernel(X@address, out@address, bandwidth) # gauss_kernel.cpp
  
  return(out)
}

bNeffective <- function(X, check_platform = FALSE){
  
  if(check_platform) check_platform()
  # Neffective.cpp
  return(BigNeffective(X@address))
}

bTempKernel <- function(X_new, X_old, sigma, check_platform = FALSE){
  
  if(check_platform) check_platform()
  # temp_kernel.cpp
  out <- big.matrix(nrow=nrow(X_new), ncol=nrow(X_old), init=0)
  
  BigTempKernel(X_new@address, X_old@address, out@address, sigma)
  return(out)
}

bCrossProd <- function(X, Y=NULL, check_platform = FALSE){
  
  if(check_platform) check_platform()
  
  out_ncol <- if(is.null(Y)) ncol(X) else ncol(Y)

  out <- big.matrix(nrow = ncol(X), ncol = out_ncol, type = 'double')
  
  if(is.null(Y)){
    BigXtX(X@address, out@address)
  }else{
    BigCrossProd(X@address, Y@address, out@address)
  }
  return(out)
}

bTCrossProd <- function(X, Y = NULL, check_platform = FALSE){
  
  if(check_platform) check_platform()
  
  out_ncol <- if(is.null(Y)) nrow(X) else nrow(X)
  
  out <- big.matrix(nrow = nrow(X), ncol = out_ncol, type = 'double')
  if(is.null(Y)){
    BigXXt(X@address, out@address)
  }else{
    BigTCrossProd(X@address, Y@address, out@address)
  }
  return(out)
}

bDerivatives <- function(X, sigma, K, coeffs, vcovmatc, check_platform = FALSE){
  
  if(check_platform) check_platform()
  derivatives <- big.matrix(nrow=nrow(X), ncol=ncol(X), init=-1)
  varavgderiv <- big.matrix(nrow=1, ncol=ncol(X), init=-1)
  out <- BigDerivMat(X@address, K@address, vcovmatc@address, 
                     derivatives@address, varavgderiv@address,
                     coeffs, sigma)
  
  return(list('derivatives'= derivatives, 'varavgderiv' = varavgderiv[]))
}


make_path <- function(object, model_subfolder_name, overwrite.existing){
  
  # thanks to Peter Foley for helpful suggestions re: file and folder management!
  # see pulls 11-13 starting with https://github.com/rdrr1990/bigKRLS/pull/11
  
  if(!overwrite.existing && dir.exists(model_subfolder_name)){
    i <- 1
    tmp.name <- paste(model_subfolder_name, i, sep="")
    while(tmp.name %in% dir()){
      tmp.name <- paste(model_subfolder_name, i, sep="")
      i <- i + 1
    }
    if(model_subfolder_name %in% dir())
      warning(cat("A subfolder named", model_subfolder_name, "exists in your current working directory. Your output will be saved to", tmp.name, "instead. To turn off this safeguard, set save.bigKRLS(..., overwrite.existing = TRUE) next time.\n\n"))
    model_subfolder_name <- tmp.name
  }
  
  dir.create(model_subfolder_name, recursive = TRUE, showWarnings = FALSE)
  if(dir.exists(model_subfolder_name)) 
    cat("Saving model estimates to:\n\n", model_subfolder_name, "\n\n") else 
      stop("Unable to create directory.")
  object[["path"]] <- normalizePath(model_subfolder_name)
  object[["model_subfolder_name"]] <- model_subfolder_name
  return(object)
  
}

bSave <- function(object, noisy){
  
  is.big.mat <- unlist(lapply(object, is.big.matrix))
  
  for(i in which(is.big.mat)){
    output_path <- file.path(object[["model_subfolder_name"]], paste0(names(object)[i], ".txt"))
    if(noisy) cat("\twriting", output_path, "...\n")
    write.big.matrix(x = object[[i]], # col.names = !is.null(colnames(object[[i]])),
                     filename = output_path)
  }
  
  Nbm <- sum(is.big.mat)
  if(noisy) cat("\n", Nbm, " matrices saved as big matrices", 
                ifelse(Nbm == 0, " (base R save() may be used safely in this case too).\n",
                       ", use load.bigKRLS() on the entire directory to reconstruct the outputted object in R.\n"), sep="")
  if(Nbm > 0){
    bigKRLS_out <- object[-which(is.big.mat)]
    class(bigKRLS_out) <- class(object) 
  }else{
    bigKRLS_out <- object
  }
  remove(object)
  stopifnot(sum(unlist(lapply(bigKRLS_out, is.big.matrix))) == 0)
  save(bigKRLS_out, 
       file = file.path(bigKRLS_out[["model_subfolder_name"]], "estimates.RData"))
  if(noisy) cat("Smaller, base R elements of the outputted object saved:", 
                file.path(bigKRLS_out[["model_subfolder_name"]], "estimates.RData"), "\n")
  
}

bLoad <- function(object, path, noisy){
  
  options(bigmemory.allow.dimnames=TRUE)
  
  if(class(object) == "bigKRLS"){
    matrices <- c("K", "X", "derivatives", "vcov.est.c", "vcov.est.fitted")
  }else{
    if(class(object) == "bigKRLS_predicted"){
      matrices <- c("predicted", "se.pred", "vcov.est.pred", "newdata", "newdataK", "ytest")
    }else{
      stop("bLoad may only be used on bigKRLS objects or bigKRLS_predicted objects.")
    }
  }
  
  `%out%` <- function(x, table) match(x, table, nomatch = 0L) == 0L
  xlabs <- object$xlabs
  which.derivatives <- object$which.derivatives
  
  for(i in 1:length(matrices)){
    
    filename <- paste0(matrices[i], ".txt")
    
    if(filename %out% dir(path = path) & matrices[i] %out% names(object)){
      
      if(noisy) cat("NOTE:", matrices[i],  
                    "not found in .RData or in big matrix file,", 
                    filename,".\n\n")
    }else{
      if(filename %in% dir(path = path)){
        if(noisy) cat("\tReading from", filename, "\n")
        object[[matrices[i]]] <- read.big.matrix(file.path(path, filename), 
                                                 type = "double")
        
        if(ncol(object[[matrices[i]]]) == length(xlabs)){
          colnames(object[[matrices[i]]]) <- xlabs
        }else{
          if(!is.null(which.derivatives)){
            if(ncol(object[[matrices[i]]]) == length(xlabs[which.derivatives])){
              colnames(object[[matrices[i]]]) <- xlabs[which.derivatives]
            }
          }
        }
        stopifnot(is.big.matrix(object[[matrices[i]]]))
      }
      
    }
  }
  
  return(object)
}

# check_data() performs all the checks that bigKRLS() performs...
# it is intended for K folds crossvalidation. 
# categorical variables can be fussy when randomly partitioning...

check_data <- function (y = NULL, X = NULL, sigma = NULL, 
                        derivative = TRUE, which.derivatives = NULL,
                        vcov.est = TRUE, lambda = NULL, L = NULL, U = NULL, 
                        tol = NULL, model_subfolder_name = NULL, 
                        overwrite.existing = FALSE, Ncores = NULL, 
                        acf = FALSE, noisy = NULL, instructions = TRUE)
{
  
  # suppressing warnings from bigmatrix
  oldw <- getOption("warn")
  options(warn = -1)
  #  options(bigmemory.allow.dimnames=TRUE)
  
  stopifnot(is.matrix(X) | is.big.matrix(X))
  big.meta <- create.metadata.dir()
  
  X <- to.big.matrix(X, path = big.meta)
  X.init.sd <- colsd(X)
  y <- to.big.matrix(y, p = 1, path = big.meta)
  
  miss.ind <- colna(X)
  if (sum(miss.ind) > 0) { 
    stop(paste("the following columns in X contain missing data, which must be removed:", 
               paste((1:length(miss.ind))[miss.ind > 0], collapse = ', '), collapse=''))
  }
  n <- nrow(X)
  p <- ncol(X)
  
  if (min(X.init.sd) == 0) {
    stop(paste("The following columns in X are constant and must be removed:",
               which(X.init.sd == 0)))
  }
  
  if (n != nrow(y)) { stop("nrow(X) not equal to number of elements in y.")}
  if (colna(y) > 0) { stop("y contains missing data.") }
  if (colsd(y) == 0) { stop("y is a constant.") }
  
  unlink(big.meta, recursive = TRUE)
  
}

bDiag <- function(A){
  # returns diagonal of big.matrix as column vector (base R matrix)  
  
  d <- matrix(nrow = nrow(A), ncol = 1)
  for(i in 1:nrow(A)){
    d[i] <- deepcopy(A, cols = i, rows = i)[] # use sub.big.matrix?
  }
  
  return(d) # make to.big.matrix(d)
} 

submatrix <- function(X, rows){
  if(is.big.matrix(X)) deepcopy(X, rows = rows) else X[rows, ]
} # make sub.big.matrix the default. 
# add cols, to.big.matrix... check on behavior of deep copy with filebacked
