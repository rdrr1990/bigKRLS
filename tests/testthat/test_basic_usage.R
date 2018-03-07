# test the basic usage in the package documentation works
context("Basic usage of bigKRLS")

# prep data to use for testing
mtcars <- datasets::mtcars
y <- as.matrix(mtcars$mpg)
X <- as.matrix(mtcars[,-1])

test_equivalent_models <- function(mod1, mod2) {
  # ignore any path differences
  mod1$path <- mod2$path <- NULL
  # ignore differences in name sorting 
  expect_identical(sort(names(mod1)), sort(names(mod2)))
  mod2 <- mod2[names(mod1)]
  
  # check that most elements are identical, except for expected handful
  identical_elements <- mapply(identical, mod1, mod2[names(mod1)])
  names_to_check <- sort(names(which(!identical_elements)))
  allowed_not_identical <- c("derivatives", "K", "path", "vcov.est.c", "X")
  expect_true(all(names_to_check %in% allowed_not_identical))
  
  # allow some fields to just be nearly equal
  allowed_nearly_equal <- c("derivatives", "K", "vcov.est.c", "X")
  for(v in intersect(names_to_check, allowed_nearly_equal)) {
    expect_equivalent(as.matrix(mod1[[v]]), as.matrix(mod2[[v]]), info = paste0("v = ",v))
  }
  names_to_check <- setdiff(names_to_check, allowed_nearly_equal)
  
  # allow some fields to be different
  names_to_check <- setdiff(names_to_check, "path")
  
  # should have nothing left to check
  expect_identical(names_to_check, character(0), info=paste0("names_to_check = ", paste(names_to_check, collapse=", ")))
  
  invisible(NULL)
}


test_that("Simple example works", {
  # fitting
  reg.out <- bigKRLS(y = y, X = X, eigtrunc=0, Ncores = 1)
  summary(reg.out)
  
  # saving/loading with normal matrices
  model_subfolder <- "bigKRLS_test_results"
  save.bigKRLS(reg.out, model_subfolder, overwrite.existing = TRUE)
  reg.out2 <- load.bigKRLS(model_subfolder, pos = NULL, return_object = TRUE)
  reg.out2[["model_subfolder_name"]] <- NULL
  # test that saved object is equivalent
  test_equivalent_models(reg.out, reg.out2)
  # remove subfolder
  unlink(model_subfolder, recursive = TRUE)
  
  # prediction
  Xnew <- datasets::mtcars[,-1]
  Xnew$hp <- 200
  forecast <- predict(reg.out, as.matrix(Xnew))
  expect_equal(mean(forecast$predicted < mtcars$mpg), 0.6875)
  
  # similarity
  
  s <- reg.out$K[, grep("Corolla", rownames(mtcars))]
  names(s) <- rownames(mtcars)
  s <- s[order(names(s))]
  s_expected <- c('AMC Javelin'=0.0547298949171582, 
                  'Cadillac Fleetwood'=0.00549165470976291, 
                  'Camaro Z28'=0.0156630175526991, 
                  'Chrysler Imperial'=0.0060180975553816, 
                  'Datsun 710'=0.860610665218997, 
                  'Dodge Challenger'=0.033400030235352, 
                  'Duster 360'=0.0143264812794483, 
                  'Ferrari Dino'=0.062192422562695, 
                  'Fiat 128'=0.973400786036153, 
                  'Fiat X1-9'=0.961130622208994, 
                  'Ford Pantera L'=0.0207382308766512, 
                  'Honda Civic'=0.753451355337079, 
                  'Hornet 4 Drive'=0.19371687432462, 
                  'Hornet Sportabout'=0.0388127837578353, 
                  'Lincoln Continental'=0.00503976771060228, 
                  'Lotus Europa'=0.528183252015446, 
                  'Maserati Bora'=0.00201340749064979, 
                  'Mazda RX4'=0.239466325088983, 
                  'Mazda RX4 Wag'=0.254841103009284, 
                  'Merc 230'=0.373560094613131, 
                  'Merc 240D'=0.464081081884477, 
                  'Merc 280'=0.250345020593959, 
                  'Merc 280C'=0.262879139614823, 
                  'Merc 450SE'=0.0344532182858226, 
                  'Merc 450SL'=0.0411135560575867, 
                  'Merc 450SLC'=0.0424741434773812, 
                  'Pontiac Firebird'=0.0270102739090449, 
                  'Porsche 914-2'=0.3709635022494, 
                  'Toyota Corolla'=1, 
                  'Toyota Corona'=0.468060548244946, 
                  'Valiant'=0.146078393891752, 
                  'Volvo 142E'=0.78179636900690)
  expect_identical(names(s), names(s_expected))
  s_difference <- s - s_expected
  expect_lt(max(s_difference), 0.01)
})


test_that("bigmemory example works", {
  model_subfolder <- "bigKRLS_test_bigmemory_results"
  big <- bigKRLS(
    y = as.big.matrix(y),
    X = as.big.matrix(X), Neig = nrow(X),
    Ncores = 1, 
    model_subfolder_name = model_subfolder,
    overwrite.existing = TRUE
  )
  
  # compare saved model and loaded model
  big2 <- load.bigKRLS(model_subfolder, pos = NULL, return_object = TRUE)
  
  Kdiff <- max(big$K - big2$K)
  expect_lt(Kdiff, 0.000001)
  
  expect_equal(big$yfitted, big2$yfitted)
  
  unlink(model_subfolder, recursive = TRUE)
  
})

test_that("crossvalidation function works", {
  
  cv <- crossvalidate.bigKRLS(y, X, ptesting = 20, seed = 123, Neig = .8*nrow(X), Ncores = 1)
  summary(cv)
  expect_lt(cv$pseudoR2_oos, cv$pseudoR2_is)
  
  cv_noderivs <- crossvalidate.bigKRLS(y, X, ptesting = 20, seed = 123,
                                       Neig = .8*nrow(X), Ncores = 1, derivative = FALSE)
  expect_equal(cv$pseudoR2_oos, cv_noderivs$pseudoR2_oos)
  
})

set.seed(1234)
X <- matrix(runif(1000), nrow = 250, ncol = 4)
y <- X %*% 1:4 + rnorm(250)

kcv <- crossvalidate.bigKRLS(y, X, Kfolds = 4, seed = 1234, Neig = 0.75*nrow(X), Ncores = 1)
kcv_noderivs <- crossvalidate.bigKRLS(y, X, Kfolds = 4, seed = 1234, Neig = 0.75*nrow(X),Ncores = 1, derivative = FALSE)
kcvbig <- crossvalidate.bigKRLS(y, as.big.matrix(X), Kfolds = 4, seed = 1234, Neig = 0.75*nrow(X),Ncores = 1)

test_that("Kfolds crossvalidation works", {
  
  summary(kcv)
  expect_equal(kcv$folds, kcvbig$folds)
  expect_equal(sum(kcv$fold_1$tested$newdata), sum(kcvbig$fold_1$tested$newdata[]))
  expect_equal(sum(kcv$fold_2$tested$predicted), sum(kcvbig$fold_2$tested$predicted[]))
  expect_equal(kcv$pseudoR2_oos, kcv_noderivs$pseudoR2_oos)
  
})

test_that("Kfolds test stats, big vs base (batch 1)", {
  
  expect_equal(kcv$R2_is, kcvbig$R2_is)
  expect_equal(kcv$R2_oos, kcvbig$R2_oos)
  expect_equal(kcv$R2AME_is, kcvbig$R2AME_is)
  expect_equal(kcv$MSE_is, kcvbig$MSE_is)

})

test_that("Kfolds test stats, big vs base (batch 2)", {

  expect_equal(kcv$MSE_oos, kcvbig$MSE_oos)
  expect_equal(kcv$MSE_AME_is, kcvbig$MSE_AME_is)
  expect_equal(kcv$R2AME_oos, kcvbig$R2AME_oos)
  expect_equal(kcv$MSE_AME_oos, kcvbig$MSE_AME_oos)
  
})
