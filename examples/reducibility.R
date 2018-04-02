reduciblility <- function(object, loss=2, q = 0.05){
  
  if(class(object) != "bigKRLS")
    warning("reducibility() was written for bigKRLS regression objects.")
  
  if(!(loss %in% 1:2))
    warning("Only L1 and L2 loss implemented, defaulting to L2.")
  
  N <- nrow(object$X)
  P <- ncol(object$X)
  
  # predicted values using only Average Marginal Effects
  yhat_ame <- object$X %*% t(object$avgderivatives)
  
  # each column of yhat_p contains a separate set of predictions
  # each uses the AMEs for all but one column of X
  # for the other, the 'actually' marginals are used (i.e., dy/dX_{i,p})
  yhat_p <- matrix(nrow=N, ncol=P)
  for(i in 1:P)
    yhat_p[,i] <- object$X[,-i] %*% object$avgderivatives[-i] + 
                  object$X[,i] * object$derivatives[,i]

  # the above code loads X into R's memory and should be redone with 
  # sub.big.matrix for big matrices if(object$has.big.matrices == TRUE)
  
  p_vals <- matrix(nrow=P)
  loss_null <- if(loss == 1) abs(object$yfitted - yhat_ame) else (object$yfitted - yhat_ame)^2
  
  for(i in 1:P){
    loss_alternative <- if(loss == 1) abs(object$y - yhat_p[,i]) else (object$y - yhat_p[,i])^2 
    p_vals[i] <- wilcox.test(loss_alternative, 
                             loss_null, paired = TRUE, 
                             alternative = "less")[["p.value"]] 
  }
     
  p_sorted <- sort(p_vals)
  i <- 1
  while(p_sorted[i] <= i*q/P)
    i <- i + 1

  # benjamini hochberg
  BH <- ifelse(p_vals <= p_sorted[i], "Reject Null", "Accept Null")
  rownames(BH) <- colnames(object$X)
  colnames(BH) <- "BH Decision"
  attr(BH, "description") <- "Null Hypothesis: the Average Marginal Effect approximates the regularized target function (y* = Kc*) as well as the Marginal Effects measured at each point (Benjamini-Hochberg decision based on Wilcox Rank Sum Test)."
  attr(BH, "loss") <- paste0("L", loss)
  object[["Reducibility"]] <- BH
  
  return(object)
}

library(MASS)
y <- as.matrix(Boston$medv)
X <- as.matrix(Boston %>% dplyr::select(-medv))

boston_out <- bigKRLS(y, X)
boston_out <- reduciblility(boston_out)
table(boston_out$Reducibility)
# Accept Null Reject Null 
#     9           4 

replication_data <- read.table("https://raw.githubusercontent.com/rdrr1990/bigKRLS/master/examples/data2016GEcensus.csv", header = TRUE, sep=",")
# reading from github seems to read row numbers as the first column (called 'X')

y <- as.matrix(replication_data$GOPdelta)
X <- as.matrix(replication_data)[, 3:ncol(replication_data)]

election_out <- bigKRLS(y, X, Neig=50)
election_out <- reduciblility(election_out)
election_out$Reducibility # rejects all

election_out <- reduciblility(election_out, loss=1)
election_out$Reducibility # rejects all



library(mvtnorm)
N <- 500
P <- 10

Sigma <- function(P, unit_variance = TRUE){
  A <- matrix(runif(P^2)*2-1, ncol=P)
  Sigma <- t(A) %*% A
  return(if(unit_variance) cov2cor(Sigma) else Sigma)
}

Nsims <- 10
for(i in 1:Nsims){
  X <- rmvnorm(N, sigma = Sigma(P))
  y <- X %*% rnorm(P) + X[,1]^2 + X[,1]^3 + rnorm(N)
  out <- reduciblility(bigKRLS(y,X,instructions = FALSE))
  print(table(out$Reducibility))
}

  
Nsims <- 25
results005 <- data.frame(true_disc = numeric(),
                  false_disc = numeric(),
                  false_accept = numeric(),
                  true_accept = numeric())

results010 <- data.frame(true_disc = numeric(),
                         false_disc = numeric(),
                         false_accept = numeric(),
                         true_accept = numeric())

which.cubic <- 1:5
for(i in 1:Nsims){
  
  X <- matrix(rnorm(N*P), ncol=P)
  y <- rnorm(N)
  for(i in which.cubic)
    y <- y + X[,i] + X[,i]^2 + X[,i]^3
  y <- y + X[,-which.cubic] %*% 1:5
  
  out <- bigKRLS(y,X,instructions = FALSE)
  out <- reduciblility(out)

  results005 <- rbind(results005, 
                   data.frame(true_disc = sum(out$Reducibility[which.cubic] == "Reject Null")/P,
                              false_disc = sum(out$Reducibility[-which.cubic] == "Reject Null")/P,
                              false_accept = sum(out$Reducibility[which.cubic] == "Accept Null")/P,
                              true_accept = sum(out$Reducibility[-which.cubic] == "Accept Null")/P))
  
  out <- reduciblility(out, q=0.1)
  results010 <- rbind(results010, 
                      data.frame(true_disc = sum(out$Reducibility[which.cubic] == "Reject Null")/P,
                                 false_disc = sum(out$Reducibility[-which.cubic] == "Reject Null")/P,
                                 false_accept = sum(out$Reducibility[which.cubic] == "Accept Null")/P,
                                 true_accept = sum(out$Reducibility[-which.cubic] == "Accept Null")/P))
  
}

colMeans(results005)
colMeans(results010)
