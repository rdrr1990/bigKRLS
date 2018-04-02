library(bigKRLS)
edata <- read.csv("https://raw.githubusercontent.com/rdrr1990/bigKRLS/master/examples/data2016GE.csv")

y <- edata[,1]
X <- cbind(edata[,2:16], 
           princomp(edata[,17:69])[["scores"]])
# rotate state dummies along with county lat/lon to avoid constant columns when subsampling


# runtime without derivative calculations
Neigen <- c(1, 2, 4, 8, 16, 32, ncol(X), nrow(X)) 
N <- c(500, 1000, 1500, 2000, 2500, 3000)
runtime <- matrix(nrow=length(Neigen)*length(N), ncol=5)
colnames(runtime) <- c("N", "Neigen", "user", "system", "elapsed")

for(i in 1:length(N)){
  for(j in 1:length(Neigen)){
    
    r <- j + length(Neigen)*(i-1)
    runtime[r, 1] <- N[i]
    runtime[r, 2] <- if(Neigen[j] < N[i]) Neigen[j] else N[i]
    # when Neigen[j] is nrow(X) == 3106,  
    # bigKRLS defaults to calculating whatever N[i] is
    
    subsamp <- sample(nrow(X), N[i])
    runtime[r, 3:5] <- system.time(bigKRLS(as.matrix(y[subsamp]), as.matrix(X[subsamp,]), 
                                        Neig=Neigen[j], 
                                        instructions = FALSE, derivative = FALSE))[1:3]
    cat("\n")
    print(runtime[1:r,])
  }
}

library(mvtnorm)
N <- 1000
P <- 50

Sigma <- function(P, unit_variance = TRUE){
  A <- matrix(runif(P^2)*2-1, ncol=P)
  Sigma <- t(A) %*% A
  return(if(unit_variance) cov2cor(Sigma) else Sigma)
}

X <- rmvnorm(N, sigma = Sigma(P))
K <- KRLS::gausskernel(X, P)
c <- rnorm(N)
lambda <- rexp(1)
beta <- rnorm(P)
y <- K %*% c + drop(lambda * crossprod(c, K %*% c))

out_eig <- bigKRLS(y,X)

out_N_eig <- bigKRLS(y,X,Neig=N)  
out_N_eig_001 <- bigKRLS(y,X)  
out_P2_eig_001 <- bigKRLS(y,X,Neig=P^2, eigtrunc = 0.01)  

system.time(bigKRLS(y,X,instructions = FALSE))
system.time(bigKRLS(y,X,eigtrunc=0,instructions = FALSE))

out_N_eig$R2
out_P_eig$R2

out_P_eig$R2AME
out_N_eig$R2AME
out_N_eig_001$R2
out_P2_eig_001$R2


cor(cbind(c, out_P_eig$coeffs, out_N_eig$coeffs))



out_P_eig$avgderivatives
out_N_eig$avgderivatives
beta

cor(t(out_P_eig$avgderivatives),beta)
cor(t(out_N_eig$avgderivatives),beta)

