library(bigKRLS)

edata <- read.csv("https://raw.githubusercontent.com/rdrr1990/bigKRLS/master/examples/data2016GE.csv")
y <- as.matrix(edata[,1])
X <- as.matrix(edata[,-1])
dim(X)
eig <- bigKRLS:::bEigen(bigKRLS:::bGaussKernel(as.big.matrix(X)))
sum(eig$values)/nrow(X)
max(which(eig$values > 0.001*eig$values[1]))
# 2895

full <- bigKRLS(y,X[,1:18],eigtrunc=0, derivative = FALSE)

M <- 4
group <- sample(M, nrow(X), replace=TRUE)
coeffs.sub <- matrix(ncol=1, nrow=nrow(X))
eigvals <-list()
for(i in 1:M){
  out <- bigKRLS(y[group == i], X[group == i,1:18], 
                eigtrunc = 0, 
                derivative = FALSE, noisy=FALSE, instructions = FALSE)
  coeffs.sub[group == i] <- out[["coeffs"]]
  eigvals[[i]] <- out[["K.eigenvalues"]]
}
cor(coeffs.sub, full$coeffs)
cor(coeffs.sub, full$coeffs, method = "kendall")
cor(coeffs.sub, full$coeffs, method = "spearman")

for(i in 1:M){
  k <- length(eigvals[[i]])
  print(cor(full$K.eigenvalues[1:k], eigvals[[i]]))
}
