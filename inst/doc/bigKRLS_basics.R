## ---- eval=FALSE---------------------------------------------------------
#  library(devtools)
#  install_github('rdrr1990/bigKRLS')

## ---- eval=FALSE---------------------------------------------------------
#  library(bigKRLS)

## ---- results='hide', message=FALSE, warning=FALSE, echo=FALSE-----------
library(bigKRLS)
library(knitr)

## ---- echo=FALSE---------------------------------------------------------
kable(head(mtcars, 6))

## ------------------------------------------------------------------------
reg.out <- bigKRLS(y = as.matrix(mtcars$mpg), X = mtcars[,-1], noisy=FALSE)

## ------------------------------------------------------------------------
summary(reg.out)

## ------------------------------------------------------------------------
reg.out$K

## ---- fig.width = 7------------------------------------------------------
s <- reg.out$K[which(mtcars$cyl == 4), grep("Corolla", rownames(mtcars))]
barplot(s, main = "Similarity to a Toyota Corolla", 
        ylab = "Kernel", sub="Toy Data from mtcars",  cex.names = .7,
        col = colorRampPalette(c("red", "blue"))(length(s))[rank(s)],
        names.arg = lapply(strsplit(rownames(mtcars), split=" "), 
                           function(x) x[2])[which(mtcars$cyl == 4)])

## ---- fig.height=6, fig.width=7.5----------------------------------------

scatter.smooth(mtcars$hp, reg.out$derivatives[,3], ylab="HP's Effect", xlab="Horsepower", pch = 19, bty = "n",
               main="Horsepower's Marginal Effect on Fuel Efficiency",
               sub="Toy Data from mtcars",
               col = colorRampPalette(c("blue", "red"))(nrow(mtcars))[rank(reg.out$coeffs^2)], 
               ylim = c(-0.042, 0.015), xlim = c(50, 400))

fields::image.plot(legend.only = T, zlim=c(1/nrow(mtcars), 1), legend.cex = 0.9,   
           col = colorRampPalette(c("red", "blue"))(nrow(mtcars)), 
           legend.shrink = .75)
text(x = 380, y = 0.015, "Relative Fit")
text(x = 380, y = 0.012, "in Full Model")


## ------------------------------------------------------------------------
 cor(reg.out$coeffs, reg.out$y - reg.out$fitted)

## ---- eval=FALSE---------------------------------------------------------
#  reg.out$K <- reg.out$vcov.c <- reg.out$vcov.fitted <- NULL

## ------------------------------------------------------------------------
set.seed(1776)
N <- 1000  
P <- 4
X <- matrix(rnorm(N*P), ncol=P)
X[,P] <- ifelse(X[,P] > 0.12345, 1, 0)
b <- runif(ncol(X))
y <- X %*% b + rnorm(nrow(X))
bigKRLS.out <- bigKRLS(X = X, y = y, noisy = F)
summary(bigKRLS.out, digits=5)

## ------------------------------------------------------------------------
KRLS.out <- KRLS::krls(X = X, y = y, print.level = 0)
max(abs(bigKRLS.out$derivatives - KRLS.out$derivatives)) < 0.00000001

## ---- echo=FALSE---------------------------------------------------------
tmp.exp <- ceiling(log(max(abs(bigKRLS.out$derivatives - KRLS.out$derivatives)), base=10))

