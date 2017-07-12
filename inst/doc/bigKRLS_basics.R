## ---- echo=F, message=F, warning=F---------------------------------------
library(bigKRLS)

## ------------------------------------------------------------------------
mtcars[1:5,]

## ---- warning=F, message=F-----------------------------------------------
reg.out <- bigKRLS(y = as.matrix(mtcars$mpg), 
                   X = as.matrix(mtcars[,-1]), Ncores = 1)

## ------------------------------------------------------------------------
summary(reg.out)

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


## ------------------------------------------------------------------------
CV.out <- crossvalidate.bigKRLS(y = as.matrix(mtcars$mpg), seed = 123, Kfolds = 4,
                   X = as.matrix(mtcars[,-1]), Ncores = 1)
cor(CV.out$fold_3$tested$predicted, CV.out$fold_3$tested$ytest)

## ---- eval = FALSE-------------------------------------------------------
#  summary(CV.out$fold_1$trained) # not run

## ------------------------------------------------------------------------
CV.out$MSE_is
CV.out$MSE_oos
CV.out$R2_oos
CV.out$R2AME_oos

## ---- eval=F-------------------------------------------------------------
#  shiny.bigKRLS(reg.out)         # not run

## ---- eval=F-------------------------------------------------------------
#  shiny.bigKRLS(reg.out, export = T)         # not run

## ------------------------------------------------------------------------
Xnew <- mtcars[,-1]
Xnew$hp <- 200
forecast <- predict(reg.out, as.matrix(Xnew))
mean(forecast$predicted < mtcars$mpg)

## ---- eval=F-------------------------------------------------------------
#  out <- bigKRLS(y, X, model_subfolder_name = "my_results") # not run
#  save.bigKRLS(out, "my_results") # not run

## ---- eval=F-------------------------------------------------------------
#  load.bigKRLS("my_results") # not run

## ------------------------------------------------------------------------
Z <- big.matrix(nrow=5, ncol=5, init=1)
Z

## ------------------------------------------------------------------------
Z[]

