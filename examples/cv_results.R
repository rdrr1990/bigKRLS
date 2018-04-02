raw <- list()
for(i in dir("examples/cv_replicates", pattern = "csv", full.names = TRUE))
  raw[[i]] <- read.csv(i)

results <- raw[[i]]
for(i in 2:length(results))
  results <- rbind(results, raw[[i]])
colnames(results)[1] <- "stat"

r <- array(dim=c(5, 2, 100))

for(i in 1:100){
  r[,,i] <- results[(i + 5*(i - 1)):(5*i),]
}

library(dplyr)
library(tidyverse)

cbind(by(results$In.Sample, results$stat, mean), 
      by(results$Out.of.Sample, results$stat, mean))


xbar.in <- by(results$In.Sample, results$stat, mean)
xbar.out <- by(results$Out.of.Sample, results$stat, mean)

xbar.in[1] <- xbar.in[1]^.5 # switching to RMSE
xbar.in[2] <- xbar.in[2]^.5

se.in <- by(results$In.Sample, results$stat, sd) / (2485^.5)
se.out <- by(results$Out.of.Sample, results$stat, sd) / (621^.5)

uppers.in <- 1.96*se.in + xbar.in
uppers.out <- -1.96*se.out + xbar.out

uppers.in <- uppers.in[-3]
uppers.out <- uppers.out[-3]

lowers.in <- 1.96*se.in + xbar.in
lowers.out <- -1.96*se.out + xbar.out

lowers.in <- uppers.in[-3]
lowers.out <- uppers.out[-3]

cat(paste0("& ", round(xbar.in[-3], 2), " & ", round(xbar.out[-3], 2), 
           " \\tabularnewline \n",
  "& (", round(lowers.in, 2), ", ", round(uppers.in, 2), ") & ",
           "(", round(lowers.out, 2), ", ", 
  round(uppers.out, 2), ") \\tabularnewline \n"
           ))
