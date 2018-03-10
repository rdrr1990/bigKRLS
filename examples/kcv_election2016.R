library(bigKRLS)

replication_data <- read.table("https://raw.githubusercontent.com/rdrr1990/bigKRLS/master/examples/data2016GEcensus.csv", header = TRUE, sep=",")
# reading from github seems to read row numbers as the first column (called 'X')

y <- as.matrix(replication_data$GOPdelta)
X <- as.matrix(replication_data)[, 3:ncol(replication_data)]

kcv <- function(startseed, Nreps){
  
  seeds <- startseed:(startseed + Nreps - 1)
  
  for(i in 1:Nreps){
    
    overview <- summary(crossvalidate.bigKRLS(y, X, Kfolds = 5, 
                                              Neig = 40, # Neig assumes bigKRLS 2.1.0 or newer 
                                              seed = seeds[i]))[["overview"]]
    
    overview <- cbind(overview, seeds[i])
    
    if(i == 1){
      overviews <- overview
      colnames(overviews)[ncol(overviews)] <- "seed"
    }else{
      overviews <- rbind(overviews, overview)
    }
    cat("\n\nfinished", i, "\n\n")
    write.csv(overviews, 
         file = paste0("kcv_seeds_", startseed, "_to_", startseed + Nreps - 1, ".csv")) 
  }
}


