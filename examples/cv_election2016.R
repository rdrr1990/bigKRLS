library(bigKRLS)

# run 
# source("https://raw.githubusercontent.com/rdrr1990/bigKRLS/master/examples/kcv_election2016.R")
# cv() with starting seed and desired number of reps

replication_data <- read.table("https://raw.githubusercontent.com/rdrr1990/bigKRLS/master/examples/data2016GEcensus.csv", header = TRUE, sep=",")
# reading from github seems to read row numbers as the first column (called 'X')

y <- as.matrix(replication_data$GOPdelta)
X <- as.matrix(replication_data)[, 3:ncol(replication_data)]

summaries <- list()

cv <- function(startseed, Nreps){
  
  seeds <- startseed:(startseed + Nreps - 1)
  
  for(i in 1:Nreps){
    
    out <- crossvalidate.bigKRLS(y, X, ptesting = 20, 
                                 Neig = 50, # Neig assumes bigKRLS 2.1.0 or newer 
                                 seed = seeds[i])
    summaries[[i]] <- summary(out)
    overview <- summaries[[i]][["overview"]]
    
    overview <- cbind(overview, seeds[i])

    if(i == 1){
      overviews <- overview[-5,]
      colnames(overviews)[ncol(overviews)] <- "seed"
    }else{
      overviews <- rbind(overviews, overview[-5,])
    }
    
    cat("\n\nfinished", i, "\n\n")
    write.csv(overviews, 
         file = paste0("cv_seeds_", startseed, "_to_", startseed + Nreps - 1, ".csv")) 
    save(summaries, 
        file = paste0("cv_seeds_summaries_", startseed, "_to_", startseed + Nreps - 1, ".RData")) 
  }
}


