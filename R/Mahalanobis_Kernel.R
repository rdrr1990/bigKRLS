###################################################
# Auxilliary functions for the Mahalanobis kernel #
###################################################
mahalanobis_kernel <- function(Z, sigma){
  Sinv <- solve(cov(Z))
  K <- matrix(NA, nrow=nrow(Z), ncol=nrow(Z))
  for(i in 1:nrow(Z)){
    for(j in 1:i){
      d <- exp(-((Z[i,] - Z[j,]) %*% Sinv %*% (Z[i,] - Z[j,]))/sigma)
      K[i,j] <- d
      K[j,i] <- d
    }
  }
  return(K)
}

mahalanobis_temp <- function(Z_new, Z_old, sigma){
  Sinv <- solve(cov(rbind(Z_new, Z_old)))
  K <- matrix(NA, nrow=nrow(Z_new), ncol=nrow(Z_old))
  
  for(i in 1:nrow(Z_new)){
    for(j in 1:nrow(Z_old)){
      d <- exp(-((Z_new[i,] - Z_old[j,]) %*% Sinv %*% (Z_new[i,] - Z_old[j,]))/sigma)
      K[i,j] <- d
    }
  }
  return(K)
}
