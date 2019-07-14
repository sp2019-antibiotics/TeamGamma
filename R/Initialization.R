#Initialization
Initialization <- function(Cluster, ob){
  p <- vector("numeric", length = length(Cluster))
  alpha <- vector("numeric", length = length(Cluster))
  beta <- vector("numeric", length = length(Cluster))
  
  for (i in 1:length(Cluster)) {
    y <- vector("numeric")
    for (j in 1:ncol(Cluster[[i]])) {
      if (Cluster[[i]][1,j] == 1) {
        y <- c(y, Cluster[[i]][2,j])
      } else{
        y <- c(y, seq(Cluster[[i]][2,j]-0.5, Cluster[[i]][2,j]+0.5, length.out = Cluster[[i]][1,j]))
      }
    }
    obc <- sum(Cluster[[i]][1,])
    
    p[i] <- obc/ob
    alpha[i] <- (mean(y)^2)/(mean(y^2-mean(y)^2))
    beta[i] <- (mean(y))/(mean(y^2-mean(y)^2))
  }
  return(list(Pi = p, Alpha = alpha, Beta = beta))
}
