#Quantiles
ClustQuant <- function(Data, k, dia){
  #checking if the sum of all observations is divisible by 0
  if(sum(Data)%/%k != 0){
    edge <- ceiling(sum(Data)/k)
  } else{
    edge <- sum(Data)/k
  }

  s <- 0
  ind <- 1
  newData <- list()
  i <- 1
  #spliting the Data into k cluster
  while (i <= length(Data)) {
    s <- s + Data[i]
    if(s >= edge){
      if (s > edge) {
        Data[i] <- Data[i] + edge - s
        newData[[ind]] <- Data[1:i]
        newData[[ind]] <- rbind(newData[[ind]], dia[1:i])
        Data[i] <- s - edge
        if (i > 1) {
          Data <- Data[-c(1:i-1)]
          dia <- dia[-c(1:i-1)]
        }
        s <- 0
        i <- 0
        ind <- ind + 1
      }else{
        newData[[ind]] <- Data[1:i]
        newData[[ind]] <- rbind(newData[[ind]], dia[1:i])
        Data <- Data[-c(1:i)]
        dia <- dia[-c(1:i)]
        s <- 0
        i <- 0
        ind <- ind + 1
      }
    }
    if(ind == k){
      newData[[ind]] <- Data
      newData[[ind]] <- rbind(newData[[ind]], dia[1:length(dia)])
    }
    i <- i + 1
  }
  if (sum(newData[[k]][1,]) <= 1) stop(paste("not enough observations for", k, "cluster, please try a different number for k"))
  return(newData) #a list of all clusters
}

