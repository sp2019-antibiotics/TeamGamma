#Ecoff
ecoff <- function(est, quant){
  ecoff <- vector("numeric", length = length(est$Pi))
  for (i in 1:length(est$Pi)) {
    if (est$Pi[i] > 0.3) {
      ecoff[i] <- qgamma(quant, shape = est$Alpha[i], rate = est$Beta[i])
    }
  }
  return(max(ecoff))
}