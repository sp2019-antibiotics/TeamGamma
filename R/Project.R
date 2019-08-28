globalVariables(c("Data", "Antimicrobial", "Bacterium"))

#' EM Algorithm for Gamma Mixture Models
#'
#' @param ZoneData data frame of the zone diameter data from web scraping.
#' @param Anti the chosen Antimicrobial for the analysis.
#' @param Bac the chosen Bacterium for the analysis.
#' @param k the chosen number of clusters.
#' @param a the scale parameter of the inverse gamma distribution used for the penalty. Default is 2.
#' @param b the rate parameter of the inverse gamma distribution used for the penalty. Default is 1.
#' @param epsilon the stopping criterion for the algorithm. Default is 0.00001.
#' @param quant the quantile used for the ECOFF. Default is 0.01.
#'
#' @description
#' EM_Gamma is used to perform an EM-Algorithm for binned data, using Gamma Mixture Models. It computes an ECOFF value as well as an AIC-value and BIC-value for model comparison.
#'
#' @details
#' For initialization clustering with the quantiles method is used. The data is split into clusters with equal amount of data values. Afterwards the parameters pi, alpha and beta are initialized using moment estimation.
#' The EM-Algorithm consists of an E-Step and a M-Step. The algorithm is an iterative process that updates the parameters and increases the log-likelihood at each iteration, until a stopping criterion epsilon is reached.
#' In the E-Step the current parameter values are used to compute p_j^(i)
#'            p_j^(i)= sum_(k=1)^K pi_k^(i) * integral_(a_j)^(b_j) f_k(x|vartheta_k^(i))dx
#' and p_jk^(i)
#'            p_jk^(i) = [pi_k^(i) * integral_(a_j)^(b_j) f_k(x|vartheta_k^(i))dx] / p_j^(i),
#' with f_k(x) the density of the k-th gamma mixture component.
#' In each M-Step the expectation of the complete log-likelihood at the i-th step has to be maximized. The maximization with respect to pi yields 
#'            pi_k^(i+1)= [sum_(j=1)^J n_j*p_jk^(i)] / [tilde(N)]
#' Using optim as the general purpose optimizer, the estimators of alpha and beta are obtained.
#' This algorithm is repeated until
#'            |l^(i+1)-l^(i)| < epsilon,
#' with l^(i) being the complete data log-likelihood at the i-th step.
#' Finally, AIC, BIC and the ECOFF-value get computed using the calculated parameters alpha, beta and pi.
#'
#' @return A list with the components:
#'
#' \itemize{
#'  \item{Pi (the estimated mixture weights)}
#'  \item{Alpha (the estimated scale parameters for the gamma distribution)}
#'  \item{Beta (the estimated rate parameters for the gamma distribution)}
#'  \item{AIC (the Akaike information criterion)}
#'  \item{BIC (the Bayesian information criterion)}
#'  \item{Ecoff (the value for the Ecoff)}
#' }
#'
#' @export
#'
#' @examples
#' data("ZD", package = "EUCASTData") #loading the Zone Diameter-Data from a CSV-file
#'
#' # example 1
#' # plotting the frequencies of the Zone Diameters for visual checking
#' # selecting the cells with Zone Diameters of a specific combination of Antimicrobial and Bacterium
#' ZDs <- subset(ZD, Antimicrobial == "Ampicillin" & Bacterium == "Escherichia coli",
#'               grepl("^Z", colnames(ZD)))
#'
#' example <- data.frame(ZD = as.integer(gsub("^Z", "", colnames(ZDs))), # creating a dataframe
#'                       Freq = unname(unlist(ZDs)))
#'
#' plot(Freq ~ ZD, data = example, type = "h") #plotting the frequencies of the Zone Diameters
#'
#' #as the visual checking of the plot suggests a distribution with 1 component k will be chosen as 1
#'
#' # calling the EM_Gamma-function with a specific combination of
#' # Antimicrobial, Bacterium, k and epsilon as default
#' est <- EM_Gamma(ZoneData = ZD, Anti = "Ampicillin", Bac = "Escherichia coli",
#'                 k = 1, epsilon = 0.001)
#' est # showing the parameters of the distribution, AIC, BIC and Ecoff
#'
#' # geting the theoretical data values for the histogramm
#' Data <- unname(unlist(ZDs))[-1]
#' dia <- as.integer(gsub("^Z", "", colnames(ZDs)))[-1]
#'
#' y <- vector("numeric")
#'
#' for (i in 1:length(Data)) {
#'   if (Data[i] == 1) {
#'     y <- c(y, dia[i])
#'   } else{
#'     y <- c(y, seq(dia[i]-0.5, dia[i]+0.5, length.out = Data[i]))
#'   }
#' }
#'
#' # plotting the histogramm of the data
#' hist(y, breaks=seq(0,50,length=50), freq = FALSE, main="Ampicillin & Escherichia coli",
#'      xlab="Zone Diameter", ylab="Density")
#'
#' # plotting the fitted density
#' lines(seq(0, 50, length = 1000),
#'       est$Pi[1]*dgamma(seq(0, 50, length = 1000), shape = est$Alpha[1], rate = est$Beta[1]),
#'       col = 2, lwd = 2)
#'
#' # example 2
#' # plotting the frequencies of the Zone Diameters for visual checking
#' # selecting the cells with Zone Diameters of a specific combination of Antimicrobial and Bacterium
#' ZDs <- subset(ZD, Antimicrobial == "Piperacillin" & Bacterium == "Escherichia coli",
#'               grepl("^Z", colnames(ZD)))
#'
#' example <- data.frame(ZD = as.integer(gsub("^Z", "", colnames(ZDs))), # creating a dataframe
#'                       Freq = unname(unlist(ZDs)))
#'
#' plot(Freq ~ ZD, data = example, type = "h") #plotting the frequencies of the Zone Diameters
#'
#' # as the visual checking of the plot suggests a distribution with 2 components k will be chosen as 2
#' # calling the EM_Gamma-function with a specific combination of
#' # Antimicrobial, Bacterium, k and epsilon as default
#' est <- EM_Gamma(ZoneData = ZD, Anti = "Piperacillin", Bac = "Escherichia coli",
#'                 k = 2, epsilon = 0.001)
#' est #showing the parameters of the distributions, AIC, BIC and Ecoff
#'
#' # geting the theoretical data values for the histogramm
#' Data <- unname(unlist(ZDs))[-1]
#' dia <- as.integer(gsub("^Z", "", colnames(ZDs)))[-1]
#'
#' y <- vector("numeric")
#'
#' for (i in 1:length(Data)) {
#'   if (Data[i] == 1) {
#'     y <- c(y, dia[i])
#'   } else{
#'     y <- c(y, seq(dia[i]-0.5, dia[i]+0.5, length.out = Data[i]))
#'   }
#' }
#'
#' # plotting the histogramm of the data
#' hist(y, breaks=seq(0,50,length=50), freq = FALSE, main="Piperacillin & Escherichia coli",
#'      xlab="Zone Diameter", ylab="Density")
#'
#' # plotting the fitted density
#' lines(seq(0, 50, length = 1000),
#'       est$Pi[1]*dgamma(seq(0, 50, length = 1000), shape = est$Alpha[1], rate = est$Beta[1])
#'       +est$Pi[2]*dgamma(seq(0, 50, length = 1000), shape = est$Alpha[2], rate = est$Beta[2]),
#'       col = 2, lwd = 2)
#'
#' # example 3
#' # plotting the frequencies of the Zone Diameters for visual checking
#' # selecting the cells with Zone Diameters of a specific combination of Antimicrobial and Bacterium
#' ZDs <- subset(ZD, Antimicrobial == "Mecillinam" & Bacterium == "Escherichia coli",
#'               grepl("^Z", colnames(ZD)))
#'
#' example <- data.frame(ZD = as.integer(gsub("^Z", "", colnames(ZDs))), # creating a dataframe
#'                       Freq = unname(unlist(ZDs)))
#'
#' plot(Freq ~ ZD, data = example, type = "h") #plotting the frequencies of the Zone Diameters
#'
#' # after the visual checking of the plot it cannot be decided clearly how many components there are,
#' # so several k will be tried out
#' # 1 component
#' est1 <- EM_Gamma(ZoneData = ZD, Anti = "Mecillinam", Bac = "Escherichia coli",
#'                  k = 1, epsilon = 0.001)
#' # 2 components
#' est2 <- EM_Gamma(ZoneData = ZD, Anti = "Mecillinam", Bac = "Escherichia coli",
#'                  k = 2, epsilon = 0.001)
#' # 3 components
#' est3 <- EM_Gamma(ZoneData = ZD, Anti = "Mecillinam", Bac = "Escherichia coli",
#'                  k = 3, epsilon = 0.001)
#'
#' # comparing the BIC and the AIC to figure out the best solution
#' est1
#' est2
#' est3 # according to BIC and AIC the distribution with 3 components is the best solution
#'
#' # geting the theoretical data values for the histogramm
#' Data <- unname(unlist(ZDs))[-1]
#' dia <- as.integer(gsub("^Z", "", colnames(ZDs)))[-1]
#'
#' y <- vector("numeric")
#'
#' for (i in 1:length(Data)) {
#'   if (Data[i] == 1) {
#'     y <- c(y, dia[i])
#'   } else{
#'     y <- c(y, seq(dia[i]-0.5, dia[i]+0.5, length.out = Data[i]))
#'   }
#' }
#'
#' # plotting the histogramm of the data
#' hist(y, breaks=seq(0,50,length=50), freq = FALSE, main="Mecillinam & Escherichia coli",
#'      xlab="Zone Diameter", ylab="Density")
#'
#' # plotting the fitted densities
#' legend('topright', legend=c("k = 1", "k = 2", "k = 3"), lty=1, col=c('black', 'red', 'green'),
#'        bty='n', lwd = 2)
#' lines(seq(0, 50, length = 1000),
#'       est1$Pi[1]*dgamma(seq(0, 50, length = 1000), shape = est1$Alpha[1], rate = est1$Beta[1]),
#'       lwd = 2)
#'
#' lines(seq(0, 50, length = 1000),
#'       est2$Pi[1]*dgamma(seq(0, 50, length = 1000), shape = est2$Alpha[1], rate = est2$Beta[1])
#'       +est2$Pi[2]*dgamma(seq(0, 50, length = 1000), shape = est2$Alpha[2], rate = est2$Beta[2]),
#'       type = "l", col="red", lwd = 2)
#'
#' lines(seq(0, 50, length = 1000),
#'       est3$Pi[1]*dgamma(seq(0, 50, length = 1000), shape = est3$Alpha[1], rate = est3$Beta[1])
#'       +est3$Pi[2]*dgamma(seq(0, 50, length = 1000), shape = est3$Alpha[2], rate = est3$Beta[2])
#'       +est3$Pi[3]*dgamma(seq(0, 50, length = 1000), shape = est3$Alpha[3], rate = est3$Beta[3]),
#'       type = "l", col="green", lwd = 2)

EM_Gamma <- function(ZoneData, Anti, Bac, k, a = 2, b = 1, epsilon = 0.00001, quant = 0.01){
  if(is.data.frame(ZoneData) == FALSE) stop("ZoneData must be the Dataframe from EUCAST for the ZoneDiameter.")
  if (is.character(Anti) == FALSE | is.character(Bac) == FALSE) stop("Anti or Bac must be a character.")
  if (k < 1) stop("k should be at least 1")
  DSub <- subset(ZoneData, Antimicrobial == Anti & Bacterium == Bac,
                 grepl("^Z", colnames(ZoneData)))

  Data <- unname(unlist(DSub))[-1]
  if (sum(is.na(Data)) > 0) stop("at least one bin is NA")
  dia <- as.integer(gsub("^Z", "", colnames(DSub)))[-1]
  if (sum(Data) < 2*k) stop(paste("not enough observations for", k, "clusters"))
  ob <- sum(Data)

  Cluster <- ClustQuant(Data = Data, k = k, dia = dia)

  est <- Initialization(Cluster = Cluster, ob = ob)

  e <- 1
  while(e > epsilon){
    pjk <- Pjk(dia, est$Pi, est$Alpha, est$Beta)

    newP <- newPi(ob, Data, pjk)

    newAlpha <- vector("numeric", length = k)
    newBeta <- vector("numeric", length = k)

    for (i in 1:k) {
      par <- optim(par = c(est$Alpha[i], est$Beta[i]), fn = Q, a = a, b = b, dia = dia, pjk = pjk, Data = Data, k = i,
                   control = list(fnscale = -1), method = "L-BFGS-B", lower = rep(.Machine$double.xmin, 2))$par
      newAlpha[i] <- par[1]
      newBeta[i] <- par[2]
    }
    newEst <- list(Pi = newP, Alpha = newAlpha, Beta = newBeta)


    e <- abs(l(newEst$Pi, newEst$Alpha, newEst$Beta, a, b, Data, dia)-l(est$Pi, est$Alpha, est$Beta, a, b, Data, dia))

    est <- newEst
  }

  f <- vector("numeric", length = length(dia))
  for (g in 1:length(est$Pi)) {
    f <- f + est$Pi[g]*((pgamma(dia+0.5, shape = est$Alpha[g], rate = est$Beta[g]) - pgamma(dia-0.5, shape = est$Alpha[g], rate = est$Beta[g])))
  }

  f[f == 0] <- .Machine$double.xmin

  loglik <- sum(Data*log(f))

  G <- 3*k

  AIC <- AIC_Gamma(loglik = loglik, G = G)

  est[["AIC"]] <- AIC

  BIC <- BIC_Gamma(loglik = loglik, G = G, n = sum(Data))

  est[["BIC"]] <- BIC

  Ecoff <- ecoff(est = est, quant = quant)

  est[["Ecoff"]] <- Ecoff

  return(est)
}
