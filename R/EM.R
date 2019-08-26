#' @importFrom invgamma dinvgamma
#' @importFrom stats optim pgamma qgamma
#Pj
Pj <- function(dia, p, alpha, beta){
  pj <- vector("numeric", length = length(dia))
  for (k in 1:length(p)) {
    pj <- pj + (pgamma(dia+0.5, shape = alpha[k], rate = beta[k])
                - pgamma(dia-0.5, shape = alpha[k], rate = beta[k]))*p[k]
  }
  pj[pj == 0] <- .Machine$double.xmin
  return(pj)
}

#Pjk
Pjk <- function(dia, p, alpha, beta){
  pjk <- matrix(0, nrow = length(p), ncol = length(dia))
  for (k in 1:length(p)) {
    pjk[k,] <- (p[k]*(pgamma(dia+0.5, shape = alpha[k], rate = beta[k])
                      - pgamma(dia-0.5, shape = alpha[k], rate = beta[k])))/Pj(dia = dia, p = p, alpha = alpha, beta = beta)
  }
  return(pjk)
}

#Gamma
G <- function(par, dia){
  G <- vector("numeric", length = length(dia))
  G <- pgamma(dia+0.5, shape = par[1], rate = par[2]) - pgamma(dia-0.5, shape = par[1], rate = par[2])
  for (i in 1:length(dia)) {
    if(G[i] == 0){
      G[i] <- .Machine$double.xmin
    }
  }
  return(G)
}

#new Pi
newPi <- function(ob, n, pjk){
  newPi <- vector("numeric", length = nrow(pjk))
  for (k in 1:nrow(pjk)) {
    newPi[k] <- sum(n*pjk[k,])/ob
  }
  return(newPi)
}

#for optim
Q <- function(par, a, b, dia, Data, pjk, k){
  q <- sum(Data*pjk[k,]*log(G(par, dia)))
  ig <- dinvgamma(par[1]/(par[2]^2), shape = a, rate = b)
  if (ig == 0) {
    ig <- .Machine$double.xmin
  }
  q <- sum(q) + log(ig)
}

#LH
GLo <- function(alpha, beta, dia){
  G <- matrix(0, nrow = length(alpha), ncol = length(dia))
  for (k in 1:length(alpha)) {
    G[k,] <- pgamma(dia+0.5, shape = alpha[k], rate = beta[k]) - pgamma(dia-0.5, shape = alpha[k], rate = beta[k])
  }
  for (i in 1:length(alpha)) {
    for (j in 1:length(dia)){
      if(G[i,j] == 0){
        G[i,j] <- .Machine$double.xmin
      }
    }
  }
  return(G)
}

#L
l <- function(p, alpha, beta, a, b, Data, pjk, dia){
  sum(colSums(Data*pjk*log(p))) + sum(colSums(Data*pjk*log(GLo(alpha, beta, dia)))) + sum(dinvgamma(alpha/(beta^2), shape = a, rate = b))
}

