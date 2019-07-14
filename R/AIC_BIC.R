#AIC
AIC_Gamma <- function(loglik, G){
  aic <- -2*loglik + 2*(3*G - 1)
  return(aic)
}

#BIC
BIC_Gamma <- function(loglik, G, n = sum(Data)){
  bic <- -2*loglik + 2*(3*G - 1) * log(n)
  return(bic)
}
