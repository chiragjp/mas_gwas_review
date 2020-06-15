###compute max AUC

#code adopted from https://github.com/kn3in/genRoc

compute_h2_liab <- function(k, h_2, p, z) {
  h2l=((h_2*k*(1-k))/(z^2))*(k*(1-k)/(p*(1-p)))
 return(h2l)
}

final_results <-  function(k, h_2, p) {
  
  T0 <- qnorm(1 - k)
  z  <- dnorm(T0)
  i  <- z / k
  
  h_2_l=compute_h2_liab(k, h_2, p ,z)
  
  print(h_2_l)

  T1 <- (2 * T0 - i * h_2_l) / sqrt(4 - h_2_l^2 * i * (i - T0))
  lambda_s <- (1 - pnorm(T1)) / k
  
  v  <- -i * (k / (1-k))
  
  dn <- function(rho) {
    (i - v) * h_2_l * rho / sqrt(h_2_l * rho * (1 - h_2_l * rho * i * (i - T0) + 1 - h_2_l * rho * v * (v - T0)))
  }
  
  auc <- function(rho) pnorm(dn(rho))
  
  return(auc(1))
  
}
