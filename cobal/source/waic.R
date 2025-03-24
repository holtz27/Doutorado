############### waic
waic = function(data, mu_t, h){
  n = length(data)
  M = matrix(nrow = dim(h)[1], ncol = n)
  for(j in 1:n){
    M[, j] = dnorm(data[j], 
                   mean = mu_t[,j], 
                   sd = exp(0.5*h[,j]), 
                   log = TRUE)
  }
  return(loo::waic(M))
}
