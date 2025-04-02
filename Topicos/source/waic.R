waic = function(data, mu_t, sigma_t){
  n = length(data)
  M = matrix(nrow = dim(sigma_t)[1], ncol = n)
  for(j in 1:n){
    M[, j] = dnorm(data[j], 
                   mean = mu_t[,j], 
                   sd = sigma_t[,j], 
                   log = TRUE)
  }
  return(loo::waic(M))
}

