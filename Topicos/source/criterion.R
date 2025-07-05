waic = function(data, mu_t, sigma_t){
  n = length(data)
  M = matrix(nrow = dim(sigma_t)[1], ncol = n)
  for(j in 1:n){
    M[,j] = dnorm(data[j], 
                  mean = mu_t[,j], 
                  sd = sigma_t[,j], 
                  log = TRUE)
  }
  return(loo::waic(M))
}
loo = function(data, mu_t, sigma_t){
  n = length(data)
  M = matrix(nrow = dim(sigma_t)[1], ncol = n)
  for(j in 1:n){
    M[,j] = dnorm(data[j], 
                  mean = mu_t[,j], 
                  sd = sigma_t[,j], 
                  log = TRUE)
  }
  return(loo::loo(M, cores=5))
}
dic = function(data, mu_t, sigma_t){
  n_iter = dim(mu_t)[1]
  n_obs  = length(data)
  
  # 1. Log-verossimilhança para cada iteração e tempo
  loglik_matrix = matrix(NA, nrow = n_iter, ncol = n_obs)
  for(j in 1:n_obs){
    loglik_matrix[,j] = dnorm(data[j],
                              mean = mu_t[,j],
                              sd   = sigma_t[,j],
                              log  = TRUE)
  }
  
  # 2. Deviance média: \bar{D}
  D_bar = -2 * mean(rowSums(loglik_matrix))
  
  # 3. Parâmetros médios: \hat{\mu}_t, \hat{\sigma}_t
  mu_hat     = colMeans(mu_t)
  sigma_hat  = colMeans(sigma_t)
  
  # 4. Deviance na média posterior: D(\hat{\theta})
  D_hat = -2 * sum(dnorm(data,
                         mean = mu_hat,
                         sd   = sigma_hat,
                         log  = TRUE))
  # 5. Penalização de complexidade efetiva
  p_D = D_bar - D_hat
  # 6. Cálculo final do DIC
  dic = D_bar + p_D
  
  return(dic)
}
