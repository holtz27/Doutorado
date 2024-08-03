forecast.vol = function(h_T, mu, phi, sigma, k=1){
  
  M = length(mu)
  h_plus_1 = matrix(nrow = M, ncol = k + 1)
  h_plus_1[ ,1] = h_T
  h_lmts = matrix(nrow = 2, ncol = k)
  h_new = NULL
  
  for( j in 2:(k+1) ){
    for( m in 1:M ){
      h_plus_1[m, j] = rnorm(1, 
                             mean = mu[m] + phi[m] * (h_plus_1[m, j-1] - mu[m]), 
                             sd = sigma[m]
      )
    }
    
    h_lmts[, j-1] = t( quantile(h_plus_1[, j], probs = c(0.025, 0.975) ) )
    # Adiciona a previsão ao vetor h
    h_new = rbind( h_new, mean(h_plus_1[, j]) )
  }
  return( list(h_new = h_new, h_lmts = h_lmts) )
}
