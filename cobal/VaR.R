VaR = function( h_T, a_T, theta_hmc, p = 0.95, k=1e3 ){
  
  theta = theta_hmc
  
  M = dim( theta )[ 2 ]
  # h_new = mu + phi * ( h_T - mu ) + sigma_h * eps_h
  eps_h = rnorm( M, sd = sqrt(theta[4, ]) )
  h_new = theta[2, ] + theta[3, ] * ( h_T - theta[2, ] ) + eps_h
  # a_new = a_T + eps_a
  eps_a = rnorm( M, sd = sqrt(theta[5, ]) )
  a_new = a_T + eps_a
  
  y_new = matrix( nrow = k, ncol = M)
  for(j in 1:k){
    eps_y = rnorm( M, sd = exp( 0.5 * h_new ) )
    y_new[j, ] = theta[1, ] + a_new * exp( h_new ) + eps_y 
  }
  
  VaR_hat = apply( y_new, MARGIN = 2, quantile, probs = c(1 - p) )
  
  return( mean( VaR_hat ) )
}
