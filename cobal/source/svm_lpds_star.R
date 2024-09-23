svm_lpds_star = function( h_T, theta_hmc, yobs, kernel = 'epanechnikov', seed = NULL ){
  
  if( is.null(seed) ) seed = sample(1:1e6, 1)
  set.seed(seed)
  
  # theta = (b, mu, phi, sigma, a1)
  theta = theta_hmc
  M = dim( theta )[ 2 ]
  
  log.pred = 0
  
  # h_new = mu + phi * ( h_T - mu ) + sigma_h * eps_h
  eps_h = rnorm( M, sd = sqrt(theta[4, ]) )
  h_new = theta[2, ] + theta[3, ] * ( h_T - theta[2, ] ) + eps_h
  # y_new = mu + a1 * exp{ h_T } + eps_y
  y_new = matrix( nrow = 1, ncol = M )
  eps_y = rnorm( M, sd = exp( 0.5 * h_new ) )
  y_new = theta[1, ] + theta[5, ] * exp( h_new ) + eps_y
  
  # Calcula a densidade usando a função density do R
  densidade = density(y_new, kernel = kernel)
  # Aproxima a densidade no ponto yobs
  densidade_yobs = approx(densidade$x, densidade$y, xout = yobs)$y
  
  x = log( densidade_yobs )
  if( !is.na(x) && !is.nan(x) ) log.pred = x
  
  return( log.pred ) 
}
