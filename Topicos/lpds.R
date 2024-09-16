lpds = function( h_T, a_T, theta_hmc, yobs, kernel = 'epanechnikov', seed = NULL ){
  
  if( is.null(seed) ) seed = sample(1:1e6, 1)
  set.seed(seed)
  theta = theta_hmc
  M = dim( theta )[ 2 ]
  
  densidade_yobs = rep( 0,length( yobs ) )
  
  for( h in 1:length( yobs ) ){
    # h_new = mu + phi * ( h_T - mu ) + sigma_h * eps_h
    eps_h = rnorm( M, sd = sqrt(theta[4, ]) )
    h_new = theta[2, ] + theta[3, ] * ( h_T - theta[2, ] ) + eps_h
    # a_new = a_T + eps_a
    eps_a = rnorm( M, sd = sqrt(theta[5, ]) )
    a_new = a_T + eps_a
    # y_new = mu + a_T * exp{ h_T } + eps_y
    y_new = matrix( nrow = 1, ncol = M )
    eps_y = rnorm( M, sd = exp( 0.5 * h_new ) )
    y_new = theta[1, ] + a_new * exp( h_new ) + eps_y
    
    # Calcula a densidade usando a função density do R
    densidade = density(y_new, kernel = kernel)
    # Aproxima a densidade no ponto yobs
    densidade_yobs[ h ] = approx(densidade$x, densidade$y, xout = yobs[ h ])$y
    h_T = h_new
    a_T = a_new
  }
  
  return( densidade_yobs ) 
}
