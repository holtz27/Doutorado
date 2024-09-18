pred.rw = function( a_T, xi_hmc, a.test, kernel = 'epanechnikov', seed = NULL ){
  
  if( is.null(seed) ) seed = sample(1:1e6, 1)
  set.seed(seed)
  xi = xi_hmc
  M = length( xi )
  
  densidade_a.test = rep( 0, length( a.test ) )

  for( h in 1:length( a.test ) ){
    # Sampling from predictive density
    a.new = rnorm(n = M, mean = a_T, sd = sqrt(xi))
    # Calcula a densidade usando a função density do R
    densidade = density(a.new, kernel = kernel)
    # Aproxima a densidade no ponto yobs
    densidade_a.test[ h ] = approx(densidade$x, densidade$y, xout = a.test[ h ])$y
    a_T = a.new
  }

  return( -sum( log(densidade_a.test ) ) ) 
}
