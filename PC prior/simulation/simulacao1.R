# Stan model
model_stan = rstan::stan_model(file = 'rw.stan')

# Settings
m = 2
err = hpd = prob.cov = rep(0, m)
T = 1500 #c(50, 500, 1500)
xi = c(1e-6, 0.05, 0.1, 1.0)
lambda = c(1, 5, 15)

# Initial Data frame
Data = matrix(nrow = 6, ncol = 4)
a1 = 0
r = c = 1

# run simulation
for( x in xi ){
  if( x == xi[1] ) time = Sys.time()
  for( l in lambda ){
    for( i in 1:m ){
      
      # Model
      a = rep(0, T)
      a[ 1 ] = a1
      for( t in 2:T ) a[ t ] = a[ t-1 ] + sqrt( x ) * rnorm( 1 ) 
      
      # RStan Sampling
      draws = rstan::sampling(model_stan, 
                              data = list(T = length( a ), 
                                          a = as.numeric( a ),
                                          lambda = l),
                              chains = 4,
                              warmup = 5e2,
                              iter = 5e2 + 250,
                              cores = 4
      )
      draw = rstan::extract( draws )
      draws_xi = as.numeric( draw$xi )
      draws_xi = coda::mcmc( data = draws_xi)
      # evaluation metrics
      xi_hat = mean( draws_xi )
      err[ i ] = xi_hat - x
      #piv.hpd = coda::HPDinterval( draws_xi )[1, ]
      #amp.hpd[ i ] = piv.hpd[ 2 ] - piv.hpd[ 1 ]
      #if( (piv.hpd[ 1 ] < xi_hat) && (piv.hpd[ 2 ] > xi_hat) ) prob.cov[ i ] = 1
      
    }
    vies = mean( err )
    reqm = sqrt( mean( err**2 ) )
    Data[ r, c ] = vies
    Data[ r + 1, c ] = reqm
    r = r + 2
  }
  c = c + 1
  r = 1
  if( x == xi[4] ) time = Sys.time() - time
}

time
row.names( Data ) = c( rep('lambda1', 2),
                       rep('lambda2', 2),
                       rep('lambda3', 2) )
colnames( Data ) = xi  
Data
