# Stan model
model_stan1 = rstan::stan_model(file = 'classic_tvp_svm.stan')
model_stan2 = rstan::stan_model(file = 'pcp_tvp_svm.stan')

# Data Settings
b = 0.1
mu = -9
phi = 0.985
s2_h = 0.025
T = 1e3

# Simulation Settings
m = 50
c.err = pcp.err = rep(0, m)
xi = c(0, 0.05, 0.1, 1.0)
U = invgamma::qinvgamma(0.975, shape = 4.5, rate = 0.065)
lambda = -log(0.025) / sqrt( U )

# Initial Data frame
Data = matrix(nrow = 4, ncol = 4)
c = 1

# run simulation
for( x in xi ){
  if( x == xi[1] ) time = Sys.time()
  for( i in 1:m ){
    
    # Model
    y = h = a = matrix(0, nrow = T, ncol = 1)
    a[1] = -0.05 
    h[1] = mu + sqrt( s2_h / (1 - phi * phi) ) * rnorm( 1 )
    y[1] = b + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )
    for( t in 2:T ){
      a[t] = a[t-1] + sqrt( xi ) * rnorm( 1 )
      h[t] = mu + phi * (h[t-1] - mu) + sqrt( s2_h ) * rnorm( 1 )
      y[t] = b + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
    }
    ##########################
    # Classic RStan Sampling #
    ##########################
    draws = rstan::sampling(model_stan1, 
                            data = list(T = length( y ), 
                                        y = as.numeric( y )
                            ),
                            chains = 2,
                            warmup = 1e3,
                            iter = 1e3 + 500,
                            cores = 2
    )
    draw = rstan::extract( draws )
    draws_xi = as.numeric( draw$s2_a )
    # evaluation metrics
    xi_hat = mean( draws_xi )
    if( x == 0){
      c.err[ i ] = xi_hat - x
    }else{
      c.err[ i ] = (xi_hat - x) / x
    }
    
    ######################
    # PCP RStan Sampling #
    ######################
    draws = rstan::sampling(model_stan2, 
                            data = list(T = length( y ), 
                                        y = as.numeric( y ),
                                        lambda = lambda 
                            ),
                            chains = 2,
                            warmup = 1e3,
                            iter = 1e3 + 500,
                            cores = 2
    )
    draw = rstan::extract( draws )
    draws_xi = as.numeric( draw$s2_a )
    # evaluation metrics
    xi_hat = mean( draws_xi )
    if( x == 0){
      pcp.err[ i ] = xi_hat - x
    }else{
      pcp.err[ i ] = (xi_hat - x) / x
    }
    
  }
  c.vies = mean( c.err )
  pcp.vies = mean( pcp.err )
  c.reqm = sqrt( mean( c.err**2 ) )
  pcp.reqm = sqrt( mean( pcp.err**2 ) )
  
  Data[ 1, c ] = c.vies
  Data[ 2, c ] = c.reqm
  Data[ 3, c ] = pcp.vies
  Data[ 4, c ] = pcp.reqm
  c = c + 1
  if( x == xi[4] ) time = Sys.time() - time
}

time
row.names( Data ) = c('c.vies', 'c.reqm', 'pcp.vies', 'pcp.reqm')
colnames( Data ) = xi  
round( Data, 4 )

#save(Data, file = 'Data.RData')
