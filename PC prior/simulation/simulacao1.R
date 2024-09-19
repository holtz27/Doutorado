source( 'https://raw.githubusercontent.com/holtz27/Doutorado/main/PC%20prior/simulation/pred.rw.R' )
# Stan model
model_stan = rstan::stan_model(file = 'rw.stan')

# Settings
h = 10
m = 2
err = pred = rep(0, m)
T = 1500 + h #c(50, 500, 1500)
xi = c(1e-6, 0.05, 0.1, 1.0)
lambda = c(1, 5, 15)

# Initial Data frame
Pred = matrix(nrow = 3, ncol = 4)
Data = matrix(nrow = 6, ncol = 4)
a1 = 0
r = c = p = 1

# run simulation
for( x in xi ){
  if( x == xi[1] ) time = Sys.time()
  for( l in lambda ){
    for( i in 1:m ){
      
      # Model
      a = rep(0, T)
      a[ 1 ] = a1
      for( t in 2:T ) a[ t ] = a[ t-1 ] + sqrt( x ) * rnorm( 1 ) 
      a.train = a[1:(T-h)]
      a.test = a[(T-h+1):T]
      
      # RStan Sampling
      draws = rstan::sampling(model_stan, 
                              data = list(T = length( a.train ), 
                                          a = as.numeric( a.train ),
                                          lambda = l),
                              chains = 4,
                              warmup = 5e2,
                              iter = 5e2 + 250,
                              cores = 4
      )
      draw = rstan::extract( draws )
      draws_xi = as.numeric( draw$xi )
      # evaluation metrics
      xi_hat = mean( draws_xi )
      err[ i ] = (xi_hat - x) / x
      ### out of sample metric
      pred[ i ] = pred.rw(a_T = a[T-h], xi_hmc = draws_xi, a.test = a.test)
    }
    vies = mean( err )
    reqm = sqrt( mean( err**2 ) )
    Data[ r, c ] = vies
    Data[ r + 1, c ] = reqm
    Pred[ p, c ] = mean( pred )
    p = p + 1
    r = r + 2
  }
  c = c + 1
  r = 1
  p = 1
  if( x == xi[4] ) time = Sys.time() - time
}

time
row.names( Data ) = c( rep('lambda1', 2),
                       rep('lambda2', 2),
                       rep('lambda3', 2) )
row.names( Pred ) = c( 'lambda1',
                       'lambda2',
                       'lambda3' )
colnames( Data ) = xi  
colnames( Pred ) = xi  
Data
Pred
