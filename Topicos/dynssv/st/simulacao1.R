# function
rtnorm = function(n){
  u = runif(n)
  return( qnorm( 0.5 * (u + 1) ) )
}
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
localenvir = rstudioapi::getActiveDocumentContext()
path = localenvir$path
dir = dirname(path)
path = paste0( dir, '/dynssv_st.stan')
model_stan = rstan::stan_model(file = path)

set.seed(164872)
T = 1e3
# log-volatility
mu_h = 0
phi_h = 0.99
sigma_h = 0.1
# dynskew
mu_a = -0.5
phi_a = 0.99
sigma_a = 0.1
v = 20

theta_vdd = matrix(c(mu_h, phi_h, sigma_h, 
                     mu_a, phi_a, sigma_a, 
                     v), ncol = 1)
summary = list()
prob = matrix(0, nrow = 7, ncol = 1)

M = 2
err = matrix(0, nrow = 7, ncol = M)
vies = reqm = matrix(0, nrow = 1, ncol = 3)

for(it in 1:M){
  if( it == 1 ) time = Sys.time()
  #Data
  a = delta = omega = rep(0 , T)
  a[1] = mu_a + sigma_a / sqrt( 1 - phi_a * phi_a ) * rnorm(1)
  for(t in 2:T) a[t] = mu_a + phi_a * (a[t-1] - mu_a) + sigma_a * rnorm(1)
  delta = a / sqrt(1 + a*a)
  k1 = sqrt(0.5 * v) * gamma(0.5*(v-1)) / gamma(0.5 * v)
  k2 = v / (v-2);
  omega = 1 / sqrt( k2 - 2 * ( delta * k1 )^2 / pi )
  mean_st = - sqrt(2 / pi) * delta * omega * k1
  W = rtnorm(T)
  U = rgamma(T, shape = 0.5 * v, rate = 0.5* v)
  y = h = rep(0, T)
  h[ 1 ] = mu_h + sigma_h / sqrt( 1 - phi_h * phi_h ) * rnorm(1)
  for(t in 2:T) h[ t ] = mu_h + phi_h * (h[ t-1 ] - mu_h) + sigma_h * rnorm(1)
  mu_t = mean_st + omega * delta * W * exp(0.5 * h) / sqrt( U )
  sigma_t = omega * sqrt(1 - ( delta )^2) * exp(0.5 * h) / sqrt( U )
  y = mu_t + sigma_t * rnorm(T)
  
  ### Sampling
  draws = rstan::sampling(model_stan, 
                          data = list(T = length( y ), 
                                      y = as.numeric( y ),
                                      lambda1 = 1.55,
                                      lambda2 = 15
                          ),
                          chains = 1,
                          warmup = 1e2,
                          iter = 1e2 + 2e1,
                          cores = 1
  )
  x = rstan::extract( draws )
  draws_mu_h = x$mu_h
  draws_phi_h = x$phi_h
  draws_s_h = x$s_h
  draws_mu_a = x$mu_a
  draws_phi_a = x$phi_a
  draws_s_a = x$s_a
  draws_v = x$v
  theta = matrix( draws_mu_h, nrow = 1 )
  theta = rbind( theta, draws_phi_h )
  theta = rbind( theta, draws_s_h )
  theta = rbind( theta, draws_mu_a)
  theta = rbind( theta, draws_phi_a )
  theta = rbind( theta, draws_s_a )
  theta = rbind( theta, draws_v )
  
  summary[[ it ]] = num_analisys(draws = theta, 
                                 names = c('mu_h', 'phi_h', 's_h', 
                                           'mu_a', 'phi_a', 's_a', 
                                           'v'),
                                 digits = 4,
                                 hdp = TRUE
  )
  # error
  err[ ,it ] = matrix(summary[[ it ]][, 1], ncol = 1) - theta_vdd
  # prob
  x1 = as.numeric( summary[[ it ]][, 3] < theta_vdd )
  x2 = as.numeric( summary[[ it ]][, 4] < theta_vdd )
  prob.piv = matrix( round( 0.5 * (x1 + x2), 0 ), ncol = 1)
  prob = prob + prob.piv
  
  if( it == M ) time = Sys.time() - time
}

time
vies = matrix( apply( err, MARGIN = 1, mean ), ncol = 1 )
reqm = matrix( apply( err^2, MARGIN = 1, mean ), ncol = 1 )
prob = prob / M
data = cbind( vies, reqm, prob )
data = data.frame( data )
row.names( data ) = c('mu_h', 'phi_h', 's_h', 
                      'mu_a', 'phi_a', 's_a', 
                      'v')
colnames( data) = c('vies', 'reqm', 'prob.cob')
round( data, 4 )
