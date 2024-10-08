# function
rtnorm = function(n){
  u = runif(n)
  return( qnorm( 0.5 * (u + 1) ) )
}
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
localenvir = rstudioapi::getActiveDocumentContext()
path = localenvir$path
dir = dirname(path)
path = paste0( dir, '/dynssv_st1.stan')
model_stan1 = rstan::stan_model(file = path)
path = paste0( dir, '/dynssv_st2.stan')
model_stan2 = rstan::stan_model(file = path)
path = paste0( dir, '/dynssv_st3.stan')
model_stan3 = rstan::stan_model(file = path)

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
                     mu_a, phi_a, log(sigma_a), 
                     v), ncol = 1)
summary1 = summary2 = summary3 = list()
prob1 = prob2 = prob3 = matrix(0, nrow = 7, ncol = 1)

M = 2
err1 = err2 = err3= matrix(0, nrow = 7, ncol = M)
vies1 = reqm1 = vies2 = reqm2 = vies3 = reqm3 = matrix(0, nrow = 1, ncol = 3)

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
  
  ### Sampling1
  draws = rstan::sampling(model_stan1, 
                          data = list(T = length( y ), 
                                      y = as.numeric( y ),
                                      lambda1 = 1.55,
                                      lambda2 = 14.34
                          ),
                          chains = 1,
                          warmup = 2e2,
                          iter = 2e2 + 2e2,
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
  
  summary1[[ it ]] = num_analisys(draws = theta, 
                                  names = c('mu_h', 'phi_h', 's_h', 
                                            'mu_a', 'phi_a', 'ls_a', 
                                            'v'),
                                  digits = 4,
                                  hdp = TRUE
                                  )
  # error
  err1[ ,it ] = matrix(summary1[[ it ]][, 1], ncol = 1) - theta_vdd
  # prob
  x1 = as.numeric( summary1[[ it ]][, 3] < theta_vdd )
  x2 = as.numeric( summary1[[ it ]][, 4] > theta_vdd )
  prob.piv = matrix( round( 0.5 * (x1 + x2), 0 ), ncol = 1)
  prob1 = prob1 + prob.piv
  
  ### Sampling2
  draws = rstan::sampling(model_stan2, 
                          data = list(T = length( y ), 
                                      y = as.numeric( y ),
                                      lambda1 = 1.55,
                                      lambda2 = 14.34
                          ),
                          chains = 1,
                          warmup = 2e2,
                          iter = 2e2 + 2e2,
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
  
  summary2[[ it ]] = num_analisys(draws = theta, 
                                  names = c('mu_h', 'phi_h', 's_h', 
                                            'mu_a', 'phi_a', 'ls_a', 
                                            'v'),
                                  digits = 4,
                                  hdp = TRUE
  )
  # error
  err2[ ,it ] = matrix(summary2[[ it ]][, 1], ncol = 1) - theta_vdd
  # prob
  x1 = as.numeric( summary2[[ it ]][, 3] < theta_vdd )
  x2 = as.numeric( summary2[[ it ]][, 4] > theta_vdd )
  prob.piv = matrix( round( 0.5 * (x1 + x2), 0 ), ncol = 1)
  prob2 = prob2 + prob.piv
  
  ### Sampling3
  draws = rstan::sampling(model_stan3, 
                          data = list(T = length( y ), 
                                      y = as.numeric( y ),
                                      lambda1 = 1.55,
                                      lambda2 = 14.34
                          ),
                          chains = 1,
                          warmup = 2e2,
                          iter = 2e2 + 2e2,
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
  
  summary3[[ it ]] = num_analisys(draws = theta, 
                                  names = c('mu_h', 'phi_h', 's_h', 
                                            'mu_a', 'phi_a', 'ls_a', 
                                            'v'),
                                  digits = 4,
                                  hdp = TRUE
  )
  # error
  err3[ ,it ] = matrix(summary3[[ it ]][, 1], ncol = 1) - theta_vdd
  # prob
  x1 = as.numeric( summary1[[ it ]][, 3] < theta_vdd )
  x2 = as.numeric( summary1[[ it ]][, 4] > theta_vdd )
  prob.piv = matrix( round( 0.5 * (x1 + x2), 0 ), ncol = 1)
  prob3 = prob3 + prob.piv
  
  if( it == M ) time = Sys.time() - time
}

time

vies1 = matrix( apply( err1, MARGIN = 1, mean ), ncol = 1 )
reqm1 = matrix( apply( err1^2, MARGIN = 1, mean ), ncol = 1 )
prob1 = prob1 / M
data1 = cbind( vies1, reqm1, prob1 )
data1 = data.frame( data1 )
colnames( data1 ) = c('vies1', 'reqm1', 'prob.cob1')

vies2 = matrix( apply( err2, MARGIN = 1, mean ), ncol = 1 )
reqm2 = matrix( apply( err2^2, MARGIN = 1, mean ), ncol = 1 )
prob2 = prob2 / M
data2 = cbind( vies2, reqm2, prob2 )
data2 = data.frame( data2 )
colnames( data2 ) = c('vies2', 'reqm2', 'prob.cob2')

vies3 = matrix( apply( err3, MARGIN = 1, mean ), ncol = 1 )
reqm3 = matrix( apply( err3^2, MARGIN = 1, mean ), ncol = 1 )
prob3 = prob3 / M
data3 = cbind( vies3, reqm3, prob3 )
colnames( data3 ) = c('vies3', 'reqm3', 'prob.cob3')

data = cbind( data1, data2, data3 )
row.names( data ) = c('mu_h', 'phi_h', 's_h', 
                       'mu_a', 'phi_a', 'ls_a', 
                       'v')
Summary = list( summary1 = summary1, 
                summary2 = summary2, 
                summary3 = summary3 )

save( data, Summary, file = paste0( dir, '/simulacao1.RData') )
