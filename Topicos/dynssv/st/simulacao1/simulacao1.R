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

set.seed( 8936381 ) # 164872
T = 1e3
# log-volatility
mu_h = 0
phi_h = 0.99
sigma_h = 0.15
# dynskew
mu_a = 0
phi_a = 0.99
sigma_a = 0.01
v = 10

theta_vdd = matrix(c(mu_h, phi_h, sigma_h, 
                     mu_a, phi_a, log(sigma_a), 
                     v), ncol = 1)
summary1 = summary2 = summary3 = list()

M = 100
m = 5
saver = M / m

warmup = 2e3
iters = 1e3

for(it in 1:M){
  
  if( (M %% m) != 0 ) stop('M precisar ser divisivel por m!')
  
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
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$mu_h, nrow = 1)
  theta = rbind(theta, x$phi_h, x$s_h, x$mu_a, x$phi_a, x$ls_a, x$v)
  summary1[[ it ]] = num_analisys(draws = theta, 
                                  names = c('mu_h', 'phi_h', 's_h', 
                                            'mu_a', 'phi_a', 'ls_a', 
                                            'v'),
                                  digits = 4,
                                  hdp = TRUE
                                  )

  ### Sampling2
  draws = rstan::sampling(model_stan2, 
                          data = list(T = length( y ), 
                                      y = as.numeric( y ),
                                      lambda1 = 1.55,
                                      lambda2 = 14.34
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$mu_h, nrow = 1)
  theta = rbind(theta, x$phi_h, x$s_h, x$mu_a, x$phi_a, x$ls_a, x$v)
  summary2[[ it ]] = num_analisys(draws = theta, 
                                  names = c('mu_h', 'phi_h', 's_h', 
                                            'mu_a', 'phi_a', 'ls_a', 
                                            'v'),
                                  digits = 4,
                                  hdp = TRUE
  )
  
  ### Sampling3
  draws = rstan::sampling(model_stan3, 
                          data = list(T = length( y ), 
                                      y = as.numeric( y ),
                                      lambda1 = 1.55,
                                      lambda2 = 14.34
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract( draws )
  theta = matrix(x$mu_h, nrow = 1)
  theta = rbind(theta, x$phi_h, x$s_h, x$mu_a, x$phi_a, x$ls_a, x$v)
  summary3[[ it ]] = num_analisys(draws = theta, 
                                  names = c('mu_h', 'phi_h', 's_h', 
                                            'mu_a', 'phi_a', 'ls_a', 
                                            'v'),
                                  digits = 4,
                                  hdp = TRUE
  )
  # saver
  if( it == saver ){
    
    Summary = list( summary1 = summary1, 
                    summary2 = summary2, 
                    summary3 = summary3 )
    
    save( Summary, theta_vdd, file = paste0( dir, '/simulacao1.1.RData') )
    saver = saver + M / m
    
  }
  if( it == M ) time = Sys.time() - time
  
}

time

# Analazing...
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/source/res.sim.R')
load('~/dynssv/st/simulacao1/simulacao1.1.RData')

s = Summary$summary3
y = res.sim( s, theta_vdd, med.abs = FALSE )
y$errors
y$metricas
