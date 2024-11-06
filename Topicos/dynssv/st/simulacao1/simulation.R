library( doParallel )

source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')

rtnorm = function(n){
  u = runif(n)
  return( qnorm( 0.5 * (u + 1) ) )
}

model.dir = '~/Doutorado/Pesquisa/Papers/Topicos/DS-SV/Code/st/models'
path = paste0( model.dir, '/ig.stan')
model_stan1 = rstan::stan_model(file = path)
path = paste0(model.dir, '/pcp.stan')
model_stan2 = rstan::stan_model(file = path)
path = paste0(model.dir, '/exp.stan')
model_stan3 = rstan::stan_model(file = path)

out.dir = '~/Doutorado/Pesquisa/Papers/Topicos/DS-SV/Code/st/simulation'

seed = sample(1:1e6, 1)
set.seed(seed)
T = 1.5e3

# log-volatility
mu = 0
phi_h = 0.99
s_h = 0.1

# dynskew
phi_a = 0.99
s_a = 0.0 # { 0, 0.05, 0.1 }

v = 10

theta_vdd = matrix(c(mu, phi_h, s_h, phi_a, s_a, v), ncol = 1)
lambda = -log( 0.5 ) / sqrt( 0.5 )

M = 3
warmup = 2e2
iters = 1e2

num_cores = detectCores() - 1
cl = makeCluster( num_cores )
registerDoParallel( cl )

result = foreach(it = 1:M, .packages = c('rstan')) %dopar% {
  
  time = Sys.time()
  
  #Data
  a = delta = omega = rep(0 , T)
  a[1] = 0
  for(t in 2:T) a[t] = phi_a * a[t-1] + s_a * rnorm(1)
  delta = a / sqrt(1 + a*a)
  
  k1 = sqrt(0.5 * v) * gamma(0.5*(v-1)) / gamma(0.5 * v)
  k2 = v / (v-2);
  omega = 1 / sqrt( k2 - 2 * ( delta * k1 )^2 / pi )
  mean_st = - sqrt(2 / pi) * delta * omega * k1
  W = rtnorm(T)
  U = rgamma(T, shape = 0.5 * v, rate = 0.5* v)
  y = h = rep(0, T)
  h[ 1 ] = mu + s_h / sqrt( 1 - phi_h * phi_h ) * rnorm(1)
  for(t in 2:T) h[ t ] = mu + phi_h * (h[ t-1 ] - mu) + s_h * rnorm(1)
  mu_t = mean_st + omega * delta * W * exp(0.5 * h) / sqrt( U )
  sigma_t = omega * sqrt(1 - ( delta )^2) * exp(0.5 * h) / sqrt( U )
  y = mu_t + sigma_t * rnorm(T)
  #####
  
  ### Sampling1
  test = 1
  while( any(test > 0, is.nan(test), is.na(test)) ){
    draws = rstan::sampling(model_stan1, 
                            data = list(T = length( y ), 
                                        y = as.numeric( y ),
                                        lambda1 = 9.473
                            ),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iters,
                            cores = 1
    )
    x = rstan::extract(draws, 
                       pars = c('mu', 'phi_h', 's_h', 
                                'phi_a', 's_a', 'ls_a', 'a1', 'h', 'a', 'v')
    )
    theta = matrix(x$mu, nrow = 1)
    theta = rbind(theta, x$phi_h, x$s_h, x$phi_a, x$s_a, x$ls_a, x$a1, x$v)
    s = num_analisys(draws = theta, 
                     names = c('mu', 'phi_h', 's_h', 
                               'phi_a', 's_a', 'ls_a', 'a1',
                               'v'),
                     digits = 4,
                     hdp = TRUE
    )
    test = sum( abs( s[ , 'CD'] ) > 1.96 )
  }
  summary1 = s
  
  ### Sampling2
  test = 1
  while( any(test > 0, is.nan(test), is.na(test)) ){
    draws = rstan::sampling(model_stan2, 
                            data = list(T = length( y ), 
                                        y = as.numeric( y ),
                                        lambda1 = 9.473,
                                        lambda2 = lambda
                            ),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iters,
                            cores = 1
    )
    x = rstan::extract(draws, 
                       pars = c('mu', 'phi_h', 's_h', 
                                'phi_a', 's_a', 'ls_a', 'a1', 'h', 'a', 'v')
    )
    theta = matrix(x$mu, nrow = 1)
    theta = rbind(theta, x$phi_h, x$s_h, x$phi_a, x$s_a, x$ls_a, x$a1, x$v)
    s = num_analisys(draws = theta, 
                     names = c('mu', 'phi_h', 's_h', 
                               'phi_a', 's_a', 'ls_a', 'a1',
                               'v'),
                     digits = 4,
                     hdp = TRUE
    )
    test = sum( abs( s[ , 'CD'] ) > 1.96 )
  }
  summary2 = s
  
  ### Sampling3
  test = 1
  while( any(test > 0, is.nan(test), is.na(test)) ){
    draws = rstan::sampling(model_stan3, 
                            data = list(T = length( y ), 
                                        y = as.numeric( y ),
                                        lambda1 = 9.473
                            ),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iters,
                            cores = 1
    )
    x = rstan::extract(draws, 
                       pars = c('mu', 'phi_h', 's_h', 
                                'phi_a', 's_a', 'ls_a', 'a1', 'h', 'a', 'v')
    )
    theta = matrix(x$mu, nrow = 1)
    theta = rbind(theta, x$phi_h, x$s_h, x$phi_a, x$s_a, x$ls_a, x$a1, x$v)
    s = num_analisys(draws = theta, 
                     names = c('mu', 'phi_h', 's_h', 
                               'phi_a', 's_a', 'ls_a', 'a1',
                               'v'),
                     digits = 4,
                     hdp = TRUE
    )
    test = sum( abs( s[ , 'CD'] ) > 1.96 )
  }
  summary3 = s
  
  time = Sys.time() - time
  
  list(summary1 = summary1,
       summary2 = summary2,
       summary3 = summary3,
       time = time)
  
}

stopCluster( cl )
result

save( result, theta_vdd,
      file = paste0( out.dir, '/simulacao_sa_', s_a, '.RData')
      )

