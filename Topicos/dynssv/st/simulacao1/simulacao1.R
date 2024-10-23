# function
rtnorm = function(n){
  u = runif(n)
  return( qnorm( 0.5 * (u + 1) ) )
}
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/forecast_vol.R')
localenvir = rstudioapi::getActiveDocumentContext()
path = localenvir$path
dir = dirname(path)
path = paste0( dir, '/dynssv_st1.stan')
model_stan1 = rstan::stan_model(file = path)
path = paste0( dir, '/dynssv_st2.stan')
model_stan2 = rstan::stan_model(file = path)
path = paste0( dir, '/dynssv_st3.stan')
model_stan3 = rstan::stan_model(file = path)

set.seed(1648723)
T = 1e3
# log-volatility
mu_h = 0
phi_h = 0.99
sigma_h = 0.1
# dynskew
mu_a = 0
phi_a = 0.99
sigma_a = 0.1
v = 10

theta_vdd = matrix(c(mu_h, phi_h, sigma_h,
                     mu_a, phi_a, sigma_a,
                     v), 
                   ncol = 1)
summary1 = summary2 = summary3 = list()
Time = 0

M = 100
m = 4
saver = M / m

k = 5
err1 = err2 = err3 = matrix(nrow = k, ncol = M)

warmup = 3e3
iters = 1e3

for(it in 1:M){
  if( (M %% m) != 0 ) stop('M precisar ser divisivel por m!')
  time = Sys.time()
  
  #cat( '\014' )
  #cat('Réplica: ', it, '/', M, '\n', '\n')
  #if( it == 1 ){
  #  cat('Tempo restante total estimado: Calculando...', '\n' )
  #}else{
  #  cat('Tempo restante total estimado: ', 
  #      round((M - it + 1) * (Time / (it-1)), 1), 'hs',
  #      '\n' )
  #}
  #####
  # Data
  a = delta = omega = rep(0 , T + k)
  a[1] = mu_a + sigma_a / sqrt( 1 - phi_a * phi_a ) * rnorm(1)
  for(t in 2:(T + k)) a[t] = mu_a + phi_a * (a[t-1] - mu_a) + sigma_a * rnorm(1)
  delta = a / sqrt(1 + a*a)
  k1 = sqrt(0.5 * v) * gamma(0.5*(v-1)) / gamma(0.5 * v)
  k2 = v / (v-2);
  omega = 1 / sqrt( k2 - 2 * ( delta * k1 )^2 / pi )
  mean_st = - sqrt(2 / pi) * delta * omega * k1
  W = rtnorm(T + k)
  U = rgamma(T + k, shape = 0.5 * v, rate = 0.5* v)
  y = h = rep(0, T + k)
  h[ 1 ] = mu_h + sigma_h / sqrt( 1 - phi_h * phi_h ) * rnorm(1)
  for(t in 2:(T + k)) h[ t ] = mu_h + phi_h * (h[ t-1 ] - mu_h) + sigma_h * rnorm(1)
  mu_t = mean_st + omega * delta * W * exp(0.5 * h) / sqrt( U )
  sigma_t = omega * sqrt(1 - ( delta )^2) * exp(0.5 * h) / sqrt( U )
  y = mu_t + sigma_t * rnorm(T + k)
  
  y.train = y[1:T]
  h.test = matrix( h[(T + 1):(T + k)], ncol = 1 )
  #####
  
  ### Sampling1
  test = tent = 1
  while( test > 0 || is.na(test) || is.nan(test) ){
    cat( '\014' )
    cat('Réplica: ', it, '/', M, '\n', '\n')
    if( it == 1 ){
      cat('Tempo restante total estimado: Calculando...', '\n', '\n' )
    }else{
      cat('Tempo restante total estimado: ', 
          round((M - it + 1) * (Time / (it-1)), 1), 'hs',
          '\n', '\n' )
    }
    cat( 'Modelo 1, tentetiva ', tent, '\n', '\n' )
    draws = rstan::sampling(model_stan1, 
                            data = list(T = length( y.train ), 
                                        y = as.numeric( y.train ),
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
    s = num_analisys(draws = theta, 
                     names = c('mu_h', 'phi_h', 's_h', 
                               'mu_a', 'phi_a', 'ls_a', 
                               'v'),
                     digits = 4,
                     hdp = TRUE
    )
    test = sum( abs( s[ , 'CD'] ) > 1.96 )
    tent = tent + 1
  }
  summary1[[ it ]] = s
  h_new = forecast.vol(h_T = apply(x$h, MARGIN = 1, tail, 1),
                       mu = x$mu_h,
                       phi = x$phi_h,
                       sigma = x$s_h,
                       k = k
                      )$h_new
  err1[, it] = exp( h_new ) - exp( h.test ) 
  cat( '\014' )
  
  ### Sampling2
  test = tent = 1
  while( test > 0 || is.na(test) || is.nan(test) ){
    cat( '\014' )
    cat('Réplica: ', it, '/', M, '\n', '\n')
    if( it == 1 ){
      cat('Tempo restante total estimado: Calculando...', '\n', '\n' )
    }else{
      cat('Tempo restante total estimado: ', 
          round((M - it + 1) * (Time / (it-1)), 1), 'hs',
          '\n', '\n' )
    }
    cat( 'Modelo 2, tentetiva ', tent, '\n', '\n' )
    draws = rstan::sampling(model_stan2, 
                            data = list(T = length( y.train ), 
                                        y = as.numeric( y.train ),
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
    s = num_analisys(draws = theta, 
                     names = c('mu_h', 'phi_h', 's_h', 
                               'mu_a', 'phi_a', 'ls_a', 
                               'v'),
                     digits = 4,
                     hdp = TRUE
    )
    test = sum( abs( s[ , 'CD'] ) > 1.96 )
    tent = tent + 1
  }
  summary2[[ it ]] = s
  h_new = forecast.vol(h_T = apply(x$h, MARGIN = 1, tail, 1),
                       mu = x$mu_h,
                       phi = x$phi_h,
                       sigma = x$s_h,
                       k = k
  )$h_new
  err2[, it] = exp( h_new ) - exp( h.test )
  cat( '\014' )
  
  ### Sampling3
  test = tent = 1
  while( test > 0 || is.na(test) || is.nan(test) ){
    cat( '\014' )
    cat('Réplica: ', it, '/', M, '\n', '\n')
    if( it == 1 ){
      cat('Tempo restante total estimado: Calculando...', '\n', '\n' )
    }else{
      cat('Tempo restante total estimado: ', 
          round((M - it + 1) * (Time / (it-1)), 1), 'hs',
          '\n', '\n' )
    }
    cat( 'Modelo 3, tentetiva ', tent, '\n', '\n' )
    draws = rstan::sampling(model_stan3, 
                            data = list(T = length( y.train ), 
                                        y = as.numeric( y.train ),
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
    s = num_analisys(draws = theta, 
                     names = c('mu_h', 'phi_h', 's_h', 
                               'mu_a', 'phi_a', 'ls_a', 
                               'v'),
                     digits = 4,
                     hdp = TRUE
    )
    test = sum( abs( s[ , 'CD'] ) > 1.96 )
    tent = tent + 1
  }
  summary3[[ it ]] = s
  h_new = forecast.vol(h_T = apply(x$h, MARGIN = 1, tail, 1),
                       mu = x$mu_h,
                       phi = x$phi_h,
                       sigma = x$s_h,
                       k = k
  )$h_new
  err3[, it] = exp( h_new ) - exp( h.test )
  
  # saver
  if( it == saver ){
    
    save( summary1,
          summary2,
          summary3,
          theta_vdd, 
          err1,
          err2,
          err3,
          file = 'simulacao_sa_0.1.RData'
          )
    saver = saver + M / m
    
  }
  time = Sys.time() - time
  Time = Time + as.numeric(time, units = 'hours')
}

Time




