# function
rtnorm = function(n){
  u = runif(n)
  return( qnorm( 0.5 * (u + 1) ) )
}
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/forecast_vol.R')
localenvir = rstudioapi::getActiveDocumentContext()
path = localenvir$path
dir = dirname( path )
path = paste0( dir, '/static.stan')
model_stan1 = rstan::stan_model(file = path)
path = paste0( dir, '/ig.stan')
model_stan2 = rstan::stan_model(file = path)
path = paste0( dir, '/pcp.stan')
model_stan3 = rstan::stan_model(file = path)

set.seed(1648723)
T = 1e2
b = 0.1
mu = 0
phi = 0.95 #0.95, 0.99
sigma = 0.15 #0.1, 0.15
xi = 0 # 0, 0.1, 0.25

theta_vdd = matrix(c(mu, phi, sigma, xi), ncol = 1)
summary1 = summary2 = summary3 = list()
Time = 0

M = 2
m = 1
saver = M / m

k = 5
err1 = err2 = err3 = matrix(nrow = k, ncol = M)
lambda = -log( 0.5 ) / sqrt( 0.5 )

warmup = 3e2
iters = 1e2

for(it in 1:M){
  if( (M %% m) != 0 ) stop('M precisar ser divisivel por m!')
  time = Sys.time()
  
  #Data
  y = h = a = matrix(0, nrow = T, ncol = 1)
  a[1] = 0 
  h[1] = mu + (sigma / sqrt( (1 - phi * phi) )) * rnorm( 1 )
  y[1] = b + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )
  for( t in 2:T ){
    a[t] = a[t-1] + sqrt( xi ) * rnorm( 1 )
    h[t] = mu + phi * (h[t-1] - mu) + sigma * rnorm( 1 )
    y[t] = b + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
  }
  
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
                                        y = as.numeric( y.train )
                            ),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iters,
                            cores = 1
    )
    x = rstan::extract(draws)
    theta = matrix(x$b, nrow = 1)
    theta = rbind(theta, x$mu, x$phi, x$s_h, x$a)
    s = num_analisys(draws = theta, 
                     names = c('b', 'mu', 'phi', 's_h', 'a'),
                     digits = 4,
                     hdp = TRUE
    )
    test = sum( abs( s[ , 'CD'] ) > 1.96 )
    tent = tent + 1
  }
  summary1[[ it ]] = s
  h_new = forecast.vol(h_T = apply(x$h, MARGIN = 1, tail, 1),
                       mu = x$mu,
                       phi = x$phi,
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
                                        y = as.numeric( y.train )
                                        ),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iters,
                            cores = 1
    )
    x = rstan::extract(draws)
    theta = matrix(x$b, nrow = 1)
    theta = rbind(theta, x$mu, x$phi, x$s_h, x$s_a)
    s = num_analisys(draws = theta, 
                     names = c('b', 'mu', 'phi', 's_h', 's_a'),
                     digits = 4,
                     hdp = TRUE
    )
    test = sum( abs( s[ , 'CD'] ) > 1.96 )
    tent = tent + 1
  }
  summary2[[ it ]] = s
  h_new = forecast.vol(h_T = apply(x$h, MARGIN = 1, tail, 1),
                       mu = x$mu,
                       phi = x$phi,
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
                                        lambda = lambda
                            ),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iters,
                            cores = 1
    )
    x = rstan::extract(draws)
    theta = matrix(x$b, nrow = 1)
    theta = rbind(theta, x$mu, x$phi, x$s_h, x$s_a)
    s = num_analisys(draws = theta, 
                     names = c('b', 'mu', 'phi', 's_h', 's_a'),
                     digits = 4,
                     hdp = TRUE
    )
    test = sum( abs( s[ , 'CD'] ) > 1.96 )
    tent = tent + 1
  }
  summary3[[ it ]] = s
  h_new = forecast.vol(h_T = apply(x$h, MARGIN = 1, tail, 1),
                       mu = x$mu,
                       phi = x$phi,
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
          file = paste0( dir,
                         '/simulacao_xi_',
                         xi,
                         '.RData')
    )
    saver = saver + M / m
    
  }
  time = Sys.time() - time
  Time = Time + as.numeric(time, units = 'hours')
}

Time




