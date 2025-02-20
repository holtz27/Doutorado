library( doParallel )

source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/forecast_vol.R')
model.dir = paste0('C:/Users/8936381/Documents/Doutorado/Eventos/cobal/models')
path = paste0(model.dir, '/ig.stan')
model_stan1 = rstan::stan_model(file = path)
path = paste0(model.dir, '/pcp.stan')
model_stan2 = rstan::stan_model(file = path)

out.dir = 'C:/Users/8936381/Documents/Doutorado/Eventos/cobal/simulation'

seed = sample(1:1e6, 1)
set.seed( seed ) 
T = 2e2
b = 0.1
mu = -1
# (phi, sigma) \in { (0.95, 0.225), (0.99, 0.1) }
phi = 0.95 
sigma = 0.225
xi = 0.00 # 0, 0.005, 0.05
a1 = -0.1
theta_vdd = matrix(c(mu, phi, sigma, xi), ncol = 1)

M = 3
k = 5
err1 = err2 = matrix(nrow = k, ncol = M)
lambda = -log(0.5)/0.5

warmup = 2e2
iters = 1e2

num_cores = detectCores() - 1
cl = makeCluster( num_cores )
registerDoParallel( cl )

result = foreach(it = 1:M, .packages = c('rstan')) %dopar% {
  
  time = Sys.time()
  
  #Data
  y = h = a = numeric(T + k)
  a[1] = a1
  h[1] = mu + (sigma / sqrt((1 - phi * phi))) * rnorm(1)
  y[1] = b + a[1] * exp(h[1]) + exp(0.5 * h[1]) * rnorm(1)
  for( t in 2:(T+k) ){
    N = 50
    A = a[t-1] + sqrt(xi) * rnorm(N)
    a[t] = mean(A)
    #a[t] = a[t-1] + sqrt( xi ) * rnorm( 1 )
    h[t] = mu + phi * (h[t-1] - mu) + sigma * rnorm( 1 )
    y[t] = b + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
  }
  
  y.train = y[1:T]
  h.test = matrix( h[(T + 1):(T + k)], ncol = 1 )
  #####
  
  ### Sampling1
  test = 1
  while( any(test > 0, is.nan(test), is.na(test)) ){
    draws = rstan::sampling(model_stan1, 
                            data = list(T = length( y.train ), 
                                        y = as.numeric( y.train )),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iters,
                            cores = 1)
    x = rstan::extract( draws, pars = c('b','mu','phi','s_h','s_a','h') )
    theta = rbind(x$b, x$mu, x$phi, x$s_h, x$s_a)
    s = num_analisys(draws = theta, 
                     names = c('b','mu','phi','s_h','s_a'),
                     digits = 4,
                     hdp = TRUE)
    test = sum( abs( s[ , 'CD'] ) > 1.96 )
  }
  summary1 = s
  h_new = forecast.vol(h_T = x$h[, T], #apply(x$h, MARGIN = 1, tail, 1),
                       mu = x$mu,
                       phi = x$phi,
                       sigma = x$s_h,
                       k = k
  )$h_new
  err1 = exp( h_new ) - exp( h.test ) 
  
  ### Sampling2
  test = 1
  while( any(test > 0, is.nan(test), is.na(test)) ){
    draws = rstan::sampling(model_stan2, 
                            data = list(T = length( y.train ), 
                                        y = as.numeric( y.train ),
                                        lambda = lambda),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iters,
                            cores = 1)
    x = rstan::extract( draws, pars = c('b','mu','phi','s_h','s_a','h') )
    theta = rbind(x$b, x$mu, x$phi, x$s_h, x$s_a)
    s = num_analisys(draws = theta, 
                     names = c('b','mu','phi','s_h','s_a'),
                     digits = 4,
                     hdp = TRUE)
    test = sum( abs( s[ , 'CD'] ) > 1.96 )
  }
  summary2 = s
  h_new = forecast.vol(h_T = x$h[, T], #apply(x$h, MARGIN = 1, tail, 1),
                       mu = x$mu,
                       phi = x$phi,
                       sigma = x$s_h,
                       k = k
  )$h_new
  err2 = exp( h_new ) - exp( h.test )
  
  time = Sys.time() - time
  
  list(summary1 = summary1,
       summary2 = summary2,
       err1 = err1,
       err2 = err2,
       time = time)

}

stopCluster( cl )
result

#save( result, theta_vdd,file = paste0( out.dir, '/simulacao_xi_', xi, '.RData') )

