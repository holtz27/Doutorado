# Install packages if necessary
if( !require(invgamma) ) install.packages('invgamma', repos = 'http://cran.us.r-project.org')
if( !require(rstan) ) install.packages('rstan', repos = 'http://cran.us.r-project.org')
if( !require(coda) ) install.packages('coda', repos = 'http://cran.us.r-project.org')

model_stan1 = rstan::stan_model(file = 'pcp.stan')
model_stan2 = rstan::stan_model(file = 'ig.stan')
model_stan3 = rstan::stan_model(file = 'exp.stan')

# Data Settings
T = 1e3
b = 0.1
mu = 0
phi = 0.99 #0.95, 0.99
sigma = 0.1
xi = 0.01 # 0.01, 0.25, 1.0

theta_vdd = matrix(c(b, mu, phi, sigma, log(xi) ), ncol = 1)
summary1 = summary2 = summary3 = list()

u = 0.95
q = invgamma::qinvgamma(u, shape = 4.5, rate = 0.065)
lambda = -log(1 - u) / sqrt( q )

M = 1
m = 1
saver = M / m

warmup = 1e3
iters = 1e3

for(it in 1:M){
  if( (M %% m) != 0 ) stop('M precisar ser divisivel por m!')
  if( it == 1 ) time = Sys.time()
  #Data
  y = h = a = matrix(0, nrow = T, ncol = 1)
  a[1] = -0.05 
  h[1] = mu + (sigma / sqrt( (1 - phi * phi) )) * rnorm( 1 )
  y[1] = b + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )
  for( t in 2:T ){
    a[t] = a[t-1] + sqrt( xi ) * rnorm( 1 )
    h[t] = mu + phi * (h[t-1] - mu) + sigma * rnorm( 1 )
    y[t] = b + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
  }
  
  ### Sampling1
  draws = rstan::sampling(model_stan1, 
                          data = list(T = length( y ), 
                                      y = as.numeric( y ),
                                      lambda = lambda
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$ls_a)
  summary1[[ it ]] = num_analisys(draws = theta, 
                                  names = c('b', 
                                            'mu', 
                                            'phi_h', 
                                            's_h', 
                                            'ls_a'),
                                  digits = 4,
                                  hdp = TRUE
  )
  
  ### Sampling2
  draws = rstan::sampling(model_stan2, 
                          data = list(T = length( y ), 
                                      y = as.numeric( y )
                                      ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$ls_a)
  summary2[[ it ]] = num_analisys(draws = theta, 
                                  names = c('b', 
                                            'mu', 
                                            'phi_h', 
                                            's_h', 
                                            'ls_a'),
                                  digits = 4,
                                  hdp = TRUE
  )
  
  ### Sampling3
  draws = rstan::sampling(model_stan3, 
                          data = list(T = length( y ), 
                                      y = as.numeric( y )
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$ls_a)
  summary3[[ it ]] = num_analisys(draws = theta, 
                                  names = c('b', 
                                            'mu', 
                                            'phi_h', 
                                            's_h', 
                                            'ls_a'),
                                  digits = 4,
                                  hdp = TRUE
  )
  
  # saver
  if( it == saver ){
    
    Summary = list( summary1 = summary1, 
                    summary2 = summary2, 
                    summary3 = summary3 )
    
    save( Summary, theta_vdd, file = paste0( dir, '/simulacao1.RData') )
    saver = saver + M / m
    
  }
  if( it == M ) time = Sys.time() - time
}

time
