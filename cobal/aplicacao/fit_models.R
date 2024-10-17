source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/waic.R')
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')

# Compilando o modelo Stan
dir = paste0(getwd(), '/Doutorado/Eventos/cobal/modelos/')
path = paste0(dir, 'static.stan')
model_stan1 = rstan::stan_model(file = path)
path = paste0(dir, 'ig.stan')
model_stan2 = rstan::stan_model(file = path)
path = paste0(dir, 'pcp.stan')
model_stan3 = rstan::stan_model(file = path)
path = paste0(dir, 'exp.stan')
model_stan4 = rstan::stan_model(file = path)

U = invgamma::qinvgamma(0.975, shape = 4.5, rate = 0.065)
lambda = -log(0.025) / sqrt( U )

summary = list('static', 'ig', 'pcp', 'exp')
waic.criterion = matrix(0, nrow = 1, ncol = 4)
colnames( waic.criterion ) = c('static', 'ig', 'pcp', 'exp')

warmup = 1e3
iter = 2e3

for( i in 1:1 ){
  # fitting model1
  draws = rstan::sampling(model_stan1, 
                          data = list(T = length( log.ret ), 
                                      y = as.numeric( log.ret )
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iter,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$a)
  summary$static = num_analisys(draws = theta, 
                                  names = c('b', 
                                            'mu', 
                                            'phi_h', 
                                            's_h', 
                                            'a'),
                                  digits = 4,
                                  hdp = FALSE
  )
  # Plots
  pdf( 'static.pdf', width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'a') )
  dev.off()
  WAIC = waic(data = log.ret, 
              draws = rbind(x$b, t(x$h), x$a ),
              dyn = FALSE)
  waic.criterion[1] = WAIC$estimates['waic', 1]
  
  # fitting model2
  draws = rstan::sampling(model_stan2, 
                          data = list(T = length( log.ret ), 
                                      y = as.numeric( log.ret )
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iter,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$ls_a)
  # Numeric Analysis
  summary$ig = num_analisys(draws = theta, 
                            names = c('b', 
                                      'mu', 
                                      'phi_h', 
                                      's_h', 
                                      'ls_a'),
                            digits = 4,
                            hdp = FALSE
  )
  # Plots
  pdf( 'ig.pdf', width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'ls_a')
              )
  dev.off()
  WAIC = waic(data = log.ret, 
              draws = rbind(x$b, t(x$h), t(x$a)),
              dyn = FALSE)
  waic.criterion[ 2 ] = WAIC$estimates['waic', 1]
  
  # fitting model3
  draws = rstan::sampling(model_stan3, 
                          data = list(T = length( log.ret ), 
                                      y = as.numeric( log.ret ),
                                      lambda = lambda
                                      ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iter,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$ls_a)
  # Numeric Analysis
  summary$pcp = num_analisys(draws = theta, 
                             names = c('b', 
                                       'mu', 
                                       'phi_h', 
                                       's_h', 
                                       'ls_a'),
                             digits = 4,
                             hdp = FALSE
  )
  # Plots
  pdf( 'pcp.pdf', width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'ls_a') )
  dev.off()
  WAIC = waic(data = log.ret, 
              draws = rbind(x$b, t(x$h), t(x$a)),
              dyn = FALSE)
  waic.criterion[ 3 ] = WAIC$estimates['waic', 1]
  
  # fitting model3
  draws = rstan::sampling(model_stan4, 
                          data = list(T = length( log.ret ), 
                                      y = as.numeric( log.ret )
                                      ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iter,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$ls_a)
  # Numeric Analysis
  summary$exp = num_analisys(draws = theta, 
                             names = c('b', 
                                       'mu', 
                                       'phi_h', 
                                       's_h', 
                                       'ls_a'),
                             digits = 4,
                             hdp = FALSE
  )
  # Plots
  pdf( 'exp.pdf', width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'ls_a') )
  dev.off()
  WAIC = waic(data = log.ret, 
              draws = rbind(x$b, t(x$h), t(x$a)),
              dyn = FALSE)
  waic.criterion[ 4 ] = WAIC$estimates['waic', 1]
  
}

summary
waic.criterion
