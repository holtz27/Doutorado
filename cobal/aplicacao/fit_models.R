source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R')
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/waic.R')
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')

# compiling stan models
dir = paste0(getwd(), '/Doutorado/Eventos/cobal/')
path = paste0(dir, 'models/static.stan')
model_stan1 = rstan::stan_model(file = path)
path = paste0(dir, 'models/ig.stan')
model_stan2 = rstan::stan_model(file = path)
path = paste0(dir, 'models/pcp.stan')
model_stan3 = rstan::stan_model(file = path)
path = paste0(dir, 'models/jeffrey.stan')
model_stan4 = rstan::stan_model(file = path)
path = paste0(dir, 'models/exp.stan')
model_stan5 = rstan::stan_model(file = path)

# outs directory
dir_out = paste0( dir, 'aplication/data/cobre/' )

summary = list('static', 'ig', 'pcp', 'jeffrey', 'exp')
waic.criterion = matrix(0, nrow = 1, ncol = 5)
colnames( waic.criterion ) = c('static', 'ig', 'pcp', 'jeffrey', 'exp')

warmup = 1e2
iter = 1e1

for( i in 1:1 ){
  ##################
  # fitting model1 #
  ##################
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
  pdf( paste0( dir_out, 'static.pdf'), 
       width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'a') )
  dev.off()
  h_hat1 = apply( x$h, MARGIN = 2, mean )
  # waic
  WAIC = waic(data = log.ret, 
              draws = rbind(x$b, t(x$h), x$a ),
              dyn = FALSE)
  waic.criterion[1] = WAIC$estimates['waic', 1]
  
  ##################
  # fitting model2 #
  ##################
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
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$ls_a, x$a1)
  # Numeric Analysis
  summary$ig = num_analisys(draws = theta, 
                            names = c('b', 
                                      'mu', 
                                      'phi_h', 
                                      's_h', 
                                      'ls_a',
                                      'a1'),
                            digits = 4,
                            hdp = FALSE
  )
  # Plots
  pdf( paste0( dir_out, 'ig.pdf'), 
       width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'ls_a', 'a1')
              )
  dev.off()
  
  pdf( paste0( dir_out, 'ig_a.pdf'), 
       width = 20, height = 10 )
  a_hat = apply( x$a, MARGIN = 2, mean )
  a_min = apply( x$a, MARGIN = 2, quantile, probs = 0.025 )
  a_max = apply( x$a, MARGIN = 2, quantile, probs = 0.975 )
  plot( a_hat, type = 'l' , ylim = c(min(a_min), max(a_max)) )
  lines( a_min, type = 'l', lty = 2 )
  lines( a_max, type = 'l', lty = 2 )
  abline( h = 0, col = 'grey', lty = 2)
  dev.off()
  
  h_hat2 = apply( x$h, MARGIN = 2, mean )
  # waic
  WAIC = waic(data = log.ret, 
              draws = rbind(x$b, t(x$h), t(x$a)),
              dyn = FALSE)
  waic.criterion[ 2 ] = WAIC$estimates['waic', 1]
  
  ##################
  # fitting model3 #
  ##################
  draws = rstan::sampling(model_stan3, 
                          data = list(T = length( log.ret ), 
                                      y = as.numeric( log.ret ),
                                      lambda = -log(0.5) / sqrt(0.5)
                                      ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iter,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$ls_a, x$a1)
  # Numeric Analysis
  summary$pcp = num_analisys(draws = theta, 
                             names = c('b', 
                                       'mu', 
                                       'phi_h', 
                                       's_h', 
                                       'ls_a',
                                       'a1'),
                             digits = 4,
                             hdp = FALSE
  )
  # Plots
  pdf( paste0( dir_out, 'pcp.pdf'), 
       width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'ls_a', 'a1') )
  dev.off()
  
  pdf( paste0( dir_out, 'pcp_a.pdf'), 
       width = 20, height = 10 )
  a_hat = apply( x$a, MARGIN = 2, mean )
  a_min = apply( x$a, MARGIN = 2, quantile, probs = 0.025 )
  a_max = apply( x$a, MARGIN = 2, quantile, probs = 0.975 )
  plot( a_hat, type = 'l' , ylim = c(min(a_min), max(a_max)) )
  lines( a_min, type = 'l', lty = 2 )
  lines( a_max, type = 'l', lty = 2 )
  abline( h = 0, col = 'grey', lty = 2)
  dev.off()
  h_hat3 = apply( x$h, MARGIN = 2, mean )
  # waic
  WAIC = waic(data = log.ret, 
              draws = rbind(x$b, t(x$h), t(x$a)),
              dyn = FALSE)
  waic.criterion[ 3 ] = WAIC$estimates['waic', 1]
  
  ##################
  # fitting model4 #
  ##################
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
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$ls_a, x$a1)
  # Numeric Analysis
  summary$jeffrey = num_analisys(draws = theta, 
                                 names = c('b', 
                                           'mu', 
                                           'phi_h', 
                                           's_h', 
                                           'ls_a',
                                           'a1'),
                                 digits = 4,
                                 hdp = FALSE
  )
  # Plots
  pdf( paste0( dir_out, 'jeffrey.pdf'), 
       width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'ls_a', 'a1') )
  dev.off()
  
  pdf( paste0( dir_out, 'jeffrey_a.pdf' ), 
       width = 20, height = 10 )
  a_hat = apply( x$a, MARGIN = 2, mean )
  a_min = apply( x$a, MARGIN = 2, quantile, probs = 0.025 )
  a_max = apply( x$a, MARGIN = 2, quantile, probs = 0.975 )
  plot( a_hat, type = 'l' , ylim = c(min(a_min), max(a_max)) )
  lines( a_min, type = 'l', lty = 2 )
  lines( a_max, type = 'l', lty = 2 )
  abline( h = 0, col = 'grey', lty = 2)
  dev.off()
  
  h_hat4 = apply( x$h, MARGIN = 2, mean )
  # waic
  WAIC = waic(data = log.ret, 
              draws = rbind(x$b, t(x$h), t(x$a)),
              dyn = FALSE)
  waic.criterion[ 4 ] = WAIC$estimates['waic', 1]
  
  ##################
  # fitting model5 #
  ##################
  draws = rstan::sampling(model_stan5, 
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
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$ls_a, x$a1)
  # Numeric Analysis
  summary$jeffrey = num_analisys(draws = theta, 
                                 names = c('b', 
                                           'mu', 
                                           'phi_h', 
                                           's_h', 
                                           'ls_a',
                                           'a1'),
                                 digits = 4,
                                 hdp = FALSE
  )
  # Plots
  pdf( paste0( dir_out, 'exp.pdf'), 
       width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'ls_a', 'a1') )
  dev.off()
  
  pdf( paste0( dir_out, 'exp_a.pdf' ), 
       width = 20, height = 10 )
  a_hat = apply( x$a, MARGIN = 2, mean )
  a_min = apply( x$a, MARGIN = 2, quantile, probs = 0.025 )
  a_max = apply( x$a, MARGIN = 2, quantile, probs = 0.975 )
  plot( a_hat, type = 'l' , ylim = c(min(a_min), max(a_max)) )
  lines( a_min, type = 'l', lty = 2 )
  lines( a_max, type = 'l', lty = 2 )
  abline( h = 0, col = 'grey', lty = 2)
  dev.off()
  
  h_hat5 = apply( x$h, MARGIN = 2, mean )
  # waic
  WAIC = waic(data = log.ret, 
              draws = rbind(x$b, t(x$h), t(x$a)),
              dyn = FALSE)
  waic.criterion[ 5 ] = WAIC$estimates['waic', 1]
  
  # Save
  save(
    summary, 
    waic.criterion,
    h_hat1, h_hat2, h_hat3, h_hat4, h_hat5,
    file = paste0(dir_out, 'out.RData')
       )  
}

summary
waic.criterion = rbind(waic.criterion,
                       as.numeric(waic.criterion[1, ] == min( waic.criterion)))
waic.criterion

# Figure Volatilities
plot(dates[-1], abs(log.ret), 
     col = 'gray', 
     xlab = '', ylab = '|Return|',
     type = 'l', cex.axis = 1.5, cex.lab = 1.5)
lines(dates[-1], exp(0.5 * h_hat1), lwd = 2, lty = 1 )
lines(dates[-1], exp(0.5 * h_hat2), lwd = 2, lty = 2, col = 'green' )
lines(dates[-1], exp(0.5 * h_hat3), lwd = 2, lty = 1, col = 'red' )
lines(dates[-1], exp(0.5 * h_hat4), lwd = 2, lty = 2, col = 'blue' )
lines(dates[-1], exp(0.5 * h_hat4), lwd = 2, lty = 2, col = 'orange' )
legend('topright', 
       legend = c('Static', 'IG', 'PCP', 'Jeffrey', 'Exp'), 
       col = c('black', 'green', 'red', 'blue', 'orange'), 
       lty = c(1, 2, 1, 2, 1), 
       lwd = 2, 
       bty = 'n')
