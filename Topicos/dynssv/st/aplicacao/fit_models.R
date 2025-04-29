source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R')
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/Topicos/waic.R')

# compiling stan models
model.dir = '~/topicos/st/models'
out.dir = '~/topicos/st/aplication'

summary = list('static','ig', 'pcp', 'exp')
waic.criterion = matrix(0, nrow = 1, ncol = 4)
colnames(waic.criterion) = c('static', 'ig', 'pcp', 'exp')

warmup = 3e2
iter = 2e2

for(i in 1:1){
  
  ##################
  # fitting model0 #
  ##################
  path = paste0(model.dir, '/static.stan')
  model_stan = rstan::stan_model(file=path)
  draws = rstan::sampling(model_stan, 
                          data = list(T=length(log.ret), 
                                      y=as.numeric(log.ret)),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iter,
                          cores = 1)
  x = rstan::extract(draws, 
                     pars=c('mu','phi_h','s_h','h','a','v','mu_t','sigma_t'))
  theta = rbind(x$mu, x$phi_h, x$s_h, x$a, x$v)
  summary$static = num_analisys(draws = theta, 
                                names = c('mu','phi_h','s_h','a','v'),
                                digits = 4,
                                hdp = TRUE)
  # Plots
  pdf(paste0(out.dir,'/static.pdf'), width=20, height=10)
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c(expression(mu), 
                        expression(phi[h]),
                        expression(sigma[h]),
                        expression(alpha),
                        expression(nu)))
  dev.off()
  h_hat1=apply(x$h, MARGIN=2, mean)
  # waic
  WAIC=waic(data=log.ret, mut = x$mu_t, st = x$sigma_t)$estimates['waic', 1]
  waic.criterion[1]=WAIC
  
  ##################
  # fitting model1 #
  ##################
  path = paste0( model.dir, '/ig.stan')
  model_stan = rstan::stan_model(file = path)
  draws = rstan::sampling(model_stan, 
                          data = list(T = length(log.ret), 
                                      y = as.numeric(log.ret)),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iter,
                          cores = 1
  )
  x = rstan::extract(draws, 
                     pars = c('mu','phi_h','s_h','h','a','s_a','v','mu_t','sigma_t'))
  theta = rbind(x$mu, x$phi_h, x$s_h, x$s_a, x$v)
  summary$ig = num_analisys(draws = theta, 
                            names = c('mu', 'phi_h', 's_h', 's_a','v'),
                            digits=4,
                            hdp=TRUE)
  # Plots
  pdf(paste0(out.dir,'/ig.pdf'), width=20, height=10)
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c(expression(mu), 
                        expression(phi[h]),
                        expression(sigma[h]),
                        expression(sigma[a]),
                        expression(nu)) )
  dev.off()
  
  pdf(paste0(out.dir, '/ig_a.pdf'), 
      width = 20, height = 10 )
  a_hat = apply( x$a, MARGIN = 2, mean )
  a_min = apply( x$a, MARGIN = 2, quantile, probs = 0.025 )
  a_max = apply( x$a, MARGIN = 2, quantile, probs = 0.975 )
  plot( a_hat, type = 'l' , ylim = c(min(a_min), max(a_max)) )
  lines( a_min, type = 'l', lty = 2 )
  lines( a_max, type = 'l', lty = 2 )
  abline( h = 0, col = 'grey', lty = 2)
  dev.off()
  
  h_hat2 = apply(x$h, MARGIN = 2, mean)
  # waic
  WAIC = waic(data = log.ret, mut = x$mu_t, st = x$sigma_t)$estimates['waic', 1]
  waic.criterion[2] = WAIC
  
  ##################
  # fitting model2 #
  ##################
  path = paste0(model.dir, '/pcp.stan')
  model_stan = rstan::stan_model(file = path)
  draws = rstan::sampling(model_stan, 
                          data = list(T = length( log.ret ), 
                                      y = as.numeric( log.ret ),
                                      lambda = -log(0.5)/0.5),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iter,
                          cores = 1)
  x = rstan::extract(draws, 
                     pars = c('mu','phi_h','s_h','h','s_a','a','v','mu_t','sigma_t'))
  theta = rbind(x$mu, x$phi_h, x$s_h, x$s_a, x$v)
  # Numeric Analysis
  summary$pcp = num_analisys(draws = theta, 
                             names = c('mu','phi_h','s_h','s_a','v'),
                             digits = 4,
                             hdp = TRUE)
  # Plots
  pdf(paste0(out.dir,'/pcp.pdf'), width=20, height=10)
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c(expression(mu), 
                        expression(phi[h]),
                        expression(sigma[h]),
                        expression(sigma[a]),
                        expression(nu)))
  dev.off()
  
  pdf(paste0(out.dir, '/pcp_a.pdf'), 
      width = 20, height = 10 )
  a_hat = apply( x$a, MARGIN = 2, mean )
  a_min = apply( x$a, MARGIN = 2, quantile, probs = 0.025 )
  a_max = apply( x$a, MARGIN = 2, quantile, probs = 0.975 )
  plot( a_hat, type = 'l' , ylim = c(min(a_min), max(a_max)) )
  lines( a_min, type = 'l', lty = 2 )
  lines( a_max, type = 'l', lty = 2 )
  abline( h = 0, col = 'grey', lty = 2)
  dev.off()
  
  h_hat3 = apply(x$h, MARGIN=2, mean)
  # waic
  WAIC = waic(data=log.ret, mut=x$mu_t, st=x$sigma_t)$estimates['waic', 1]
  waic.criterion[3] = WAIC
  
  ##################
  # fitting model3 #
  ##################
  path = paste0(model.dir, '/exp.stan')
  model_stan = rstan::stan_model(file = path)
  draws = rstan::sampling(model_stan, 
                         data = list(T=length(log.ret), 
                                     y=as.numeric(log.ret)),
                         chains = 1,
                         warmup = warmup,
                         iter = warmup + iter,
                         cores = 1)
  x = rstan::extract(draws, 
                     pars = c('mu','phi_h','s_h','h','s_a','a','v','mu_t','sigma_t'))
  theta = rbind(x$mu, x$phi_h, x$s_h, x$s_a, x$v)
  # Numeric Analysis
  summary$exp =  num_analisys(draws = theta, 
                              names = c('mu','phi_h','s_h','s_a','v'),
                              digits = 4,
                              hdp = TRUE)
  # Plots
  pdf(paste0(out.dir,'/exp.pdf'), width=20, height=10)
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c(expression(mu), 
                        expression(phi[h]),
                        expression(sigma[h]),
                        expression(sigma[a]),
                        expression(nu)
              ) 
  )
  dev.off()
  
  pdf(paste0(out.dir, '/exp_a.pdf'), 
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
  WAIC = waic(data=log.ret, mut=x$mu_t, st=x$sigma_t)$estimates['waic', 1]
  waic.criterion[ 4 ] = WAIC
  
  # Save
  save(
    summary, 
    waic.criterion,
    h_hat1, h_hat2, h_hat3, h_hat4,
    file = paste0(dir_out, 'nasdaq.RData')
  )  
}

summary
#waic.criterion = rbind(waic.criterion,
#                       as.numeric(waic.criterion[1, ] == min( waic.criterion)))
waic.criterion

# Figure Volatilities
plot(dates[-1], abs(log.ret), 
     col = 'gray', 
     xlab = '', ylab = '|Return|',
     type = 'l', cex.axis = 1.5, cex.lab = 1.5)
lines(dates[-1], exp(0.5 * h_hat1), lwd = 2, lty = 1 )
lines(dates[-1], exp(0.5 * h_hat2), lwd = 2, lty = 2, col = 'green' )
lines(dates[-1], exp(0.5 * h_hat3), lwd = 2, lty = 1, col = 'red' )
lines(dates[-1], exp(0.5 * h_hat4), lwd = 2, lty = 2, col = 'orange' )
legend('topright', 
       legend = c('Static', 'IG', 'PCP', 'Exp'), 
       col = c('black', 'green', 'red', 'orange'), 
       lty = c(1, 2, 1, 2), 
       lwd = 2, 
       bty = 'n')
