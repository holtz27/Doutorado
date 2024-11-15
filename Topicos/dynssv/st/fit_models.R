source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R')
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/Topicos/waic.R')

# compiling stan models
model.dir = '~/Doutorado/Pesquisa/Papers/Topicos/DS-SV/Code/st/models'
path = paste0( model.dir, '/ig.stan')
model_stan1 = rstan::stan_model(file = path)
path = paste0(model.dir, '/pcp.stan')
model_stan2 = rstan::stan_model(file = path)
path = paste0(model.dir, '/exp.stan')
model_stan3 = rstan::stan_model(file = path)

out.dir = '~/Doutorado/Pesquisa/Papers/Topicos/DS-SV/Code/st/aplicacao'

summary = list('ig', 'pcp', 'exp')
waic.criterion = matrix(0, nrow = 1, ncol = 3)
colnames(waic.criterion) = c('ig', 'pcp', 'exp')

warmup = 3e3
iter = 2e3

for( i in 1:1 ){
  ##################
  # fitting model1 #
  ##################
  test = 1
  while( test > 0 || is.na(test) || is.nan(test) ){
    draws = rstan::sampling(model_stan1, 
                            data = list(T = length( log.ret ), 
                                        y = as.numeric( log.ret ),
                                        lambda1 = 9.473
                            ),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iter,
                            cores = 1
    )
    x = rstan::extract(draws, 
                       pars = c('mu', 'phi_h', 's_h', 'h',
                                'phi_a', 's_a', 'a', 'ls_a',
                                'v', 'mu_t', 'sigma_t')
    )
    theta = matrix(x$mu, nrow = 1)
    theta = rbind(theta, x$phi_h, x$s_h, x$phi_a, x$s_a, x$v)
    summary$ig = num_analisys(draws = theta, 
                                  names = c('mu', 'phi_h', 's_h', 
                                            'phi_a', 's_a', 
                                            'v'),
                                  digits = 4,
                                  hdp = TRUE
    )
    test = sum( abs( summary$ig[ , 'CD'] ) > 1.96 )
  }
  # Plots
  pdf( paste0( out.dir, 'ig.pdf'), width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c(expression(mu), 
                        expression(phi[h]),
                        expression(sigma[h]),
                        expression(phi[a]),
                        expression(sigma[a]),
                        expression(nu)
              ) 
  )
  dev.off()
  h_hat1 = apply( x$h, MARGIN = 2, mean )
  # waic
  WAIC = waic(data = log.ret, mut = x$mu_t, st = x$sigma_t)$estimates['waic', 1]
  waic.criterion[1] = WAIC
  
  ##################
  # fitting model2 #
  ##################
  test = 1
  while( test > 0 || is.na(test) || is.nan(test) ){
    draws = rstan::sampling(model_stan2, 
                            data = list(T = length( log.ret ), 
                                        y = as.numeric( log.ret ),
                                        lambda1 = 9.473,
                                        lambda2 = -log( 0.5 )/sqrt(0.5)
                            ),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iter,
                            cores = 1
    )
    x = rstan::extract(draws, 
                       pars = c('mu', 'phi_h', 's_h', 'h',
                                'phi_a', 's_a', 'a', 'ls_a',
                                'v', 'mu_t', 'sigma_t')
    )
    theta = matrix(x$mu, nrow = 1)
    theta = rbind(theta, x$phi_h, x$s_h, x$phi_a, x$s_a, x$v)
    # Numeric Analysis
    summary$pcp = num_analisys(draws = theta, 
                               names = c('mu', 'phi_h', 's_h', 
                                         'phi_a', 's_a', 
                                         'v'),
                               digits = 4,
                               hdp = TRUE
    )
    test = sum( abs( summary$ig[ , 'CD'] ) > 1.96 )
  }
  # Plots
  pdf( paste0( out.dir, 'pcp.pdf'), width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c(expression(mu), 
                        expression(phi[h]),
                        expression(sigma[h]),
                        expression(phi[a]),
                        expression(sigma[a]),
                        expression(nu)
              ) 
  )
  dev.off()
  
  #pdf( paste0( dir_out, 'ig_a.pdf'), 
  #     width = 20, height = 10 )
  #a_hat = apply( x$a, MARGIN = 2, mean )
  #a_min = apply( x$a, MARGIN = 2, quantile, probs = 0.025 )
  #a_max = apply( x$a, MARGIN = 2, quantile, probs = 0.975 )
  #plot( a_hat, type = 'l' , ylim = c(min(a_min), max(a_max)) )
  #lines( a_min, type = 'l', lty = 2 )
  #lines( a_max, type = 'l', lty = 2 )
  #abline( h = 0, col = 'grey', lty = 2)
  #dev.off()
  
  h_hat2 = apply( x$h, MARGIN = 2, mean )
  # waic
  WAIC = waic(data = y, mut = x$mu_t, st = x$sigma_t)$estimates['waic', 1]
  waic.criterion[ 2 ] = WAIC
  
  ##################
  # fitting model3 #
  ##################
  test = 1
  while( test > 0 || is.na(test) || is.nan(test) ){
    draws = rstan::sampling(model_stan3, 
                            data = list(T = length( log.ret ), 
                                        y = as.numeric( log.ret ),
                                        c = 0.1,
                                        d = 0.1
                            ),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iter,
                            cores = 1
    )
    x = rstan::extract(draws, 
                       pars = c('mu', 'phi_h', 's_h', 'h',
                                'phi_a', 's_a', 'a', 'ls_a',
                                'v', 'mu_t', 'sigma_t')
    )
    theta = matrix(x$mu, nrow = 1)
    theta = rbind(theta, x$phi_h, x$s_h, x$phi_a, x$s_a, x$v)
    # Numeric Analysis
    summary$exp =  num_analisys(draws = theta, 
                                names = c('mu', 'phi_h', 's_h', 
                                          'phi_a', 's_a', 
                                          'v'),
                                digits = 4,
                                hdp = TRUE
    )
    test = sum( abs( summary$pcp[ , 'CD'] ) > 1.96 )
  }
  
  # Plots
  pdf( paste0( out.dir, 'exp.pdf'), width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c(expression(mu), 
                        expression(phi[h]),
                        expression(sigma[h]),
                        expression(phi[a]),
                        expression(sigma[a]),
                        expression(nu)
              ) 
  )
  dev.off()
  
  #pdf( paste0( dir_out, 'pcp_a.pdf'), 
  #     width = 20, height = 10 )
  #a_hat = apply( x$a, MARGIN = 2, mean )
  #a_min = apply( x$a, MARGIN = 2, quantile, probs = 0.025 )
  #a_max = apply( x$a, MARGIN = 2, quantile, probs = 0.975 )
  #plot( a_hat, type = 'l' , ylim = c(min(a_min), max(a_max)) )
  #lines( a_min, type = 'l', lty = 2 )
  #lines( a_max, type = 'l', lty = 2 )
  #abline( h = 0, col = 'grey', lty = 2)
  #dev.off()
  h_hat3 = apply( x$h, MARGIN = 2, mean )
  # waic
  WAIC = waic(data = y, mut = x$mu_t, st = x$sigma_t)$estimates['waic', 1]
  waic.criterion[ 3 ] = WAIC
  
  # Save
  save(
    summary, 
    waic.criterion,
    h_hat1, h_hat2, h_hat3,
    file = paste0(dir_out, 'nasdaq.RData')
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
lines(dates[-1], exp(0.5 * h_hat5), lwd = 2, lty = 2, col = 'blue' )
lines(dates[-1], exp(0.5 * h_hat4), lwd = 2, lty = 2, col = 'orange' )
legend('topright', 
       legend = c('Static', 'IG', 'PCP', 'Jeffrey', 'Exp'), 
       col = c('black', 'green', 'red', 'blue', 'orange'), 
       lty = c(1, 2, 1, 2, 1), 
       lwd = 2, 
       bty = 'n')
