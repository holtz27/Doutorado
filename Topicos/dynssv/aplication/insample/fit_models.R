source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R')
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/Topicos/waic.R')

# compiling stan models
model.dir = c('~/topicos/sn/models', 
              '~/topicos/st/models',
              '~/topicos/ss/models')

out.dir = c('~/topicos/sn/aplication/sp500',
            '~/topicos/st/aplication/sp500',
            '~/topicos/ss/aplication/sp500')

summary = list('SN'=list('ig', 'pcp', 'exp'), 
               'St'=list('ig', 'pcp', 'exp'), 
               'SS'=list('ig', 'pcp', 'exp'))


waic.criterion = matrix(0, nrow = 3, ncol = 3)
row.names(waic.criterion) = c('SN', 'St', 'SS')
colnames(waic.criterion) = c('ig', 'pcp', 'exp')
waic.criterion 


warmup=5e1
iter=2e1
log.ret = log.ret-mean(log.ret)

for(i in 1:3){
  
  if(i==1){
    pars=c('mu','phi_h','s_h','h','a','s_a','mu_t','sigma_t')
    names = c('mu','phi_h','s_h','s_a')
  }else{
    pars=c('mu','phi_h','s_h','h','a','s_a','mu_t','sigma_t','v')
    names = c('mu','phi_h','s_h','s_a','v')
  }
  
  
  ##################
  # fitting model1 #
  ##################
  path = paste0(model.dir[i], '/ig.stan')
  model_stan = rstan::stan_model(file = path)
  draws = rstan::sampling(model_stan, 
                          data = list(T = length(log.ret), 
                                      y = as.numeric(log.ret)),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iter,
                          cores = 1)
  x = rstan::extract(draws, pars=pars)
  if(i==1){
    theta = rbind(x$mu, x$phi_h, x$s_h, x$s_a)  
  }else{
    theta = rbind(x$mu, x$phi_h, x$s_h, x$s_a, x$v)
  }
  summary[i][[1]][[1]] = num_analisys(draws=theta, names=names, digits=4, hdp=TRUE)
  # Plots
  pdf(paste0(out.dir[i],'/ig.pdf'), width=20, height=10)
  trace_plots(theta, dens=FALSE, names=names)
  dev.off()
  
  pdf(paste0(out.dir[i], '/ig_a.pdf'), 
      width = 20, height = 10 )
  a_hat = apply( x$a, MARGIN = 2, mean )
  a_min = apply( x$a, MARGIN = 2, quantile, probs = 0.025 )
  a_max = apply( x$a, MARGIN = 2, quantile, probs = 0.975 )
  plot( a_hat, type = 'l' , ylim = c(min(a_min), max(a_max)) )
  lines( a_min, type = 'l', lty = 2 )
  lines( a_max, type = 'l', lty = 2 )
  abline( h = 0, col = 'grey', lty = 2)
  dev.off()
  
  h_hat1 = apply(x$h, MARGIN=2, mean)
  # waic
  WAIC = waic(data = log.ret, mut = x$mu_t, st = x$sigma_t)$estimates['waic', 1]
  waic.criterion[i,1] = WAIC
  ##################
  # fitting model2 #
  ##################
  path = paste0(model.dir[i], '/pcp.stan')
  model_stan = rstan::stan_model(file = path)
  draws = rstan::sampling(model_stan, 
                          data = list(T = length(log.ret), 
                                      y = as.numeric(log.ret),
                                      lambda = -log(0.5)/0.5),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iter,
                          cores = 1)
  x = rstan::extract(draws, pars=pars)
  if(i==1){
    theta = rbind(x$mu, x$phi_h, x$s_h, x$s_a)  
  }else{
    theta = rbind(x$mu, x$phi_h, x$s_h, x$s_a, x$v)
  }
  # Numeric Analysis
  summary[i][[1]][[2]] = num_analisys(draws=theta, names=names, digits=4, hdp=TRUE)
  # Plots
  pdf(paste0(out.dir[i],'/pcp.pdf'), width=20, height=10)
  trace_plots(theta, dens=FALSE, names=names)
  dev.off()
  
  pdf(paste0(out.dir[i], '/pcp_a.pdf'), 
      width = 20, height = 10 )
  a_hat = apply( x$a, MARGIN = 2, mean )
  a_min = apply( x$a, MARGIN = 2, quantile, probs = 0.025 )
  a_max = apply( x$a, MARGIN = 2, quantile, probs = 0.975 )
  plot( a_hat, type = 'l' , ylim = c(min(a_min), max(a_max)) )
  lines( a_min, type = 'l', lty = 2 )
  lines( a_max, type = 'l', lty = 2 )
  abline( h = 0, col = 'grey', lty = 2)
  dev.off()
  
  h_hat2 = apply(x$h, MARGIN=2, mean)
  # waic
  WAIC = waic(data=log.ret, mut=x$mu_t, st=x$sigma_t)$estimates['waic', 1]
  waic.criterion[i,2]  = WAIC
  ##################
  # fitting model3 #
  ##################
  path = paste0(model.dir[i], '/exp.stan')
  model_stan = rstan::stan_model(file = path)
  draws = rstan::sampling(model_stan, 
                         data = list(T=length(log.ret), 
                                     y=as.numeric(log.ret)),
                         chains = 1,
                         warmup = warmup,
                         iter = warmup + iter,
                         cores = 1)
  x = rstan::extract(draws, pars=pars)
  if(i==1){
    theta = rbind(x$mu, x$phi_h, x$s_h, x$s_a)  
  }else{
    theta = rbind(x$mu, x$phi_h, x$s_h, x$s_a, x$v)
  }
  # Numeric Analysis
  summary[i][[1]][[3]] = num_analisys(draws=theta, names=names, digits=4, hdp=TRUE)
  # Plots
  pdf(paste0(out.dir[i],'/exp.pdf'), width=20, height=10)
  trace_plots(theta, dens=FALSE, names=names)
  dev.off()
  
  pdf(paste0(out.dir[i], '/exp_a.pdf'), 
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
  WAIC = waic(data=log.ret, mut=x$mu_t, st=x$sigma_t)$estimates['waic', 1]
  waic.criterion[i,3] = WAIC
  # Save
  save(summary, waic.criterion, h_hat1, h_hat2, h_hat3, 
       file = paste0(out.dir[i], '/sp500.RData'))  
}

summary
waic.criterion









# Figure Volatilities
plot(dates[-1], abs(log.ret), 
     col = 'gray', 
     xlab = '', ylab = '|Return|',
     type = 'l', cex.axis = 1.5, cex.lab = 1.5)
lines(dates[-1], exp(0.5 * h_hat1), lwd = 2, lty = 1 )
lines(dates[-1], exp(0.5 * h_hat2), lwd = 2, lty = 2, col = 'green' )
lines(dates[-1], exp(0.5 * h_hat3), lwd = 2, lty = 1, col = 'red' )
#lines(dates[-1], exp(0.5 * h_hat4), lwd = 2, lty = 2, col = 'orange' )
legend('topright', 
       legend = c('Static', 'IG', 'PCP', 'Exp'), 
       col = c('black', 'green', 'red', 'orange'), 
       lty = c(1, 2, 1, 2), 
       lwd = 2, 
       bty = 'n')
