source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )
# Compilando o modelo Stan
model_stan1 = rstan::stan_model(file = 'tvp_svm.stan')
model_stan2 = rstan::stan_model(file = 'classic_tvp_svm.stan')
model_stan3 = rstan::stan_model(file = 'pcp_tvp_svm.stan')

U = invgamma::qinvgamma(0.975, shape = 4.5, rate = 0.065)
lambda = -log(0.025) / sqrt( U )

summary = list() 

for( i in 1:1 ){
  # fitting model1
  draws = rstan::sampling(model_stan1, 
                          data = list(T = length( log.ret ), 
                                      y = as.numeric( log.ret )
                          ),
                          chains = 4,
                          warmup = 2e3,
                          iter = 2e3 + 500,
                          cores = 4
  )
  x = rstan::extract( draws )
  draws_b = x$b
  draws_mu = x$mu
  draws_phi = x$phi
  draws_s2_h = x$s2_h
  draws_a1 = x$a1
  theta = matrix( draws_b, nrow = 1 )
  theta = rbind( theta, draws_mu )
  theta = rbind( theta, draws_phi )
  theta = rbind( theta, draws_s2_h )
  theta = rbind( theta, draws_a1 )
  # Plots
  pdf( 'tvp_traces.pdf', width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'a1') )
  dev.off()
  # Numeric Analysis
  summary[[1]] = rstan::summary( draws )$summary[c('b', 'mu', 'phi', 's2_h', 'a1'), 
                                                 c('mean', '2.5%', '97.5%', 'n_eff', 'Rhat')]
  
  # fitting model2
  draws = rstan::sampling(model_stan2, 
                          data = list(T = length( log.ret ), 
                                      y = as.numeric( log.ret )
                          ),
                          chains = 4,
                          warmup = 2e3,
                          iter = 2e3 + 500,
                          cores = 4
  )
  x = rstan::extract( draws )
  draws_b = x$b
  draws_mu = x$mu
  draws_phi = x$phi
  draws_s2_h = x$s2_h
  draws_a1 = x$a1
  theta = matrix( draws_b, nrow = 1 )
  theta = rbind( theta, draws_mu )
  theta = rbind( theta, draws_phi )
  theta = rbind( theta, draws_s2_h )
  theta = rbind( theta, draws_a1 )
  # Plots
  pdf( 'classic_traces.pdf', width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'ls_a', 'a1') )
  dev.off()
  # Numeric Analysis
  summary[[2]] = rstan::summary( draws )$summary[c('b', 'mu', 'phi', 's2_h', 'ls_a', 'a1'), 
                                                 c('mean', '2.5%', '97.5%', 'n_eff', 'Rhat')]
  
  # fitting model3
  draws = rstan::sampling(model_stan, 
                          data = list(T = length( log.ret ), 
                                      y = as.numeric( log.ret ),
                                      lambda = lambda
                                      ),
                          chains = 4,
                          warmup = 2e3,
                          iter = 2e3 + 500,
                          cores = 4
  )
  x = rstan::extract( draws )
  draws_b = x$b
  draws_mu = x$mu
  draws_phi = x$phi
  draws_s2_h = x$s2_h
  draws_a1 = x$a1
  theta = matrix( draws_b, nrow = 1 )
  theta = rbind( theta, draws_mu )
  theta = rbind( theta, draws_phi )
  theta = rbind( theta, draws_s2_h )
  theta = rbind( theta, draws_a1 )
  # Plots
  pdf( 'pcp_traces.pdf', width = 20, height = 10 )
  trace_plots(theta,
              burn = 0, lags = 1,
              names = c('b', 'mu', 'phi', 's2_h', 'ls_a', 'a1') )
  dev.off()
  # Numeric Analysis
  summary[[3]] = rstan::summary( draws )$summary[c('b', 'mu', 'phi', 's2_h', 'ls_a', 'a1'), 
                                                 c('mean', '2.5%', '97.5%', 'n_eff', 'Rhat')]
  
}





































