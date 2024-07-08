library(ggplot2)
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )

# Compilando o modelo Stan
model_stan = rstan::stan_model(file = 'pc_tvp_svm.stan')

# Reading data
ibovespa = read.csv('https://raw.githubusercontent.com/holtz27/svmsmn/main/aplica%C3%A7%C3%A3o/ibovespa/ibovespa.csv')
ibovespa = ibovespa[, c('Date', 'Adj.Close')]
ibovespa[, 2] = as.numeric( ibovespa[, 2] ) 
ibovespa = na.omit(ibovespa)
T = nrow(ibovespa)
log.ret = 100 * ( log( ibovespa[2:T, 2] ) - log( ibovespa[1:(T-1), 2] ) )
T = length( log.ret )
plot( log.ret, type = 'l' )

# RStan Sampling
p = 0.5
U = 0.5
draws = rstan::sampling(model_stan, 
                        data = list(T = length( log.ret ), 
                                    y = as.numeric( log.ret ),
                                    lambda = -log( p ) / sqrt( U ) ),
                        chains = 4,
                        warmup = 5e3,
                        iter = 5e3 + 500,
                        cores = 4
)

# Análise numérica
round(
  rstan::summary( draws )$summary[c('tau', 'mu', 'phi', 's2_h', 's2_a'), 
                                  c('mean', '2.5%', '97.5%', 'n_eff', 'Rhat')],
  3
)

# Análise gráfica
x = rstan::extract( draws )
draws_tau = x$tau
draws_mu = x$mu
draws_phi = x$phi
draws_s2_h = x$s2_h
draws_s2_a = log( x$s2_a )
Draws = matrix( draws_tau, nrow = 1 )
Draws = rbind( Draws, draws_mu )
Draws = rbind( Draws, draws_phi )
Draws = rbind( Draws, draws_s2_h )
Draws = rbind( Draws, draws_s2_a )
trace_plots(Draws, names = c(expression(tau), 
                             expression(mu), 
                             expression(phi), 
                             expression(sigma[h]^{2}), 
                             expression(sigma[a]^{2}))
)

# a
a_hmc = x$a
draws_a = apply( a_hmc, MARGIN = 2, mean )
a_min = apply( a_hmc, MARGIN = 2, quantile, probs = c(0.025, 0.975) )[1, ] 
a_max = apply( a_hmc, MARGIN = 2, quantile, probs = c(0.025, 0.975) )[2, ]
data = matrix(c(draws_a, a_min, a_max), ncol = 3)
data = data.frame( data )
data = cbind( c(1:T), data )
names(data) = c('time', 'rw.hat', 'a.min','a.max')
g1 = ggplot(data)
g1 = g1 + geom_line(aes(time, rw.hat), color = 'red', alpha = 1.0, linewidth = 0.5)
g1 = g1 + geom_ribbon(aes(x = time, ymax = a.max, ymin = a.min), 
                      fill = 'grey70', alpha = 0.2)
g1 = g1 + geom_line(aes(time, rw.hat), linewidth = 0.75)
g1 = g1 + geom_hline(yintercept = 0, linetype = "dashed")
g1 = g1 + theme_test() + xlab('t') + ylab('Impactos')
g1 = g1 + theme(axis.title.x = element_text(size = 26),
                axis.title.y = element_text(size = 26),
                axis.text.x = element_text(size = 24),
                axis.text.y = element_text(size = 26))
g1
plot(draws_a, type = 'l')

### VaR
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/VaR.R' )

horizon = 30
VaR_hat = rep(0, horizon)

for(it in 1:horizon){
  
  if(it == 1) time.init = Sys.time()
  
  y.now = log.ret[ 1:(T - horizon + it) ]
  draws = rstan::sampling(model_stan, 
                          data = list(T = length( y.now ), 
                                      y = as.numeric( y.now ),
                                      lambda = -log( p ) / sqrt( U ) ),
                          chains = 4,
                          warmup = 5e3,
                          iter = 5e3 + 500,
                          cores = 4
  )
  x = rstan::extract( draws )
  draws_tau = x$tau
  draws_mu = x$mu
  draws_phi = x$phi
  draws_s2_h = x$s2_h
  draws_s2_a = x$s2_a
  theta = matrix( draws_tau, nrow = 1 )
  theta = rbind( theta, draws_mu )
  theta = rbind( theta, draws_phi )
  theta = rbind( theta, draws_s2_h )
  theta = rbind( theta, draws_s2_a )
  h_hmc = x$h
  a_hmc = x$a
  # Evaluation VaR
  VaR_hat[ it ] = VaR(h_T = h_hmc[, dim( h_hmc )[2]],
                      a_T = a_hmc[, dim( h_hmc )[2]],
                      theta_hmc = theta )
  
  if(it == horizon) time.end = Sys.time() - time.init
  
}
time.end
plot( tail(log.ret, horizon),
      main = 'Pc prior',
      ylim = c(-.5, .5)
)
lines( VaR_hat )
