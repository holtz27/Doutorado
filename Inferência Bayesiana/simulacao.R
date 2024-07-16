path = '~/Doutorado/Disciplinas/Inf.Bayesiana/projeto0/simulacao/'
# Compilando o modelo Stan
model_stan = rstan::stan_model(file = 'pc_tvp_svm.stan')

# Simulando dados
tau = 0.1
mu = -9
phi = 0.985
s2_h = 0.025
T = 1e3
seeds = NULL
Summary = list()

# Cenários
s2_a = c( 0.0, 0.05, 0.5, 1.0 )
P = c( 0.01, 0.25, 0.5 )
U = c( 0.1, 0.25, 0.5 )

# Replicas por cenário
N = 2

# Running...
for( u in U ){
  for( s in s2_a ){
    for( p in P ){
      for( n in 1:N ){
        seed = sample( 1:1e6, 1 )
        set.seed( seed )
        seeds = cbind( seeds, seed )
        # Simulando dados
        y = h = a = matrix(0, nrow = T, ncol = 1)
        a[1] = -0.05 #rnorm(1, mean = -2, sd = 1)
        h[1] = mu + sqrt( s2_h / (1 - phi * phi) ) * rnorm( 1 )
        y[1] = tau + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )
        for( t in 2:T ){
          a[t] = a[t-1] + sqrt( s ) * rnorm( 1 )
          h[t] = mu + phi * (h[t-1] - mu) + sqrt( s2_h ) * rnorm( 1 )
          y[t] = tau + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
        }
        # RStan Sampling
        draws = rstan::sampling(model_stan, 
                                data = list(T = length( y ), 
                                            y = as.numeric( y ),
                                            lambda = -log( p ) / sqrt( u ) ),
                                chains = 2,
                                warmup = 1e3,
                                iter = 1e3 + 500,
                                cores = 2
        )
        Summary[[ n ]] = rstan::summary( draws )$summary[c('tau', 'mu', 'phi', 's2_h', 's2_a'), 
                                                  c('mean', '2.5%', '97.5%', 'n_eff', 'Rhat')]
      }
      save( Summary, 
            file = paste0( path, 'U_', u, '/' , 'cenario', s , '_' , p, '.RData' ) 
      )
    }
  }
}

#### Análise
library( ggplot2 )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )
# LOAD
load('~/Doutorado/Disciplinas/Inf.Bayesiana/projeto0/simulacao/U_0.1/cenario0_0.01.RData')
# Análise numérica
summary = rstan::summary( draws )$summary[c('tau', 'mu', 'phi', 's2_h', 's2_a'), 
                                          c('mean', '2.5%', '97.5%', 'n_eff', 'Rhat')]
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
theta = rbind( theta, log(draws_s2_a) )
# Análise gráfica
trace_plots(theta, names = c(expression(tau), 
                             expression(mu), 
                             expression(phi), 
                             expression(sigma[h]^{2}), 
                             expression( log(xi) )
)
)
# h
h_hmc = x$h
draws_h = apply( h_hmc, MARGIN = 2, mean )
h_min = apply( h_hmc, MARGIN = 2, quantile, probs = c(0.025, 0.975) )[1, ] 
h_max = apply( h_hmc, MARGIN = 2, quantile, probs = c(0.025, 0.975) )[2, ]
data = matrix(c(h, draws_h, h_min, h_max), ncol = 4)
data = data.frame( data )
data = cbind( c(1:length(h)), data )
names(data) = c('time', 'log.vol', 'log.vol.hat', 'h.min','h.max')
g0 = ggplot(data)
g0 = g0 + geom_line(aes(time, log.vol), color = 'red', 
                    alpha = 1.0, linewidth = 0.5)
g0 = g0 + geom_ribbon(aes(x = time, ymax = h.max, ymin = h.min), 
                      fill = 'grey70', alpha = 0.2)
g0 = g0 + geom_line(aes(time, log.vol.hat), linewidth = 0.75)
g0 = g0 + theme_test() + xlab('t') + ylab('Log-volatilidade')
g0 = g0 + theme(axis.title.x = element_text(size = 26),
                axis.title.y = element_text(size = 26),
                axis.text.x = element_text(size = 24),
                axis.text.y = element_text(size = 26))
# a
a_hmc = x$a
draws_a = apply( a_hmc, MARGIN = 2, mean )
a_min = apply( a_hmc, MARGIN = 2, quantile, probs = c(0.025, 0.975) )[1, ] 
a_max = apply( a_hmc, MARGIN = 2, quantile, probs = c(0.025, 0.975) )[2, ]
data = matrix(c(a, draws_a, a_min, a_max), ncol = 4)
data = data.frame( data )
data = cbind( c(1:length(a)), data )
names(data) = c('time', 'rw', 'rw.hat', 'a.min','a.max')
g1 = ggplot(data)
g1 = g1 + geom_line(aes(time, rw), color = 'red', 
                    alpha = 1.0, linewidth = 0.5)
g1 = g1 + geom_ribbon(aes(x = time, ymax = a.max, ymin = a.min), 
                      fill = 'grey70', alpha = 0.2)
g1 = g1 + geom_line(aes(time, rw.hat), linewidth = 0.75)
g1 = g1 + geom_hline(yintercept = 0, linetype = "dashed")
g1 = g1 + theme_test() + xlab('t') + ylab('Impactos')
g1 = g1 + theme(axis.title.x = element_text(size = 26),
                axis.title.y = element_text(size = 26),
                axis.text.x = element_text(size = 24),
                axis.text.y = element_text(size = 26))
gridExtra::grid.arrange(g0, g1, nrow = 2)
