library(ggplot2)
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )

# Compilando o modelo Stan
model_stan = rstan::stan_model(file = 'tvp_svm.stan')

# Simulando dados
tau = 0.1
mu = -5
phi = 0.98
s2_h = 0.25
s2_a = 0.25

T = 5e2

#set.seed( 8936381 )
y = h = a = matrix(0, nrow = T, ncol = 1)
a[1] = -1.0
h[1] = mu + sqrt( s2_h / (1 - phi * phi) ) * rnorm( 1 )
y[1] = tau + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )

for( t in 2:T ){
  a[t] = a[t-1] + sqrt( s2_a ) * rnorm( 1 )
  h[t] = mu + phi * (h[t-1] - mu) + sqrt( s2_h ) * rnorm( 1 )
  y[t] = tau + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
}

plot(y, type = 'l')
plot(a, type = 'l')
plot(exp(h), type = 'l')


# Smapling
p = 0.01
U = 0.2
draws = rstan::sampling(model_stan, 
                        data = list(T = length( y ), 
                                    y = as.numeric( y ),
                                    lambda = -log( p ) / sqrt( U ) ),
                        chains = 4,
                        iter = 5e2
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
draws_s2_a = x$s2_a
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
# h
draws_h = apply( x$h, MARGIN = 2, mean )
h_min = apply( x$h, MARGIN = 2, quantile, probs = c(0.025, 0.975) )[1, ] 
h_max = apply( x$h, MARGIN = 2, quantile, probs = c(0.025, 0.975) )[2, ]
data = matrix(c(h, draws_h, h_min, h_max), ncol = 4)
data = data.frame( data )
data = cbind( c(1:T), data )
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
draws_a = apply( x$a, MARGIN = 2, mean )
a_min = apply( x$a, MARGIN = 2, quantile, probs = c(0.025, 0.975) )[1, ] 
a_max = apply( x$a, MARGIN = 2, quantile, probs = c(0.025, 0.975) )[2, ]
data = matrix(c(a, draws_a, a_min, a_max), ncol = 4)
data = data.frame( data )
data = cbind( c(1:T), data )
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

gridExtra::grid.arrange(g0, g1, nrow=1, widths=c(1.5, 1))
