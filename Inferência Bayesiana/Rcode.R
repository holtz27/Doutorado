# Simulando dados
tau = 0.1
mu = -10.0
phi = 0.98
s2_h = 0.1
s2_a = 0.1

a1 = -1.0
T = 1e3

set.seed( 8936381 )
y = h = a = matrix(0, nrow = T, ncol = 1)
a[1] = a1
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

# Carregar o modelo Stan
model_stan = rstan::stan_model(file = 'tvp_svm.stan')
# Dados para o modelo
draws = rstan::sampling(model_stan, 
                        data = list(T = length( y ), 
                                    y = as.numeric( y ),
                                    a1 = 0,
                                    lambda = 1.0),
                        chains = 4,
                        iter = 5e3
                        )
# Visualizar os resultados
rstan::summary( draws )$summary[c('tau', 
                                  'mu',
                                  'phi',
                                  's2_h',
                                  's2_a'), 
                                c('mean', 
                                  '2.5%',
                                  '97.5%',
                                  'n_eff',
                                  'Rhat')] 
# Plots
x = rstan::extract( draws )
plot( x$tau, type = 'l' )
plot( x$mu, type = 'l' )
plot( x$phi, type = 'l' )
plot( x$s2_h, type = 'l' )
plot( x$s2_a, type = 'l' )
