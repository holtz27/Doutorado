# function
rtnorm = function(n){
  u = runif(n)
  return( qnorm( 0.5 * (u+1) ) )
}

source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R')
# Obtém o caminho completo do arquivo
localenvir = rstudioapi::getActiveDocumentContext()
path = localenvir$path
# Para obter apenas o diretório do arquivo
dir = dirname(path)
path = paste0( dir, '/sv_sn.stan')
model_stan = rstan::stan_model(file = path)

#set.seed(164872)
T = 1e3
# log-volatility
mu_h = 0
phi_h = 0.99
sigma_h = 0.1

# dynskew
a1 = -0.5
mu_a = 0
phi_a = 0.995
sigma_a = 0.01
a = delta = omega = rep(0 , T)
a[1] = a1
for(t in 2:T) a[t] = mu_a + phi_a * (a[t-1] - mu_a) + sigma_a * rnorm(1)
delta = a / sqrt(1 + a*a)

# drawing data
omega = 1 / sqrt( 1 - 2 * delta * delta / pi )
mean_sn = - sqrt(2 / pi) * delta * omega 
W = rtnorm(T)
y = h = rep(0, T)
h[ 1 ] = mu_h + sigma_h / sqrt( 1 - phi_h * phi_h )
for(t in 2:T) h[ t ] = mu_h + phi_h * (h[ t-1 ] - mu_h) + sigma_h * rnorm(1)
mu_t = mean_sn + omega * delta * W * exp(0.5 * h) 
sigma_t = omega * sqrt(1 - ( delta )^2) * exp(0.5 * h) 
y = mu_t + sigma_t * rnorm(T)
par( mfrow = c(1, 2))
plot( y , type = 'l' )
hist( y, breaks = 40, probability = TRUE )
abline(v = mean( y ), col = 'red', lty = 2)
par( mfrow = c(1, 1) )
data_summary = matrix(c( mean( y ),
                         sd( y ),
                         min( y ),
                         max( y ),
                         moments::skewness( y ),
                         moments::kurtosis( y ) ), nrow = 1)
colnames( data_summary ) = c('mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round( data_summary, digits = 4 )

### Sampling
draws = rstan::sampling(model_stan, 
                        data = list(T = length( y ), 
                                    y = as.numeric( y ),
                                    lambda1 = 2,
                                    lambda2 = 5
                                    ),
                        chains = 4,
                        warmup = 1e3,
                        iter = 1e3 + 0.5e3,
                        cores = 4
)
x = rstan::extract( draws )
draws_mu = x$mu
draws_phi_h = x$phi_h
draws_s_h = x$s_h
draws_a1 = x$a1
draws_phi_a = x$phi_a
draws_s_a = x$s_a
theta = matrix( draws_mu, nrow = 1 )
theta = rbind( theta, draws_phi_h )
theta = rbind( theta, draws_s_h )
theta = rbind( theta, draws_a1 )
theta = rbind( theta, draws_phi_a )
theta = rbind( theta, draws_s_a )
# Plots
#pdf( 'pcp_traces.pdf', width = 20, height = 10 )
trace_plots(theta,
            burn = 0, lags = 1,
            names = c('mu', 'phi_h', 's_h', 'a1', 'phi_a', 's_a') )
#dev.off()

# Numeric Analysis
summary = rstan::summary( draws )$summary[c('mu', 'phi_h', 's_h', 'a1', 'phi_a', 's_a'), 
                                          c('mean', '2.5%', '97.5%', 'n_eff', 'Rhat')]
round( summary, 4 )

# a
draws_a = x$a
#draws_a = draws_a[, 900:T]
a_hat = apply( draws_a, MARGIN = 2, mean )
a_min = apply( draws_a, MARGIN = 2, quantile, probs = 0.025 )
a_max = apply( draws_a, MARGIN = 2, quantile, probs = 0.975 )
plot( a_hat, type = 'l' , ylim = c(-2, 2) )
lines( a_min, type = 'l', lty = 2 )
lines( a_max, type = 'l', lty = 2 )
lines( a, type = 'l', col = 'red')

# h
draws_h = x$h
h_hat = apply( draws_h, MARGIN = 2, mean )
h_min = apply( draws_h, MARGIN = 2, quantile, probs = 0.025 )
h_max = apply( draws_h, MARGIN = 2, quantile, probs = 0.975 )
plot( h_hat, type = 'l', ylim = c(-10, 10) )
lines( h_min, type = 'l', lty = 2 )
lines( h_max, type = 'l', lty = 2 )
lines( h, type = 'l', col = 'red')



