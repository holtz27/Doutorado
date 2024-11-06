# function
rtnorm = function(n){
  u = runif(n)
  return( qnorm( 0.5 * (u + 1) ) )
}
source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R')
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
# Obtém o caminho completo do arquivo
localenvir = rstudioapi::getActiveDocumentContext()
path = localenvir$path
# Para obter apenas o diretório do arquivo
dir = dirname(path)
path = paste0( dir, '/dynssv_st.stan')
model_stan = rstan::stan_model(file = path)

set.seed(1648723)
T = 3e3
# log-volatility
mu_h = 0
phi_h = 0.99
s_h = 0.1
# dynskew
phi_a = 0.99
s_a = 0.0 # { 0, 0.05, 0.1 }

a = delta = omega = rep(0 , T)
a[1] = 0
for(t in 2:T) a[t] = phi_a * a[t-1] + s_a * rnorm(1)
delta = a / sqrt(1 + a*a)

v = 10

k1 = sqrt(0.5 * v) * gamma(0.5*(v-1)) / gamma(0.5 * v)
k2 = v / (v-2);
omega = 1 / sqrt( k2 - 2 * ( delta * k1 )^2 / pi )
mean_st = - sqrt(2 / pi) * delta * omega * k1
W = rtnorm(T)
U = rgamma(T, shape = 0.5 * v, rate = 0.5* v)
y = h = rep(0, T)
h[ 1 ] = mu_h + s_h / sqrt( 1 - phi_h * phi_h ) * rnorm(1)
for(t in 2:T) h[ t ] = mu_h + phi_h * (h[ t-1 ] - mu_h) + s_h * rnorm(1)
mu_t = mean_st + omega * delta * W * exp(0.5 * h) / sqrt( U )
sigma_t = omega * sqrt(1 - ( delta )^2) * exp(0.5 * h) / sqrt( U )
y = mu_t + sigma_t * rnorm(T)
par( mfrow = c(2, 2))
plot( y , type = 'l' )
hist( y, breaks = 40, probability = TRUE )
abline(v = mean( y ), col = 'red', lty = 2)
plot( a, type = 'l', col = 'red' )
par( mfrow = c(1, 1) )
data_summary = matrix(c( mean( y ),
                         sd( y ),
                         min( y ),
                         max( y ),
                         moments::skewness( y ),
                         moments::kurtosis( y ) ), nrow = 1)
colnames( data_summary ) = c('mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round( data_summary, digits = 4 )

################
### Sampling ###
################
draws = rstan::sampling(model_stan, 
                        data = list(T = length( y ), 
                                    y = as.numeric( y ),
                                    lambda1 = 9.473,
                                    lambda2 = -log( 0.5 ) / sqrt( 0.5 )
                        ),
                        chains = 1,
                        warmup = 3e3,
                        iter = 3e3 + 2e3,
                        cores = 1
)
x = rstan::extract(draws, 
                   pars = c('mu', 'phi_h', 's_h', 'phi_a', 'ls_a', 'a1', 'h', 'a', 'v')
                   )
theta = matrix(x$mu, nrow = 1)
theta = rbind(theta, x$phi_h, x$s_h, x$phi_a, x$ls_a, x$a1, x$v)
# Plots
#pdf( 'pcp_traces.pdf', width = 20, height = 10 )
trace_plots(theta,
            burn = 0, lags = 1,
            names = c(expression(mu), 
                      expression(phi[h]),
                      expression(sigma[h]),
                      expression(phi[a]),
                      expression(log(sigma[a])),
                      expression(a[1]),
                      expression(nu)
                      ) 
            )
#dev.off()

# Numeric Analysis
num_analisys(draws = theta, 
             names = c('mu_h', 'phi_h', 's_h', 
                       'phi_a', 'ls_a', 'a1',
                       'v'),
             digits = 4,
             hdp = TRUE
)

# h
#pdf( 'pcp_h.pdf', width = 20, height = 10 )
draws_h = x$h
h_hat = apply( draws_h, MARGIN = 2, mean )
h_min = apply( draws_h, MARGIN = 2, quantile, probs = 0.025 )
h_max = apply( draws_h, MARGIN = 2, quantile, probs = 0.975 )
plot( h, type = 'l', col = 'red', ylim = c(min(h_min), max(h_max)) )
lines( h_hat, type = 'l')
lines( h_min, type = 'l', lty = 2 )
lines( h_max, type = 'l', lty = 2 )
#dev.off()

# a
#pdf( '2ig_a.pdf', width = 20, height = 10 )
draws_a = x$a
a_hat = apply( draws_a, MARGIN = 2, mean )
a_min = apply( draws_a, MARGIN = 2, quantile, probs = 0.025 )
a_max = apply( draws_a, MARGIN = 2, quantile, probs = 0.975 )
plot( a_hat, type = 'l', ylim = c(min(a_min), max(a_max)) )
lines( a_min, type = 'l', lty = 2 )
lines( a_max, type = 'l', lty = 2 )
lines( a, type = 'l', col = 'red')
#dev.off()





