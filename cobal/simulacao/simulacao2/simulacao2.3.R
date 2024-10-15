# Install packages if necessary
if( !require(invgamma) ) install.packages('invgamma', repos = 'http://cran.us.r-project.org')
if( !require(rstan) ) install.packages('rstan', repos = 'http://cran.us.r-project.org')
if( !require(coda) ) install.packages('coda', repos = 'http://cran.us.r-project.org')
# function
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/forecast_vol.R')
dir = paste0( getwd(), '/Doutorado/Eventos/cobal/' )
path = paste0( dir, 'modelos/static.stan')
model_stan1 = rstan::stan_model(file = path)

# Data Settings
T = 1e3
b = 0.1
mu = 0
phi = 0.95 #0.95, 0.99
sigma = 0.15 #0.1, 0.15
xi = 0.1 # 0.01, 0.1, 1.0

summary = list()

u = 0.95
q = invgamma::qinvgamma(u, shape = 4.5, rate = 0.065)
lambda = -log(1 - u) / sqrt( q )

M = 1
m = 1
saver = M / m

k = 5
err = matrix(nrow = k, ncol = M)

warmup = 1e2
iters = 1e2

for(it in 1:M){
  if( (M %% m) != 0 ) stop('M precisar ser divisivel por m!')
  if( it == 1 ) time = Sys.time()
  
  #Data
  y = h = a = matrix(0, nrow = T, ncol = 1)
  a[1] = -0.05 
  h[1] = mu + (sigma / sqrt( (1 - phi * phi) )) * rnorm( 1 )
  y[1] = b + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )
  for( t in 2:T ){
    a[t] = a[t-1] + sqrt( xi ) * rnorm( 1 )
    h[t] = mu + phi * (h[t-1] - mu) + sigma * rnorm( 1 )
    y[t] = b + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
  }
  
  ### Sampling1
  y.train = y[1:(T-k)]
  h.test = matrix( h[(T-k+1):T], ncol = 1 )
  
  draws = rstan::sampling(model_stan1, 
                          data = list(T = length( y.train ), 
                                      y = as.numeric( y.train ),
                                      lambda = lambda
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s_h, x$a)
  summary[[ it ]] = num_analisys(draws = theta, 
                                 names = c('b', 
                                           'mu', 
                                           'phi_h', 
                                           's_h', 
                                           'a'),
                                 digits = 4,
                                 hdp = TRUE
  )
  
  h_new = forecast.vol(h_T = apply(x$h, MARGIN = 1, tail, 1),
                       mu = x$mu,
                       phi = x$phi,
                       sigma = x$s_h,
                       k = k
  )$h_new
  err[, it] = exp( h_new ) - exp( h.test )
  
  # saver
  if( it == saver ){
    
    save( summary, err,
          file = paste0( dir, 'simulacao/simulacao2/simulacao2.3.RData') )
    saver = saver + M / m
    
  }
  if( it == M ) time = Sys.time() - time
}

time




