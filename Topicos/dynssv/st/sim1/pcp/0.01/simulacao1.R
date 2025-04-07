library(doParallel)

source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')

rtnorm = function(n){
  u = runif(n)
  return(qnorm(0.5*(u + 1)))
}

#path = 'ig.stan'
#model_stan1 = rstan::stan_model(file = path)
path = 'pcp.stan'
model_stan2 = rstan::stan_model(file = path)
#path = 'exp.stan'
#model_stan3 = rstan::stan_model(file = path)

#seed = sample(1:1e6, 1)
#set.seed(seed)
T = 1.5e3
# log-volatility
mu = 0
phi_h = 0.99
s_h = 0.1
# dynskew
s_a = 0.01 
v = 8

theta_vdd = matrix(c(mu, phi_h, s_h, s_a, v), ncol = 1)

M = 100
warmup = 5e3
iters = 2e3

num_cores = detectCores()-1
cl = makeCluster(num_cores)
registerDoParallel(cl)

result = foreach(it = 1:M, .packages = c('rstan')) %dopar% {
  
  time = Sys.time()
  
  #Data
  a = delta = omega = numeric(T)
  a[1] = 0
  for(t in 2:T) a[t] = a[t-1] + s_a*rnorm(1)
  delta = a/sqrt(1+a*a)
  k1 = sqrt(0.5*v)*gamma(0.5*(v-1))/gamma(0.5*v)
  k2 = v/(v-2);
  omega = 1/sqrt(k2-2*(delta*k1)^2/pi)
  mean_st = -sqrt(2/pi)*delta*omega*k1
  W = rtnorm(T)
  U = rgamma(T, shape = 0.5*v, rate = 0.5*v)
  y = h = numeric(T)
  h[ 1 ] = mu + s_h/sqrt(1-phi_h*phi_h)*rnorm(1)
  for(t in 2:T) h[t] = mu + phi_h*(h[t-1] - mu) + s_h*rnorm(1)
  mu_t = mean_st + omega*delta*W*exp(0.5*h)/sqrt(U)
  sigma_t = omega*sqrt(1-(delta)^2)*exp(0.5*h)/sqrt(U)
  y = mu_t + sigma_t*rnorm(T)
  #####
  
  ### Sampling2
  test=1
  while(any(test>0, is.nan(test), is.na(test))){
    draws = rstan::sampling(model_stan2, 
                            data = list(T = length(y), 
                                        y = as.numeric(y),
                                        lambda = -log(0.5)/0.5),
                            chains = 1,
                            warmup = warmup,
                            iter = warmup + iters,
                            cores = 1
    )
    x = rstan::extract(draws, pars = c('mu','phi_h','s_h','s_a','h','v'))
    theta = rbind(x$mu, x$phi_h, x$s_h, x$s_a, x$v)
    s = num_analisys(draws = theta, 
                     names = c('mu','phi_h','s_h','s_a','v'),
                     digits = 4, hdp = TRUE)
    test = sum(abs(s[,'CD']) > 1.96)
  }
  summary2 = s
  
  
  time = Sys.time() - time
  list(#summary1 = summary1,
       summary2 = summary2,
       #summary3 = summary3,
       time = time)
  
}

stopCluster(cl)

save(result, theta_vdd, file = paste0('pcp_sa_', s_a, '.RData'))

