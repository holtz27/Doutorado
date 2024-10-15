### VaR and LPDS*
source( 'https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/tvpsvm_VaR.R' )
source( 'https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/tvp_lpds_star.R' )
source('https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R')

# compiling Stan model
dir = paste0(getwd(), '/Doutorado/Eventos/cobal/modelos/')
path = paste0(dir, 'exp.stan')
model_stan = rstan::stan_model(file = path)

horizon = 1
VaR_hat = lpds.star = rep(0, horizon)
summary = list()

for(it in 1:horizon){
  
  if(it == 1) time = Sys.time()
  
  y.now = log.ret[ 1:(T - horizon + it - 1) ]
  draws = rstan::sampling(model_stan, 
                          data = list(T = length( y.now ), 
                                      y = as.numeric( y.now )
                          ),
                          chains = 1,
                          warmup = 2e3,
                          iter = 2e3 + 1e3,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s2_h, x$s2_a)
  h_hmc = x$h
  a_hmc = x$a
  # Summary
  summary[[ it ]] = num_analisys(draws = theta, 
                                 names = c('b', 
                                           'mu', 
                                           'phi_h', 
                                           's2_h', 
                                           's2_a'),
                                 digits = 4,
                                 hdp = TRUE
  )
  # Evaluation VaR
  VaR_hat[ it ] = VaR(h_T = h_hmc[, dim( h_hmc )[2]],
                      a_T = a_hmc[, dim( h_hmc )[2]],
                      theta_hmc = theta
  )
  # Evaluation LPDS
  lpds.star[ it ] = tvp_lpds_star(h_T = h_hmc[, dim( h_hmc )[2]],
                                  a_T = a_hmc[, dim( h_hmc )[2]],
                                  theta_hmc = theta,
                                  yobs = log.ret[T - horizon + it]
  )
  
  if(it == horizon) time = Sys.time() - time
  
}

time

save(VaR_hat, lpds.star, summary, file = 'pcp_out_sample.RData')
