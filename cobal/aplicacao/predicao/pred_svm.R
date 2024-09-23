### VaR and LPDS*
source( 'https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/svm_VaR.R' )
source( 'https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/svm_lpds_star.R' )
# compiling Stan model
model_stan = rstan::stan_model(file = 'svm.stan')

horizon = 2
VaR_hat = lpds.star = rep(0, horizon)

for(it in 1:horizon){
  
  if(it == 1) time = Sys.time()
  
  y.now = log.ret[ 1:(T - horizon + it - 1) ]
  draws = rstan::sampling(model_stan, 
                          data = list(T = length( y.now ), 
                                      y = as.numeric( y.now )
                                      ),
                          chains = 2,
                          warmup = 2e2,
                          iter = 2e2 + 500,
                          cores = 2
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
  h_hmc = x$h
  # Evaluation VaR
  VaR_hat[ it ] = svm_VaR(h_T = h_hmc[, dim( h_hmc )[2]],
                          theta_hmc = theta
  )
  # Evaluation LPDS
  lpds.star[ it ] = svm_lpds_star(h_T = h_hmc[, dim( h_hmc )[2]],
                                  theta_hmc = theta,
                                  yobs = log.ret[T - horizon + it]
  )
  
  if(it == horizon) time = Sys.time() - time
  
}

time
Data = list( VaR_hat, lpds.star )
save(Data, file = 'svm_out_sample.RData')
