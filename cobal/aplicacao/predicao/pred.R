### VaR and LPDS*
source( 'https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/svm_VaR.R' )
source( 'https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/svm_lpds_star.R' )
source( 'https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/tvpsvm_VaR.R' )
source( 'https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/cobal/source/tvp_lpds_star.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/refs/heads/main/source/num_analisys.R' )

# compiling Stan model
dir = paste0(getwd(), '/cobal/')
path = paste0(dir, 'models/static.stan')
model_stan1 = rstan::stan_model(file = path)
path = paste0(dir, 'models/ig.stan')
model_stan2 = rstan::stan_model(file = path)
path = paste0(dir, 'models/pcp.stan')
model_stan3 = rstan::stan_model(file = path)
path = paste0(dir, 'models/jeffrey.stan')
model_stan4 = rstan::stan_model(file = path)
path = paste0(dir, 'models/exp.stan')
model_stan5 = rstan::stan_model(file = path)

# outs directory
dir_out = paste0( dir, 'aplication/pred/ouro/' )

horizon = 120
VaR_hat1 = VaR_hat2 = VaR_hat3 = VaR_hat4 = VaR_hat5 = rep(0, horizon)
lpds.star1 = lpds.star2 = lpds.star3 = lpds.star4 = lpds.star5 = rep(0, horizon)
summary1 = summary2 = summary3 = summary4 = summary5 = list()
Time = 0

warmup = 2e3
iters = 1e3

for(it in 1:horizon){
  
  time = Sys.time()
  
  y.now = log.ret[ 1:(T - horizon + it - 1) ]
  
  ##########
  # Model1 #
  ##########
  cat( '\014' )
  cat('Horizon: ', it, '/', horizon, '\n', '\n')
  if( it == 1 ){
    cat('Estimated total remaining time: Calculating...', '\n', '\n' )
  }else{
    cat('Estimated total remaining time: ', 
        round((horizon - it + 1) * (Time / (it-1)), 1), 'hs',
        '\n', '\n' )
  }
  cat( 'Model 1', '\n' )
  draws = rstan::sampling(model_stan1, 
                          data = list(T = length( y.now ), 
                                      y = as.numeric( y.now )
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s2_h, x$a)
  h_hmc = x$h
  # Summary
  summary1[[ it ]] = num_analisys(draws = theta, 
                                  names = c('b', 
                                            'mu', 
                                            'phi_h', 
                                            's2_h', 
                                            'a'),
                                  digits = 4,
                                  hdp = TRUE
  )
  # Evaluation VaR
  VaR_hat1[ it ] = svm_VaR(h_T = h_hmc[, dim( h_hmc )[2]],
                           theta_hmc = theta
  )
  # Evaluation LPDS
  lpds.star1[ it ] = svm_lpds_star(h_T = h_hmc[, dim( h_hmc )[2]],
                                   theta_hmc = theta,
                                   yobs = log.ret[T - horizon + it]
  )
  
  ##########
  # Model2 #
  ##########
  cat( '\014' )
  cat('Horizon: ', it, '/', horizon, '\n', '\n')
  if( it == 1 ){
    cat('Estimated total remaining time: Calculating...', '\n', '\n' )
  }else{
    cat('Estimated total remaining time: ', 
        round((horizon - it + 1) * (Time / (it-1)), 1), 'hs',
        '\n', '\n' )
  }
  cat( 'Model 2', '\n' )
  draws = rstan::sampling(model_stan2, 
                          data = list(T = length( y.now ), 
                                      y = as.numeric( y.now )
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s2_h, x$s2_a)
  h_hmc = x$h
  a_hmc = x$a
  # Summary
  summary2[[ it ]] = num_analisys(draws = theta, 
                                  names = c('b', 
                                            'mu', 
                                            'phi_h', 
                                            's2_h', 
                                            'ls_a'),
                                  digits = 4,
                                  hdp = TRUE
  )
  # Evaluation VaR
  VaR_hat2[ it ] = tvp_VaR(h_T = h_hmc[, dim( h_hmc )[2]],
                           a_T = a_hmc[, dim( h_hmc )[2]],
                           theta_hmc = theta
  )
  # Evaluation LPDS
  lpds.star2[ it ] = tvp_lpds_star(h_T = h_hmc[, dim( h_hmc )[2]],
                                   a_T = a_hmc[, dim( h_hmc )[2]],
                                   theta_hmc = theta,
                                   yobs = log.ret[T - horizon + it]
  )
  ##########
  # Model3 #
  ##########
  cat( '\014' )
  cat('Horizon: ', it, '/', horizon, '\n', '\n')
  if( it == 1 ){
    cat('Estimated total remaining time: Calculating...', '\n', '\n' )
  }else{
    cat('Estimated total remaining time: ', 
        round((horizon - it + 1) * (Time / (it-1)), 1), 'hs',
        '\n', '\n' )
  }
  cat( 'Model 3', '\n' )
  draws = rstan::sampling(model_stan3, 
                          data = list(T = length( y.now ), 
                                      y = as.numeric( y.now ),
                                      lambda = -log(0.5) / sqrt(0.5)
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s2_h, x$s2_a)
  h_hmc = x$h
  a_hmc = x$a
  # Summary
  summary3[[ it ]] = num_analisys(draws = theta, 
                                  names = c('b', 
                                            'mu', 
                                            'phi_h', 
                                            's2_h', 
                                            'ls_a'),
                                  digits = 4,
                                  hdp = TRUE
  )
  # Evaluation VaR
  VaR_hat3[ it ] = tvp_VaR(h_T = h_hmc[, dim( h_hmc )[2]],
                           a_T = a_hmc[, dim( h_hmc )[2]],
                           theta_hmc = theta
  )
  # Evaluation LPDS
  lpds.star3[ it ] = tvp_lpds_star(h_T = h_hmc[, dim( h_hmc )[2]],
                                   a_T = a_hmc[, dim( h_hmc )[2]],
                                   theta_hmc = theta,
                                   yobs = log.ret[T - horizon + it]
  )
  ##########
  # Model4 #
  ##########
  cat( '\014' )
  cat('Horizon: ', it, '/', horizon, '\n', '\n')
  if( it == 1 ){
    cat('Estimated total remaining time: Calculating...', '\n', '\n' )
  }else{
    cat('Estimated total remaining time: ', 
        round((horizon - it + 1) * (Time / (it-1)), 1), 'hs',
        '\n', '\n' )
  }
  cat( 'Model 4', '\n' )
  draws = rstan::sampling(model_stan4, 
                          data = list(T = length( y.now ), 
                                      y = as.numeric( y.now )
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s2_h, x$s2_a)
  h_hmc = x$h
  a_hmc = x$a
  # Summary
  summary4[[ it ]] = num_analisys(draws = theta, 
                                  names = c('b', 
                                            'mu', 
                                            'phi_h', 
                                            's2_h', 
                                            'ls_a'),
                                  digits = 4,
                                  hdp = TRUE
  )
  # Evaluation VaR
  VaR_hat4[ it ] = tvp_VaR(h_T = h_hmc[, dim( h_hmc )[2]],
                           a_T = a_hmc[, dim( h_hmc )[2]],
                           theta_hmc = theta
  )
  # Evaluation LPDS
  lpds.star4[ it ] = tvp_lpds_star(h_T = h_hmc[, dim( h_hmc )[2]],
                                   a_T = a_hmc[, dim( h_hmc )[2]],
                                   theta_hmc = theta,
                                   yobs = log.ret[T - horizon + it]
  )
  ##########
  # Model5 #
  ##########
  cat( '\014' )
  cat('Horizon: ', it, '/', horizon, '\n', '\n')
  if( it == 1 ){
    cat('Estimated total remaining time: Calculating...', '\n', '\n' )
  }else{
    cat('Estimated total remaining time: ', 
        round((horizon - it + 1) * (Time / (it-1)), 1), 'hs',
        '\n', '\n' )
  }
  cat( 'Model 5', '\n' )
  draws = rstan::sampling(model_stan5, 
                          data = list(T = length( y.now ), 
                                      y = as.numeric( y.now )
                          ),
                          chains = 1,
                          warmup = warmup,
                          iter = warmup + iters,
                          cores = 1
  )
  x = rstan::extract(draws)
  theta = matrix(x$b, nrow = 1)
  theta = rbind(theta, x$mu, x$phi, x$s2_h, x$s2_a)
  h_hmc = x$h
  a_hmc = x$a
  # Summary
  summary5[[ it ]] = num_analisys(draws = theta, 
                                  names = c('b', 
                                            'mu', 
                                            'phi_h', 
                                            's2_h', 
                                            'ls_a'),
                                  digits = 4,
                                  hdp = TRUE
  )
  # Evaluation VaR
  VaR_hat5[ it ] = tvp_VaR(h_T = h_hmc[, dim( h_hmc )[2]],
                           a_T = a_hmc[, dim( h_hmc )[2]],
                           theta_hmc = theta
  )
  # Evaluation LPDS
  lpds.star5[ it ] = tvp_lpds_star(h_T = h_hmc[, dim( h_hmc )[2]],
                                   a_T = a_hmc[, dim( h_hmc )[2]],
                                   theta_hmc = theta,
                                   yobs = log.ret[T - horizon + it]
  )
  
  time = Sys.time() - time
  Time = Time + as.numeric(time, units = 'hours')
  if(it == horizon){
    save(VaR_hat1, lpds.star1, summary1, 
         file = paste0( dir_out, 'static.RData' )
         )
    save(VaR_hat2, lpds.star2, summary2, 
         file = paste0( dir_out, 'ig.RData' )
         )
    save(VaR_hat3, lpds.star3, summary3, 
         file = paste0( dir_out, 'pcp.RData' )
         )
    save(VaR_hat4, lpds.star4, summary4, 
         file = paste0( dir_out, 'jeffrey.RData' )
         )
    save(VaR_hat5, lpds.star5, summary5, 
         file = paste0( dir_out, 'exp.RData' )
    )
  } 
  
}

