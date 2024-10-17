loglik = function(data_i, draws, data_, T){
  # data_ = ( y_{1}, y_{2}, ..., y_{T} )
  # data_past = ( y_{0}, y_{1}, ..., y_{T-1} )
  # draws (2 * T + 3) x M
  k = which( data_ == as.numeric( data_i ) )
  N = ncol( draws )
  log_l = rep(0, N)
  
  for(col in 1:N){
    
    b = as.numeric( draws[1, col] )
    h = as.numeric( draws[1 + k, col] )
    a = as.numeric( draws[T + 1 + k, col] )
    
    log_l[ col ] = dnorm(as.numeric( data_i ), 
                         mean = b + a * exp( h ), 
                         sd = exp( 0.5 * h ), 
                         log = TRUE)
  }
  
  return( log_l )
}

############### waic
waic = function(data, draws){
  T = length( data )
  X = sapply(X = data, 
             FUN = loglik,
             draws,
             data,
             length( data )
  ) 
  
  waic = loo::waic( X )
  
  return( waic )
}
