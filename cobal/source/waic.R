############### waic
waic = function(data, b, a, h, dyn=TRUE){
  n = length(data)
  M = matrix(nrow = dim(h)[1], ncol = n)
  if(dyn){
    for(j in 1:n){
      M[, j] = dnorm(data[j], 
                     mean = b + a[,j]*exp(h[,j]), 
                     sd = exp(0.5*h[,j]), 
                     log = TRUE)
    }
  }else{
    for(j in 1:n){
      M[, j] = dnorm(data[j], 
                     mean = b + a*exp(h[,j]), 
                     sd = exp(0.5*h[,j]), 
                     log = TRUE)
    }
  }
  return(loo::waic(M))
}
