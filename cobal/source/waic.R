############### waic
waic = function(data, mut, st){
  n = length(data)
  M = matrix(nrow = dim(mut)[1], ncol = n)
  for(j in 1:n){
    M[, j] = dnorm(y[j], 
                   mean = mut[, j], 
                   sd = st[, j], 
                   log = TRUE)
  }
  return(loo::waic(M))
}
