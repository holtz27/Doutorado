############### waic
waic = function(data, mut, st){
  n = length(data)
  M = matrix(nrow = dim(mut)[1], ncol = n)
  for(j in 1:n){
    M[, j] = dnorm(y[j], 
                   mean = mut[, 1], 
                   sd = st[, 1], 
                   log = TRUE)
  }
  return(loo::waic(M))
}
