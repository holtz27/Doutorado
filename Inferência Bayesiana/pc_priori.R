### Random Walk
T = 500
mu = matrix(0, nrow = T, ncol = 1)
Sigma = matrix(nrow = T, ncol = T)
s_a = 0.25

for(i in 1:T){
  Sigma[i, seq(1, i)] = seq(1, i)
  Sigma[seq(1, i), i] = Sigma[i, seq(1, i)]
}

#set.seed( 42 )
rw = mvtnorm::rmvnorm(n = 1, mean = mu, sigma = s_a * s_a * Sigma)
plot(1:T, rw, type = 'l')


### PC Priori
pcpriori = function(x, lambda){
  y = 0.5 * lambda * x**(-0.5) * exp( - lambda * sqrt( x ) )
  return( y )
}
a = 0.01
U = 0.2
lambda = -log( a ) / sqrt( U ) 

x = seq(0, 1, by = 0.01)
plot(x, 
     dweibull(x, shape = 0.5, scale = lambda**(-2)), 
     type = 'l', 
     xlab = expression(sigma^{2}), ylab = 'densidade',
     main = 'PC priori'
)
