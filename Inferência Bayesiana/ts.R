# Simulando dados
tau = 0.1
mu = -5
phi = 0.985
s2_h = 0.025
T = 3e3

# time-series 1
s2_a1 = 0.0
y1 = h = a = matrix(0, nrow = T, ncol = 1)
a[1] = rnorm(1, mean = -2, sd = 1)
h[1] = mu + sqrt( s2_h / (1 - phi * phi) ) * rnorm( 1 )
y1[1] = tau + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )
for( t in 2:T ){
  a[t] = a[t-1] + sqrt( s2_a1 ) * rnorm( 1 )
  h[t] = mu + phi * (h[t-1] - mu) + sqrt( s2_h ) * rnorm( 1 )
  y1[t] = tau + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
}

# time-series 2
s2_a2 = 0.01
y2 = h = a = matrix(0, nrow = T, ncol = 1)
a[1] = rnorm(1, mean = -2, sd = 1)
h[1] = mu + sqrt( s2_h / (1 - phi * phi) ) * rnorm( 1 )
y2[1] = tau + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )
for( t in 2:T ){
  a[t] = a[t-1] + sqrt( s2_a2 ) * rnorm( 1 )
  h[t] = mu + phi * (h[t-1] - mu) + sqrt( s2_h ) * rnorm( 1 )
  y2[t] = tau + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
}

# time-series 3
s2_a3 = 0.05
y3 = h = a = matrix(0, nrow = T, ncol = 1)
a[1] = rnorm(1, mean = -2, sd = 1)
h[1] = mu + sqrt( s2_h / (1 - phi * phi) ) * rnorm( 1 )
y3[1] = tau + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )
for( t in 2:T ){
  a[t] = a[t-1] + sqrt( s2_a3 ) * rnorm( 1 )
  h[t] = mu + phi * (h[t-1] - mu) + sqrt( s2_h ) * rnorm( 1 )
  y3[t] = tau + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
}

# time-series 4
s2_a4 = 0.1
y4 = h = a = matrix(0, nrow = T, ncol = 1)
a[1] = rnorm(1, mean = -2, sd = 1)
h[1] = mu + sqrt( s2_h / (1 - phi * phi) ) * rnorm( 1 )
y1[1] = tau + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )
for( t in 2:T ){
  a[t] = a[t-1] + sqrt( s2_a4 ) * rnorm( 1 )
  h[t] = mu + phi * (h[t-1] - mu) + sqrt( s2_h ) * rnorm( 1 )
  y4[t] = tau + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
}

# time-series 5
s2_a5 = 0.5
y5 = h = a = matrix(0, nrow = T, ncol = 1)
a[1] = rnorm(1, mean = -2, sd = 1)
h[1] = mu + sqrt( s2_h / (1 - phi * phi) ) * rnorm( 1 )
y1[1] = tau + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )
for( t in 2:T ){
  a[t] = a[t-1] + sqrt( s2_a5 ) * rnorm( 1 )
  h[t] = mu + phi * (h[t-1] - mu) + sqrt( s2_h ) * rnorm( 1 )
  y5[t] = tau + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
}

# time-series 6
s2_a6 = 1.0
y6 = h = a = matrix(0, nrow = T, ncol = 1)
a[1] = rnorm(1, mean = -2, sd = 1)
h[1] = mu + sqrt( s2_h / (1 - phi * phi) ) * rnorm( 1 )
y1[1] = tau + a[1] * exp( h[1] ) + exp( 0.5 * h[1] ) * rnorm( 1 )
for( t in 2:T ){
  a[t] = a[t-1] + sqrt( s2_a6 ) * rnorm( 1 )
  h[t] = mu + phi * (h[t-1] - mu) + sqrt( s2_h ) * rnorm( 1 )
  y6[t] = tau + a[t] * exp( h[t] ) + exp( 0.5 * h[t] ) * rnorm( 1 )
}
par( mfrow = c(3, 2) )
plot(y1, type = 'l', xlab = '(a)', ylab = '')
plot(y2, type = 'l', xlab = '(b)', ylab = '')
plot(y3, type = 'l', xlab = '(c)', ylab = '')
plot(y4, type = 'l', xlab = '(d)', ylab = '')
plot(y5, type = 'l', xlab = '(e)', ylab = '')
plot(y6, type = 'l', xlab = '(f)', ylab = '')
par( mfrow = c(1, 1) )
