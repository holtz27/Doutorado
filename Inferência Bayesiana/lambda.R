### PC Priori
x = seq(0, 0.25, by = 0.005)

# high informative
p1 = 0.01
U1 = 0.1
lambda1 = -log( p1 ) / sqrt( U1 ) 
y1 = dweibull(x, shape = 0.5, scale = lambda1**(-2))

#middle informative
p2 = 0.1
U2 = 0.1
lambda2 = -log( p2 ) / sqrt( U2 ) 
y2 = dweibull(x, shape = 0.5, scale = lambda2**(-2))

# weakly informative
p3 = 0.5
U3 = 0.1
lambda3 = -log( p3 ) / sqrt( U3 ) 
y3 = dweibull(x, shape = 0.5, scale = lambda3**(-2))

plot(x, 
     y1, 
     type = 'l',
     ylim = c(0, 15),
     xlab = expression(sigma[a]^{2}), ylab = 'densidade',
     main = 'PC priori'
)
lines(x, y2, lty = 2)
lines(x, y3, lty = 3)
legend('topright', 
       legend = c(paste0('(U, p) = (', U1, ', ', p1,')' ), 
                  paste0('(U, p) = (', U2, ', ', p2,')' ), 
                  paste0('(U, p) = (', U3, ', ', p3,')' )), 
      lty = c(1, 2, 3), lwd = 2, bty = 'n')
