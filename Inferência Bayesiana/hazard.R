h = function(t, lambda) 0.5 * lambda / sqrt( t )

### h(t) hazard function
t = seq(0, 0.50, by = 0.001)
plot( t, 
      h(t, lambda1), 
      type = 'l',
      xlab = expression(xi),
      ylab = 'h(t)',
      ylim = c(0, 50),
      main = 'Função de Hazard'
      )
lines(t, h(t, lambda2), lty = 2)
lines(t, h(t, lambda3), lty = 3)
legend('topright',
       inset = c(0.05, 0),
       legend = c(paste0('p = ', p1), 
                  paste0('p = ', p2), 
                  paste0('p = ', p3)
                  ), 
       lty = c(1, 2, 3), lwd = 2, bty = 'n')

