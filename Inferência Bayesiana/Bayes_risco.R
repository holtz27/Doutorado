# Definir a função r
r = function(U, p) U**2 / log(p)**4

a = 0.5 * 9
b = 0.5 * 0.13
r_ig = b**2 * gamma( a - 2) / gamma( a )



x = seq(0, 0.25, 0.005)


plot(x, r(0.1, x), 
     main = 'Risco de Bayes PC priori', xlab = 'p', ylab = 'risco',
     type = 'l', lty = 1)
lines(x, r(0.25, x), lty = 2)
lines(x, r(0.5, x), lty = 3)
abline(h = r_ig, col = '#A9A9A9', lty = 2)

legend('topleft', legend = c('U = 0.1', 
                             'U = 0.25', 
                             'U = 0.5'), lty = 1:3, bty = "n")
