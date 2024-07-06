library(invgamma)

x = seq(0.0, 0.1, 0.0005)

y_c = dinvgamma(x, shape = 4.5, rate = 0.065)
y_j = 1/x
lambda = -log(0.25) / sqrt(.1)
y_pc = dweibull(x, shape = 0.5, scale = 1 / lambda**2)

plot(x, y_j, 
     type = 'l',
     xlab = expression(xi),
     ylab = 'Densidade',
     ylim = c(0, 100),
     lty = 3)  # tipo de linha sólida
lines(x, y_c, lty = 2)  # tipo de linha tracejada
lines(x, y_pc, lty = 1)  # tipo de linha pontilhada

# Adicionar legenda sem contorno
legend("topright", 
       legend = c("Jeffrey", "Inv. Gamma", "Pc prior"), 
       lty = c(3, 2, 1),
       bty = "n")
