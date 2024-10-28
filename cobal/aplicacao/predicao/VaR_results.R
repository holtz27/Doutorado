load("~/Doutorado/Eventos/cobal/aplicacao/prediction/static_out.RData")
lpds.star1 = lpds.star
VaR_hat1 = VaR_hat
var_out1 = sum( tail(log.ret, 100) < VaR_hat1 )



load("~/Doutorado/Eventos/cobal/aplicacao/prediction/ig_out.RData")
lpds.star2 = lpds.star
VaR_hat2 = VaR_hat
var_out2 = sum( tail(log.ret, 100) < VaR_hat2 )

 
sum( -lpds.star1 )
sum( -lpds.star2 )
var_out1
var_out2



plot( tail(log.ret, 100), 
      main = '95% value at risk (VaR) of gold returns', 
      xlab = '', ylab = 'Returns',
      ylim = c(-2, 3) )
lines( VaR_hat1 )
lines( VaR_hat2, lty = 2 )
legend('topright', 
       legend = c('Static', 'IG'), 
       col = c('black', 'black'), 
       lty = c(1, 2), 
       bty = 'n')
