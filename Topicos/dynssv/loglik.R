quant=function(x, weights, probs=c(0.025, 0.5, 0.975)){
  
  ix = order(x)
  x_sorted = x[ix]
  weights_sorted = weights[ix]
  weights_sorted = weights_sorted / sum(weights_sorted)
  cumulative_weights = cumsum(weights_sorted)
  quants = numeric(length(probs))
  for(i in seq_along(probs)){
    index = which(cumulative_weights >= probs[i])[1]
    quants[i] = x_sorted[index]
  }
  
  return(quants)
}
Rcpp::sourceCpp('~/Documentos/loglik_dynsn.cpp')
Rcpp::sourceCpp('~/Documentos/loglik_staticsn.cpp')
Rcpp::sourceCpp('~/Documentos/loglik_dynst.cpp')
Rcpp::sourceCpp('~/Documentos/loglik_staticst.cpp')
Rcpp::sourceCpp('~/Documentos/loglik_dynss.cpp')
Rcpp::sourceCpp('~/Documentos/loglik_staticss.cpp')
################################################################################
################################################################################

# static
# sn
theta=c(mu=2.3043, phi=0.8807, sh=0.5212, a1=0.0203)
X=loglik_staticsn(y=log.ret, theta=theta, N=1e4)
-X$loglik/T

# st
theta=c(mu=2.5604, phi=0.9761, sh=0.1991, a1=0.0159, v=4.1170)
X=loglik_staticst(y=log.ret, theta=theta, N=1e4)
-X$loglik/T
# ss
theta=c(mu=2.8564, phi=0.9769, sh=0.2040, a1=0.0148, v=1.3277)
X=loglik_staticss(y=log.ret, theta=theta, N=1e4)
-X$loglik/T

# dynamic
# sn
theta=c(mu=2.2868, phi=0.8837, sh=0.5114, a1=0.1331, sa=0.0176)
X=loglik_dynsn(y=log.ret, theta=theta, N=1e4)
-X$loglik/T

# st
theta=c(mu=2.6063, phi=0.9773, sh=0.1937, a1=0.1771, sa=0.0189, v=3.7915)
X=loglik_dynst(y=log.ret, theta=theta, N=1e4)
-X$loglik/T

#ss
theta=c(mu=2.8250, phi=0.9771, sh=0.2015, a1=0.2118, sa=0.0186, v=1.3371)
X=loglik_dynss(y=log.ret, theta=theta, N=1e4)
-X$loglik/T






















# Ess
#weigth=X$w
#plot(1/colSums( weigth^2 ))

# Plots pred_a, pred_h
a_hat=a_min=a_max=numeric(T)
h_hat=h_min=h_max=numeric(T)
for(t in 1:T){
  h_hat[t] = sum(X$pred_h[,t]*X$w[,t])
  a_min[t] = quant(X$pred_a[,t], weights=X$w[,t], probs=0.025)
  a_hat[t] = sum(X$pred_a[,t]*X$w[,t])
  a_max[t] = quant(X$pred_a[,t], weights=X$w[,t], probs=0.975)
}

# Pred_h
plot(abs(log.ret), type='l', col='gray60', lwd=1)
lines(exp(0.5*h_hat))

# Pred_a
plot(a_hat, type='l', col='gray60', ylim=c(-1, 1), lwd=2)
#lines(a_hat, col='purple')
lines(a_min, col='purple', lty=2)
lines(a_max, col='purple', lty=2)
abline(h=0, lty=2)
polygon(c(1:T, rev(1:T)),               
        c(a_max, rev(a_min)),       
        col = rgb(0.6, 0.2, 0.8, 0.2), 
        border = NA)                


################################################################################
################################################################################
sax=seq(-1, 1, 0.05)
logL=numeric(length(sax))

for(i in 1:length(sax)){
  theta=c(mu=sax[i], phi=phi_h, sh=s_h, sa=s_a, v=v)
  X=loglik_APF(y=y, theta=theta, N=1e4)
  logL[i]=X$loglik  
}

plot(sax, logL)

























