loglik_APF = function(y, theta, N){
  
  moms = function(a, h, v){
    
    invsqrt_U=sqrt(0.5*v)*gamma(0.5*(v-1))/gamma(0.5*v)
    #U=rgamma(length(a), shape=0.5*v, rate=0.5*v)
    #invsqrt_U=1/sqrt(U)
    W=sqrt(2/pi)
    #W=rtnorm(length(a))
    
    delta_t = a / sqrt(1 + a^2)
    k1=sqrt(0.5*v)*gamma(0.5*(v-1))/gamma(0.5*v)
    k2=v/(v-2)
    omega_t =1/sqrt(k2-(2/pi)*(delta_t*k1)^2)
    gamma_t = -sqrt(2/pi)*delta_t*omega_t*k1
    mu_t = gamma_t + omega_t * delta_t * W * exp(0.5*h) 
    mu_t = mu_t * invsqrt_U
    sigma_t = omega_t*sqrt(1-delta_t^2)*exp(0.5*h)
    sigma_t = sigma_t * invsqrt_U
    
    return(matrix(c(mu_t, sigma_t), ncol=2))
  }
  
  T = length(y)
  mu = theta['mu']
  phi = theta['phi']
  sh = theta['sh']
  sa = theta['sa']
  v = theta['v']
  
  h_par = a_par = w = matrix(0, nrow = N, ncol = T)
  
  h_par0 = rnorm(N, mean=mu, sd=sh/sqrt(1-phi^2))
  a_par0 = rep(0, N)
  w0 = rep(1/N, N)
  hbar=abar=Pk=indx=numeric(N)
  
  # t=1
  # mu_t
  hbar = mu + phi*(h_par0 - mu)
  abar = a_par0
  Ms = moms(a=abar, h=hbar, v=v)
  Pk = dnorm(y[1], mean=Ms[,1], sd=Ms[,2]) * w0
  Pk = Pk/sum(Pk)
  indx = sample(1:N, N, prob=Pk, replace=TRUE)
  
  # new states
  h_par[,1] = rnorm(N, mu + phi*(hbar[indx] - mu), sh)
  a_par[,1] = rnorm(N, abar[indx], sa)
  
  Ms=moms(a=a_par[,1], h=h_par[,1], v=v)
  dens1=dnorm(y[1], mean=Ms[,1], sd=Ms[,2])
  
  Ms=moms(a=abar, h=hbar, v=v)
  dens2=dnorm(y[1], mean=Ms[,1], sd=Ms[,2])
  
  w[,1]=dens1/dens2
  w[,1]=w[,1]/sum(w[,1])
  
  loglik = log( sum(dens1*w[,1]) )
  
  # t>1
  for(t in 2:T){
    
    # mu_t
    hbar = mu + phi*(h_par[,t-1] - mu)
    abar = a_par[,t-1]
    Ms=moms(a=abar, h=hbar, v=v)
    Pk = dnorm(y[t], mean=Ms[,1], sd=Ms[,2]) * w[,t-1]
    Pk = Pk/sum(Pk)
    indx=sample(1:N, N, prob=Pk, replace=TRUE)
    
    # new states
    h_par[,t] = rnorm(N, mu+phi*(hbar[indx]-mu), sh)
    a_par[,t] = rnorm(N, abar[indx], sa)
    
    Ms=moms(a=a_par[,t], h=h_par[,t], v=v)
    dens1=dnorm(y[t], mean=Ms[,1], sd=Ms[,2])
    
    Ms=moms(a=abar, h=hbar, v=v)
    dens2=dnorm(y[t], mean=Ms[,1], sd=Ms[,2])
    
    w[,t]=dens1/dens2
    w[,t]=w[,t]/sum(w[,t])
    
    loglik = loglik + log(sum( dens1*w[,t] ))
    
    hbar=abar=Pk=indx=numeric(N)
  }
  
  return(list(loglik=loglik, pred_h=h_par, pred_a=a_par, w=w)) 
}
quant = function(x, weights, probs=c(0.025, 0.5, 0.975)){
  
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
Rcpp::sourceCpp('~/topicos/st/loglik_APF.cpp')
################################################################################
################################################################################


theta=c(mu=mu, phi=phi_h, sh=s_h, sa=0.1, v=v)
X=loglik_APF(y=y, theta=theta, N=1e4)
X$loglik

# Ess
weigth=X$w
plot(1/colSums( weigth^2 ))

# Plots pred_a, pred_h
a_hat=a_min=a_max=numeric(T)
h_hat=h_min=h_max=numeric(T)
for(t in 1:T){
  h_hat[t] = sum(X$pred_h[,t]*X$w[, t])
  a_min[t] = quant(X$pred_a[,t], weights=X$w[,t], probs=0.025)
  a_hat[t] = sum(X$pred_a[,t]*X$w[,t])
  a_max[t] = quant(X$pred_a[,t], weights=X$w[,t], probs=0.975)
}

# Pred_h
plot(h_hat, type='l', col='purple', lwd=1)
lines(h, col='gray60', lwd=2)

# Pred_a
plot(a, type='l', col='gray60', ylim=c(-1, 1), lwd=2)
lines(a_hat, col='purple')
lines(a_min, col='purple', lty=2)
lines(a_max, col='purple', lty=2)
polygon(c(1:T, rev(1:T)),               
        c(a_max, rev(a_min)),       
        col = rgb(0.6, 0.2, 0.8, 0.2), 
        border = NA)                


################################################################################
################################################################################
sax=seq(0, 1, 0.005)
logL=numeric(length(sax))

for(i in 1:length(sax)){
  theta=c(mu=mu, phi=phi_h, sh=s_h, sa=sax[i], v=v)
  X=loglik_APF(y=y, theta=theta, N=5e3)
  logL[i]=X$loglik  
}

plot(sax, logL)

























