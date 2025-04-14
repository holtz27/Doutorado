rtnorm = function(n){
  u = runif(n)
  return(qnorm(0.5*(u+1)))
}

T = 2e1
# log-volatility
mu = 0.1
phi_h = 0.98
s_h = 0.22
# dynskew
s_a = 0.01
v = 8
#Data
#set.seed(42)
a = delta = omega = numeric(T)
a[1] = 0
for(t in 2:T) a[t] = a[t-1] + s_a*rnorm(1)
delta = a/sqrt(1+a*a)
k1 = sqrt(0.5*v)*gamma(0.5*(v-1))/gamma(0.5*v)
k2 = v/(v-2)
omega = 1/sqrt(k2-(2/pi)*(delta*k1)^2)
gamma_t = -sqrt(2/pi)*delta*omega*k1
W = rtnorm(T)
U = rgamma(T, shape=0.5*v, rate=0.5*v)
y = h = numeric(T)
h[1] = mu + s_h/sqrt(1-phi_h*phi_h)*rnorm(1)
for(t in 2:T) h[t] = mu + phi_h*(h[t-1]-mu) + s_h*rnorm(1)
mu_t = gamma_t + omega*delta*W*exp(0.5*h)/sqrt(U)
sigma_t = omega*sqrt(1-delta^2)*exp(0.5*h)/sqrt(U)
y = mu_t + sigma_t*rnorm(T)
par(mfrow=c(1,2))
plot(y, type='l', xlab='day', ylab='Return')
hist(y, breaks=40, xlab='', main='', freq=FALSE)
par(mfrow=c(1,1))
################################################################################
# W
postW = function(y, h, U, a, v){
  
  delta = a/sqrt(1+a*a)
  k1 = sqrt(0.5*v)*gamma(0.5*(v-1))/gamma(0.5*v)
  k2 = v/(v-2)
  omega = 1/sqrt(k2-(2/pi)*(delta*k1)^2)
  gamma_t = -sqrt(2/pi)*delta*omega*k1
  
  mu_t = delta*sqrt(U)*exp(-0.5*h)*(y-gamma_t*exp(0.5*h))/omega
  sigma_t = 1/sqrt(1-delta*delta)
  
  n = length(mu_t)
  nt = numeric(n)
  
  for(i in 1:n){
    z = runif(1)*(1-pnorm(0, mu_t[i], sigma_t[i]))+pnorm(0, mu_t[i], sigma_t[i])
    nt[i] = qnorm(z, mu_t[i], sigma_t[i])
  }
  
  return(nt)
}
z = postW(y=y, h=h, U=U, a=a, v=v)
plot(z)
# U
postU = function(U, y, h, W, a, v, tx){
  
  n = length(y)
  U_star = logqU_star = logqU = numeric(n)
  U_end = U
  
  delta = a/sqrt(1+a*a)
  k1 = sqrt(0.5*v)*gamma(0.5*(v-1))/gamma(0.5*v)
  k2 = v/(v-2)
  omega = 1/sqrt(k2-(2/pi)*(delta*k1)^2)
  gamma_t = -sqrt(2/pi)*delta*omega*k1
  
  c = (exp(-0.5*h)*y-gamma_t)/(omega*(1-delta*delta))
  logQ = sqrt(U)*delta*W*c
  logpU = logQ + 0.5*(v-1)*log(U) - 0.5*U*(v+(1-delta*delta)*c)
  
  for(i in 1:n){
    
    U_star[i] = rgamma(1, 
                       shape=0.5*(v + 1), 
                       rate=0.5*(v + (1+delta[i]*delta[i])*c[i]^2) )
    
    logqU_star[i] = dgamma(U_star[i], 
                           shape=0.5*(v + 1), 
                           rate=0.5*(v + (1+delta[i]*delta[i])*c[i]^2), 
                           log=TRUE)
    
    logqU[i] = dgamma(U[i], 
                      shape=0.5*(v + 1), 
                      rate=0.5*(v + (1+delta[i]*delta[i])*c[i]^2), 
                      log=TRUE)
  }  
  
  logQ_star = sqrt(U_star)*delta*W*c
  logpU_star = logQ_star + 0.5*(v-1)*log(U_star) - 0.5*U_star*(v+(1-delta*delta)*c)
  
  #MH step
  for(i in 1:n){
    
    #r = pU_star[i]*qU[i]/(pU[i]*qU_star[i])
    #pacc = min(1, r)
    
    logr = logpU_star[i] + logqU[i] - logpU[i] - logqU_star[i]
    logpacc = min(0, logr)
    
    if(log(runif(1)) < logpacc){
      U_end[i] = U_star[i]
      tx[i] = tx[i] + 1
    }
  }
  
  return(list(U=U_end, tx=tx))
  
}
z = postU(U=U, y=y, h=h, W=W, a=a, v=v, tx=rep(0, T))
plot(z$U)
sum(z$tx)/T

# h
posth=function(y, h, U, W, a, v, theta2, counth){
  
  time = Sys.time()
  hout = h
  n=length(y)
  delta = a/sqrt(1+a*a)
  k1 = sqrt(0.5*v)*gamma(0.5*(v-1))/gamma(0.5*v)
  k2 = v/(v-2)
  omega = 1/sqrt(k2-(2/pi)*(delta*k1)^2)
  gammat = -sqrt(2/pi)*delta*omega*k1
  c1 = gammat + omega*delta*W/sqrt(U)
  c2 = omega*(1-delta*delta)/sqrt(U)
  invc22 = 1/c2^2
  
  #theta2=c(mu, phi, s_h)
  mu = theta2[1]
  phi = theta2[2]
  sig2 = theta2[3]^2
  
  delta_tilde = c(mu, rep(mu * (1-phi), T-1))
  Hphi = diag(n)
  for(t in 2:n) Hphi[t, t-1] = -phi
  d_values = c((1 - phi^2) / sig2, rep(1/sig2, T-1))
  D = diag(d_values)
  HinvSH = t(Hphi) %*% D %*% Hphi
  deltah = solve(Hphi, delta_tilde)
  HinvSHdeltah = HinvSH %*% deltah
  
  # Optimization step
  ht = h
  errh = 1
  NRstep = 0.1
  while(errh > 1e-3){
    
    yexpht = y*exp(-0.5*ht)
    fh = -0.5 + 0.5*yexpht*invc22*(yexpht-c1)
    # Newton-Rapshon
    #Gh = -0.5*y^2*exp(-ht)/c2^2 + 0.25*y*exp(-0.5*ht)/c2^2
    # Fisher Score
    Gh = 0.25*(2*(c1^2 + c2^2) - c1)*invc22
    
    Gdiag = diag(Gh)
    kh_tilde = HinvSH + Gdiag
    kh = fh + Gh * ht + HinvSHdeltah 
    hp = solve(kh_tilde, kh)
    newht = (1-NRstep)*ht + NRstep*hp
    errh = max(abs(newht-ht))
    ht = newht

  }
  cat('optim step: ok!')
  cholHh = chol(kh_tilde)
  
  # AR step
  hstar = ht
  uh = hstar-deltah
  logc = -0.5*crossprod(uh, HinvSH%*%uh) -0.5*sum(hstar/c2)
  logc = logc - 0.5*crossprod(exp(-hstar)*invc22,(y-c1*exp(0.5*hstar))^2)
  logc = logc + log(3)
  flag = 0
  while(flag==0){
    
    hc = ht + cholHh %*% rnorm(n)
    vhc = hc-ht
    uhc = hc-deltah
    alpARc = -0.5*crossprod(uhc, HinvSH%*%uhc) - 0.5*sum(hc/c2)
    alpARc = alpARc - 0.5*crossprod(exp(-hc)*invc22, (y-c1*exp(0.5*hc))^2)
    alpARc = alpARc - logc
    alpARc = alpARc + 0.5*crossprod(vhc, kh_tilde%*%vhc)
    
    if(alpARc>log(runif(1))) flag = 1
    
  }
  
  # MH step
  vh = h-ht
  uh = h-deltah
  alpAR = -0.5*crossprod(uh, HinvSH%*%uh) -0.5*sum(h/c2)
  alpAR = alpAR - 0.5*crossprod(exp(-h)*invc22, (y-c1*exp(0.5*h))^2)
  alpAR = alpAR - logc
  alpAR = as.numeric(alpAR)
  
  if(alpAR < 0){
    alpMH = 1
  }else{
    if(alpARc < 0){
      alpMH = - alpAR;
    }else alpMH = alpARc - alpAR;
  } 
  
  if(alpMH > log(runif(1))){
    hout = hc
    counth = counth + 1
  }
  time = Sys.time() - time
  
  return(list(hout=hout, time=time, counth=counth))
}
z = posth(y=y, h=h, U=U, W=W, a=a, v=v, theta2=c(mu,phi_h,s_h), counth=0)

#z$time
ht = z$hout
z$counth
z$time

plot(h, type='l', col='red',ylim=c(-25,25))
lines(ht)
