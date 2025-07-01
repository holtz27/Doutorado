RiskMetrics=function(ht, at, theta){
  
  # model=c(sn, st, ss) 
  rtnorm=function(n){
    u=runif(n)
    return(qnorm(0.5*(u+1)))
  }
  
  N=length(ht)
  mu=theta[1,]
  phi=theta[2,]
  sh=theta[3,]
  sa=theta[4,]
  v=theta[5,]
  
  newy=numeric(N)
  newa=newW=newU=newh=delta=k1=k2=omega=gammat=mut=st=numeric(N)
  
  for(i in 1:N){
    newh[i] = mu[i] + phi[i]*(ht[i]-mu[i]) + sh[i]*rnorm(1)
    newW[i] = rtnorm(1)
    newa[i] = at[i]+sa[i]*rnorm(1)
    delta[i] = newa[i]/sqrt(1+newa[i]*newa[i])
    ############################################################################
    newU[i] = rgamma(1, shape=0.5*v[i], rate=0.5*v[i])
    k1[i] = sqrt(0.5*v[i])*gamma(0.5*(v[i]-1))/gamma(0.5*v[i])
    k2[i] = v[i]/(v[i]-2)
    ############################################################################
    omega[i] = 1/sqrt(k2[i]-(2/pi)*(delta[i]*k1[i])^2)
    gammat[i] = -sqrt(2/pi)*delta[i]*omega[i]*k1[i]
    mut[i] = gammat[i] + omega[i]*delta[i]*newW[i]*exp(0.5*newh[i])/sqrt(newU[i])
    st[i] = omega[i]*sqrt(1-delta[i]^2)*exp(0.5*newh[i])/sqrt(newU[i])
    # y_T+1
    newy[i]=mut[i]+st[i]*rnorm(1)
  }
  
  return(list(newa=newa,
              newh=newh,
              newy=newy))
}
