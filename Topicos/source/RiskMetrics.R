RiskMetrics=function(ht, at, theta, yobs, alpha=0.05, model){
  
  # model=c(sn, st, ss) 
  rtnorm = function(n){
    u = runif(n)
    return(qnorm(0.5*(u+1)))
  }
  
  mu=theta[1,]
  phi=theta[2,]
  sh=theta[3,]
  sa=theta[4,]  
  if(model != 'sn') v=theta[5,]
  
  N=length(ht)
  newy=newa=newW=newU=newh=delta=k1=k2=omega=gammat=mut=st=numeric(N)
  var=es=matrix(0,length(alpha),N)
  #es=numeric(N)
  lpdsstar=0
  
  for(i in 1:N){
    newh[i] = mu[i] + phi[i]*(ht[i]-mu[i]) + sh[i]*rnorm(1)
    newW[i] = rtnorm(1)
    newa[i] = at[i] + sa[i]*rnorm(1)
    delta[i] = newa[i]/sqrt(1+newa[i]*newa[i])
    ############################################################################
    if(model=='st'){
      newU[i] = rgamma(1, shape=0.5*v[i], rate=0.5*v[i])
      k1[i] = sqrt(0.5*v[i])*gamma(0.5*(v[i]-1))/gamma(0.5*v[i])
      k2[i] = v[i]/(v[i]-2)
    }else{
      if(model=='ss'){
        newU[i] = rbeta(1, shape1=v[i], shape2=1)
        k1[i] = 2*v[i]/(2*v[i]-1)
        k2[i] = v[i]/(v[i]-1)
      }else{
        newU[i] = 1.0
        k1[i] = 1.0
        k2[i] = 1.0
      }
    }
    ############################################################################
    omega[i] = 1/sqrt(k2[i]-(2/pi)*(delta[i]*k1[i])^2)
    gammat[i] = -sqrt(2/pi)*delta[i]*omega[i]*k1[i]
    mut[i] = gammat[i] + omega[i]*delta[i]*newW[i]*exp(0.5*newh[i])/sqrt(newU[i])
    st[i] = omega[i]*sqrt(1-delta[i]^2)*exp(0.5*newh[i])/sqrt(newU[i])
    # y_T+1
    newy=mut[i]+st[i]*rnorm(1e3)
    # VaR
    var[,i]=quantile(newy, probs=alpha)
    # ES
    for(j in 1:length(alpha)){
      es[j,i]=mean(newy[which(newy<var[j,i])])  
    }
    # LPDS*
    lpdsstar=lpdsstar+dnorm(yobs, mut[i], st[i])
  }
  
  # LPDS*
  lpdsstar=lpdsstar/N 
  lpdsstar=log(lpdsstar)
  
  return(list(var=apply(var,1,mean), es=apply(es,1,mean), lpdsstar=lpdsstar))
}
