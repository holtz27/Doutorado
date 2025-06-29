RiskMetrics=function(ht, at, theta, yobs, alpha=0.05){
  
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
  
  k=1e3
  newy=numeric(k)
  newa=newW=newU=newh=delta=k1=k2=omega=gammat=mut=st=numeric(N)
  var=es=matrix(0,length(alpha), 1)
  qs=matrix(0,length(alpha), 1)
  lpdsstar=err=rmse=0
  
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
    newy[i]=mean(mut[i]+st[i]*rnorm(k)) 
  }
  
  #MSE
  err=mean(newy-yobs)
  rmse=mean((newy-yobs)^2)   
  # LPDS*
  lpdsstar=mean(dnorm(yobs, mut, st))
  lpdsstar=log(lpdsstar)
  # VaR, ES
  var[,1]=quantile(newy, probs=alpha)
  # ES
  for(j in 1:length(alpha)){
    es[j,1]=mean(newy[which(newy<var[j,1])])  
  }
  # CRPS
  crps=scoringRules::crps_sample(y=yobs, dat=newy)
  # logS
  logs=scoringRules::logs_sample(y=yobs, dat=newy)
  # QS
  qs[1,1]=scoringRules::qs_quantiles(y=yobs, 
                                     x=quantile(newy, probs=alpha[1]), 
                                     alpha=alpha[1])
  qs[2,1]=scoringRules::qs_quantiles(y=yobs, 
                                     x=quantile(newy, probs=alpha[2]), 
                                     alpha=alpha[2])
  
  return(list(var=var,
              es=es, 
              lpdsstar=lpdsstar,
              crps=crps,
              logs=logs,
              qs=qs,
              err=err, 
              rmse=rmse,
              newy=newy))
}
