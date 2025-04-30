pds=function(ht, at, theta, yobs){
  
  rtnorm = function(n){
    u = runif(n)
    return(qnorm(0.5*(u+1)))
  }
  
  mu=theta[1,]
  phi=theta[2,]
  sh=theta[3,]
  sa=theta[4,]
  v=theta[5,]
  
  k=length(yobs)
  pds_star=numeric(k)
  N=length(ht)
  newy=matrix(0, N, k)
  
  for(j in 1:k){
    x=0
    newa=newW=newU=newh=delta=k1=k2=omega=gammat=mut=st=numeric(N)
    for(i in 1:N){
      newh[i] = mu[i] + phi[i]*(ht[i]-mu[i]) + sh[i]*rnorm(1)
      newU[i] = rgamma(1, shape=0.5*v[i], rate=0.5*v[i])
      if(is.na(newU[i])) cat( v[i])
      newW[i] = rtnorm(1)
      newa[i] = at[i] + sa[i]*rnorm(1)
      
      delta[i] = newa[i]/sqrt(1+newa[i]*newa[i])
      k1[i] = sqrt(0.5*v[i])*gamma(0.5*(v[i]-1))/gamma(0.5*v[i])
      k2[i] = v[i]/(v[i]-2)
      omega[i] = 1/sqrt(k2[i]-(2/pi)*(delta[i]*k1[i])^2)
      gammat[i] = -sqrt(2/pi)*delta[i]*omega[i]*k1[i]
      mut[i] = gammat[i] + omega[i]*delta[i]*newW[i]*exp(0.5*newh[i])/sqrt(newU[i])
      st[i] = omega[i]*sqrt(1-delta[i]^2)*exp(0.5*newh[i])/sqrt(newU[i])
      newy[i,j] = mut[i] + st[i]*rnorm(1)
      x = x + dnorm(newy[i,j], mut[i], st[i])
      pds_star[j] = x
    }
  }
  return(list(newy=newy, pds_star=pds_star/N) )
}
