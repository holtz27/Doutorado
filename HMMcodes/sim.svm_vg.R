################################################################################
rm(list=ls(all=TRUE))
library("Rcpp")
library("RcppArmadillo")
#library(parallel)
#library(RcppParallel)
#library(RcppNumerical)
library('mvtnorm')
library(invgamma)
sourceCpp("mLogLk_Rcpp.cpp")
sourceCpp("pdf_vg.cpp")
#sourceCpp("pdf2_vg.cpp")
#sourceCpp("parallel_pdf2_vg.cpp")
################################################################################
svm.pn2pw <- function(beta,mu,phi,sigma,nu){
  lbeta1<- beta[1]
  lbeta2<-log((1+beta[2])/(1-beta[2]))
  lbeta3<-beta[3]
  lmu<-mu
  lphi <- log((1+phi)/(1-phi))
  lsigma <- log(sigma)
  lnu = log(nu)
  parvect <- c(lbeta1,lbeta2,lbeta3,lmu,lphi,lsigma,lnu)
  return(parvect)
}
svm.pw2pn <- function(parvect){
  beta=array(0,dim=3)
  beta[1]= parvect[1]
  beta[2]=(exp(parvect[2])-1)/(exp(parvect[2])+1)
  beta[3]=parvect[3]
  mu=parvect[4]
  phi <- (exp(parvect[5])-1)/(exp(parvect[5])+1)
  sigma <- exp(parvect[6])
  nu = exp(parvect[7])
  return(list(beta=beta,mu=mu,phi=phi,sigma=sigma,nu=nu))
}
quantile <- function(x, weights, probs){
  
  ix <- order(x)
  x_sorted <- x[ix]
  weights_sorted <- weights[ix]
  weights_sorted <- weights_sorted / sum(weights_sorted)
  cumulative_weights <- cumsum(weights_sorted)
  quantiles <- numeric(length(probs))
  for (i in seq_along(probs)) {
    index <- which(cumulative_weights >= probs[i])[1]
    quantiles[i] <- x_sorted[index]
  }
  
  return(quantiles)
}
fillallprobs.vg <- function(x,beg,beta,nu,y){
  return((1/beg)*pdf_vg((x-beta[1]-beta[2]*y-beta[3]*beg^2)/beg,0.0,1.0,nu))
}
svmvg.mllk <-function(parvect,y,y0,m,gmax){
  ny <- length(y)
  p <- svm.pw2pn(parvect)
  K= m+1
  b=seq(-gmax,gmax,length=K) 
  bs=(b[-1]+b[-K])*0.5
  E=p$mu+p$phi*(bs-p$mu)
  intlen <- b[2]-b[1]
  sey= exp(bs/2) 
  Gamma=matrix(0,m,m) #06
  for (i in 1:m){
    Gamma[i,]=dnorm(bs,E[i],p$sigma)*intlen
  }
  Gamma= Gamma/apply(Gamma,1,sum) 
  xx<-y
  yy<-c(y0,y[1:(ny-1)])
  allprobs <- outer(xx,sey,"fillallprobs.vg",beta=p$beta,nu=p$nu,yy)
  delta <-dnorm(bs,p$mu,p$sigma/sqrt(1-p$phi^2))*intlen
  foo <- delta*allprobs[1,]
  lscale = mlogLk_Rcpp(allprobs,Gamma,foo,ny) #Rcpp function
  return(-lscale)
}
svmvg.prior <-function(parvect){
  lprior=log(dnorm(parvect[1],0,10))+log(dnorm(parvect[2],0.5,10))
  +log(dnorm(parvect[3],0,10))+log(dnorm(parvect[4],0,10))+log(dnorm(parvect[5],4.5,10))
  +log(dnorm(parvect[6],-1.5,10)) + log(dnorm(parvect[7],0,10))
  return(-lprior)  
}
svmvg.posterior <-function(parvect,y,y0,m,gmax){
  return(svmvg.mllk(parvect,y,y0,m,gmax)+svmvg.prior(parvect))  
}
svmvg.map <- function(y,y0,m,beta0,mu0,phi0,sigma0,nu0,gmax){
  s <- Sys.time()
  parvect <- svm.pn2pw(beta=beta0,mu=mu0,phi=phi0,sigma=sigma0,nu=nu0)
  mod <- optim(parvect,svmvg.posterior,y=y,y0=y0,m=m,gmax=gmax,hessian=T)
  mode <- mod$par
  s = Sys.time()-s
  return(list(mode=mode,lpostsvmvg=-mod$value,hessian=mod$hessian,time=s))
}
svmvg.sim <-function(mu,phi,sigma,nu,beta,y0,g_dim){
  
  y=array(0,dim=g_dim)  
  h=array(0,dim=g_dim)    
  l=rinvgamma(g_dim,0.5*nu,0.5*nu)
  
  b0=beta[1]
  b1=beta[2]
  b2=beta[3]
  
  h[1]=rnorm(1,mu,sigma/sqrt(1-phi^2))
  y[1]= b0+b1*y0+b2*exp(h[1])+(exp(0.5*h[1])/sqrt(l[1]))*rnorm(1)
  for (j in 2:g_dim){
    h[j]=mu+phi*(h[j-1]-mu)+rnorm(1,0,sigma)
    y[j]=b0+b1*y[j-1]+b2*exp(h[j])+(exp(0.5*h[j])/sqrt(l[j]))*rnorm(1)
  }
  
  return (list(y=y,h=h,l=l))}
################################################################################
betasim = c(0.2,0.07,-0.18)
mu = 0.1
phi = 0.98
sigma = 0.1
nu = 10
g_dim = 6000
y0 = 0.2
################################################################################
set.seed(7531956)
################################################################################
mu0=0.5
phi0=0.96
sigma0=0.2
nu0=8
beta0=c(0.5,0.1,-0.2)
gmax=4
################################################################################
simsvmvg1<-svmvg.sim(y0=y0,beta=betasim,mu=mu,phi=phi,sigma=sigma,nu=nu,g_dim=g_dim)
################################################################################
Results = optim_time = weigth_time = list()
i=1
s1=1500
for(m in c(50,100,150,200)){
  res.svmvg200_4_s1 <- svmvg.map(simsvmvg1$y[1:s1],m=m,
                                 beta0=beta0,mu0=mu0,phi0=phi0,
                                 sigma0=sigma0,nu0=nu0,y0=simsvmvg1$y0,
                                 gmax=gmax)
  k = -res.svmvg200_4_s1$lpostsvmvg
  ################################################################################
  H1=signif(solve(res.svmvg200_4_s1$hessian),6)
  n = 500
  X=rmvnorm(n,res.svmvg200_4_s1$mode,H1)
  Weigth <- array(0,dim=c(n,1))
  for(j in 1:n){
    if(j==1) time=Sys.time()
    Weigth[j,1]=exp(k #1000
                    -svmvg.posterior(X[j,],simsvmvg1$y[1:s1],
                                     simsvmvg1$y0,m=m,gmax=gmax)
                    -dmvnorm(X[j,],res.svmvg200_4_s1$mode,sigma=H1,log=T)
    )
    if(j==n) time = Sys.time()-time
  }
  ################################################################################
  Weigth=Weigth/sum(Weigth)
  beta0m=sum(X[,1]*Weigth)
  beta1m=sum(((exp(X[,2])-1)/(exp(X[,2])+1))*Weigth)
  beta2m=sum(X[,3]*Weigth)
  mum=sum(X[,4]*Weigth)
  phim=sum(((exp(X[,5])-1)/(exp(X[,5])+1))*Weigth)
  sigmam=sum(exp(X[,6])*Weigth)
  num=sum((exp(X[,7]))*Weigth)
  ################################################################################
  int_b0 = quantile(X[,1], Weigth, c(0.025, 0.975))
  int_b1 = quantile((exp(X[,2])-1)/(exp(X[,2])+1), Weigth, c(0.025, 0.975))
  int_b2 = quantile(X[,3], Weigth, c(0.025, 0.975))
  int_mu = quantile(X[,4], Weigth, c(0.025, 0.975))
  int_phi = quantile( (exp(X[,5])-1)/(exp(X[,5])+1), Weigth, c(0.025, 0.975))
  int_sigma = quantile(exp(X[,6]), Weigth, c(0.025, 0.975))
  int_nu = quantile( (exp(X[,7])), Weigth, c(0.025, 0.975))
  ################################################################################
  # Results
  results = matrix(c(betasim[1], beta0m, int_b0),nrow = 1)
  results = rbind(results,
                  c(betasim[2], beta1m, int_b1),
                  c(betasim[3], beta2m, int_b2),
                  c(mu, mum, int_mu),
                  c(phi, phim, int_phi),
                  c(sigma, sigmam, int_sigma),
                  c(nu, num, int_nu))
  row.names(results) = c('b0','b1','b2','mu','phi','sigma','nu')
  colnames(results) = c('true','mean','2.5%','97.5%')
  round(results, 4)
  
  Results[[i]] = round(results, 4)
  optim_time[[i]] = res.svmvg200_4_s1$time
  weigth_time[[i]] = time
  i = i + 1
}
save(Results, optim_time, weigth_time, file = 'vgouts.1500.RData')
################################################################################
Results = optim_time = weigth_time = list()
i=1
s1=3000
for(m in c(50,100,150,200)){
  res.svmvg200_4_s1 <- svmvg.map(simsvmvg1$y[1:s1],m=m,
                                 beta0=beta0,mu0=mu0,phi0=phi0,
                                 sigma0=sigma0,nu0=nu0,y0=simsvmvg1$y0,
                                 gmax=gmax)
  k = -res.svmvg200_4_s1$lpostsvmvg
  ################################################################################
  H1=signif(solve(res.svmvg200_4_s1$hessian),6)
  n = 500
  X=rmvnorm(n,res.svmvg200_4_s1$mode,H1)
  Weigth <- array(0,dim=c(n,1))
  for(j in 1:n){
    if(j==1) time=Sys.time()
    Weigth[j,1]=exp(k #1000
                    -svmvg.posterior(X[j,],simsvmvg1$y[1:s1],
                                     simsvmvg1$y0,m=m,gmax=gmax)
                    -dmvnorm(X[j,],res.svmvg200_4_s1$mode,sigma=H1,log=T)
    )
    if(j==n) time = Sys.time()-time
  }
  ################################################################################
  Weigth=Weigth/sum(Weigth)
  beta0m=sum(X[,1]*Weigth)
  beta1m=sum(((exp(X[,2])-1)/(exp(X[,2])+1))*Weigth)
  beta2m=sum(X[,3]*Weigth)
  mum=sum(X[,4]*Weigth)
  phim=sum(((exp(X[,5])-1)/(exp(X[,5])+1))*Weigth)
  sigmam=sum(exp(X[,6])*Weigth)
  num=sum((exp(X[,7]))*Weigth)
  ################################################################################
  int_b0 = quantile(X[,1], Weigth, c(0.025, 0.975))
  int_b1 = quantile((exp(X[,2])-1)/(exp(X[,2])+1), Weigth, c(0.025, 0.975))
  int_b2 = quantile(X[,3], Weigth, c(0.025, 0.975))
  int_mu = quantile(X[,4], Weigth, c(0.025, 0.975))
  int_phi = quantile( (exp(X[,5])-1)/(exp(X[,5])+1), Weigth, c(0.025, 0.975))
  int_sigma = quantile(exp(X[,6]), Weigth, c(0.025, 0.975))
  int_nu = quantile( (exp(X[,7])), Weigth, c(0.025, 0.975))
  ################################################################################
  # Results
  results = matrix(c(betasim[1], beta0m, int_b0),nrow = 1)
  results = rbind(results,
                  c(betasim[2], beta1m, int_b1),
                  c(betasim[3], beta2m, int_b2),
                  c(mu, mum, int_mu),
                  c(phi, phim, int_phi),
                  c(sigma, sigmam, int_sigma),
                  c(nu, num, int_nu))
  row.names(results) = c('b0','b1','b2','mu','phi','sigma','nu')
  colnames(results) = c('true','mean','2.5%','97.5%')
  round(results, 4)
  
  Results[[i]] = round(results, 4)
  optim_time[[i]] = res.svmvg200_4_s1$time
  weigth_time[[i]] = time
  i = i + 1
}
save(Results, optim_time, weigth_time, file = 'vgouts.3000.RData')
################################################################################
Results = optim_time = weigth_time = list()
i=1
s1=6000
for(m in c(50,100,150,200)){
  res.svmvg200_4_s1 <- svmvg.map(simsvmvg1$y[1:s1],m=m,
                                 beta0=beta0,mu0=mu0,phi0=phi0,
                                 sigma0=sigma0,nu0=nu0,y0=simsvmvg1$y0,
                                 gmax=gmax)
  k = -res.svmvg200_4_s1$lpostsvmvg
  ################################################################################
  H1=signif(solve(res.svmvg200_4_s1$hessian),6)
  n = 500
  X=rmvnorm(n,res.svmvg200_4_s1$mode,H1)
  Weigth <- array(0,dim=c(n,1))
  for(j in 1:n){
    if(j==1) time=Sys.time()
    Weigth[j,1]=exp(k #1000
                    -svmvg.posterior(X[j,],simsvmvg1$y[1:s1],
                                     simsvmvg1$y0,m=m,gmax=gmax)
                    -dmvnorm(X[j,],res.svmvg200_4_s1$mode,sigma=H1,log=T)
    )
    if(j==n) time = Sys.time()-time
  }
  ################################################################################
  Weigth=Weigth/sum(Weigth)
  beta0m=sum(X[,1]*Weigth)
  beta1m=sum(((exp(X[,2])-1)/(exp(X[,2])+1))*Weigth)
  beta2m=sum(X[,3]*Weigth)
  mum=sum(X[,4]*Weigth)
  phim=sum(((exp(X[,5])-1)/(exp(X[,5])+1))*Weigth)
  sigmam=sum(exp(X[,6])*Weigth)
  num=sum((exp(X[,7]))*Weigth)
  ################################################################################
  int_b0 = quantile(X[,1], Weigth, c(0.025, 0.975))
  int_b1 = quantile((exp(X[,2])-1)/(exp(X[,2])+1), Weigth, c(0.025, 0.975))
  int_b2 = quantile(X[,3], Weigth, c(0.025, 0.975))
  int_mu = quantile(X[,4], Weigth, c(0.025, 0.975))
  int_phi = quantile( (exp(X[,5])-1)/(exp(X[,5])+1), Weigth, c(0.025, 0.975))
  int_sigma = quantile(exp(X[,6]), Weigth, c(0.025, 0.975))
  int_nu = quantile( (exp(X[,7])), Weigth, c(0.025, 0.975))
  ################################################################################
  # Results
  results = matrix(c(betasim[1], beta0m, int_b0),nrow = 1)
  results = rbind(results,
                  c(betasim[2], beta1m, int_b1),
                  c(betasim[3], beta2m, int_b2),
                  c(mu, mum, int_mu),
                  c(phi, phim, int_phi),
                  c(sigma, sigmam, int_sigma),
                  c(nu, num, int_nu))
  row.names(results) = c('b0','b1','b2','mu','phi','sigma','nu')
  colnames(results) = c('true','mean','2.5%','97.5%')
  round(results, 4)
  
  Results[[i]] = round(results, 4)
  optim_time[[i]] = res.svmvg200_4_s1$time
  weigth_time[[i]] = time
  i = i + 1
}
save(Results, optim_time, weigth_time, file = 'vgouts.6000.RData')
################################################################################
