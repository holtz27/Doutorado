################################################################################
#rm(list=ls(all=TRUE))
library(parallel)
library(mvtnorm)
Rcpp::sourceCpp("~/Documentos/svmHMM/mlogLk_Rcpp.cpp")
Rcpp::sourceCpp("~/Documentos/svmHMM/pdf_vg.cpp")
### Viterbi
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/HMMcodes/source/svm.viterbi.R')
### Diagnostic
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/HMMcodes/source/ISdiag.R')
################################################################################
num_cores=detectCores(logical=FALSE) 
RcppParallel::setThreadOptions(numThreads=num_cores-1) 
################################################################################
svmvg.sim=function(mu,phi,sigma,nu,beta,y0,g_dim){
  
  y=array(0,dim=g_dim)  
  h=array(0,dim=g_dim)    
  l=invgamma::rinvgamma(g_dim,0.5*nu,0.5*nu)
  
  b0=beta[1]
  b1=beta[2]
  b2=beta[3]
  
  h[1]=rnorm(1,mu,sigma/sqrt(1-phi^2))
  y[1]= b0+b1*y0+b2*exp(h[1])+(exp(0.5*h[1])/sqrt(l[1]))*rnorm(1)
  for (j in 2:g_dim){
    h[j]=mu+phi*(h[j-1]-mu)+rnorm(1,0,sigma)
    y[j]=b0+b1*y[j-1]+b2*exp(h[j])+(exp(0.5*h[j])/sqrt(l[j]))*rnorm(1)
  }
  
  return (list(y=y,h=h,l=l))
}
svm.pn2pw <- function(beta, mu, phi, sigma, nu){
  lbeta1<- beta[1]
  lbeta2<-log((1+beta[2])/(1-beta[2]))
  lbeta3<-beta[3]
  lmu<-mu
  lphi <- log((1+phi)/(1-phi))
  lsigma <- log(sigma)
  # 2<nu<40
  #lnu = log(nu-2)-log(40-nu)
  alpha=0.1
  lnu=(2/alpha)*atanh( (2*nu-40-2)/(40-2))
  parvect = c(lbeta1,lbeta2,lbeta3,lmu,lphi,lsigma,lnu)
  return(parvect)
}
svm.pw2pn <- function(parvect){
  beta=array(0,dim=3)
  beta[1]= parvect[1]
  beta[2]=(exp(parvect[2])-1)/(exp(parvect[2])+1)
  beta[3]=parvect[3]
  mu=parvect[4]
  phi=(exp(parvect[5])-1)/(exp(parvect[5])+1)
  sigma=exp(parvect[6])
  #nu = exp(parvect[7]) + 2
  #nu = (40*exp(parvect[7])+2)/(1+exp(parvect[7]))
  alpha=0.1
  nu=0.5*((40-2)*tanh(0.5*alpha*parvect[7])+(40+2))
  #return(list(beta=beta,mu=mu,phi=phi,sigma=sigma,nu=nu))
  return(c(beta, mu, phi, sigma, nu))
}
## function that will be used to compute 'allprobs' in mllk below
fillallprobs=function(x,beg,beta,nu,y){
  z = (x-beta[1]-beta[2]*y-beta[3]*beg^2)/beg 
  #w = pdf_vg(y=z, mu=0, sigma=1, nu=nu)  
  w = pdf_vg(y=z, nu=nu)
  return(w/beg)  
}
svmvg.mllk=function(parvect, y, y0, m, gmax){
  ny = length(y)
  p = svm.pw2pn(parvect)
  K = m+1
  b=seq(-gmax,gmax,length=K) 
  bs=(b[-1]+b[-K])*0.5
  E=p[4]+p[5]*(bs-p[4])
  intlen <- b[2]-b[1]
  sey = exp(bs/2) 
  Gamma=matrix(0,m,m) #06
  for (i in 1:m){
    Gamma[i,]=dnorm(bs, mean=E[i], sd=p[6])
    sg = sum(Gamma[i,])
    if(sg == 0){
      stop('Built Gamma error')
    }else{
      Gamma[i,] = Gamma[i,]/sg
    } 
  }
  Gamma = intlen*Gamma
  #Gamma = Gamma/apply(Gamma,1,sum) 
  
  xx<-y
  yy<-c(y0,y[1:(ny-1)])
  allprobs = outer(xx,sey,"fillallprobs", beta=p[1:3], nu=p[7], yy)
  delta=dnorm(bs,p[4],p[6]/sqrt(1-p[5]^2))*intlen
  foo = delta*allprobs[1,]
  lscale = mlogLk_Rcpp(allprobs,Gamma,foo,ny) #Rcpp function
  return(-lscale)
}
svmvg.prior=function(parvect){
  lprior = log(dnorm(parvect[1], 0, 10))
  + log(dnorm(parvect[2], 0.5, 10))
  + log(dnorm(parvect[3], 0, 10))
  + log(dnorm(parvect[4], 0, 10))
  + log(dnorm(parvect[5], 4.5, 10))
  + log(dnorm(parvect[6], -1.5, 10)) 
  + log(dnorm(parvect[7], -10, 10))
  return(-lprior)  
}
svmvg.posterior=function(parvect, y, y0, m, gmax){
  return(svmvg.mllk(parvect,y,y0,m,gmax)+svmvg.prior(parvect))  
}
svmvg.map=function(y, m, parvect0, y0, gmax){
  
  #parvect <- svm.pn2pw(beta=beta0,mu=mu0,phi=phi0,sigma=sigma0,nu=nu0)
  mod=optim(parvect0, svmvg.posterior, y=y, y0=y0, m=m, gmax=gmax,hessian=TRUE)
  mode=mod$par
  conv=mod$convergence
  return(list(mode=mode,
              lpostsvvg=-mod$value,
              hessian=mod$hessian,
              conv=conv))
}
################################################################################
mu0=0
phi0=0.98
sigma0 = 0.15
beta0=c(0.1 , 0.1 ,-0.05)
nu0=10
y0=0.2
################################################################################
m=200
gmax=5.0

try({
  time = Sys.time()
  ############################################################################
  # Optimize
  eigenvalues=rep(0, 7)
  cont=0
  parvect0=svm.pn2pw(beta0, mu0, phi0, sigma0, nu0)
  while(any(eigenvalues <= 0) && cont<5){
    
    optim_res=svmvg.map(y=ytrain, m=m, parvect0=parvect0, y0=y0, gmax=gmax)
    H1=signif(solve(optim_res$hessian), 6)
    eigenvalues=eigen(H1, only.values=TRUE)$values
    parvect0=optim_res$mode
    cont=cont+1
    
  }
  ############################################################################
  # Test if H1 is positive definite
  if(any(eigenvalues <= 0)){
    stop('Hessian is not positive definite')
  }else{
    k=-optim_res$lpostsvvg
    map=optim_res$mode
    ##########################################################################
    # Weigths Evaluation
    n=1e3
    X=rmvnorm(n, map, H1)
    Weigth <- array(0,dim=c(n,1))
    for(j in 1:n){
      Weigth[j,1]=exp(k
                      -svmvg.posterior(parvect=X[j,],y=ytrain,y0=y0,m=m,gmax=gmax)
                      -dmvnorm(X[j,],mean=map,sigma=H1,log=TRUE))
    }
    s = sum(Weigth)
    if((s != 0) && !is.nan(s)){
      Weigth=Weigth/s  
    }else{
      stop('Error normalize constante weigths!')
    }
    ### Resample
    indx=sample(1:n, prob=Weigth, replace=TRUE)
    X=X[indx,]
    Weigth=rep(1/n, n)
    Results=ISdiag(Weigth=Weigth, X=X)
    times=as.numeric(Sys.time()-time, units='mins')
  }
})
Results
times

h_hat=svm.viterbi(y=ytrain,theta_hat=Results$Results[,1],m=m,gmax=gmax)
plot(abs(ytrain), type='l', col='gray')
lines(exp(0.5*h_hat))

### -log p(y|\theta)
# DIC
theta_hat=apply(X, 2, mean)
D=2*svmvg.mllk(parvect=theta_hat,y=ytrain,y0=y0,m=m,gmax=gmax)
D

# Evaluating \bar{D(\theta)}
Dbar=0
for(j in 1:n){
  pv=X[j,]
  Dbar=Dbar+Weigth[j]*svmvg.mllk(parvect=pv,y=ytrain,y0=y0,m=m,gmax=gmax)
}
Dbar=2*Dbar
pd=Dbar-D
DIC=D+2*pd
DIC

# Log Predictive Score
LPS=0.5*D/length(ytrain)
LPS
################################################################################
cdf_vg=function(x, nu){
  
  n=length(x)
  z=numeric(n)
  for(i in 1:n){
    z[i]=integrate(pdf_vg, lower=-Inf, upper=x[i], nu)$value  
  }
  return(z)
}
SVMvg.HMM.lalpha=function(x, m, gbmax, pvect, y0){
  
  nx <- length(x)
  yy = c(y0, x[1:(nx-1)])
  p<-svm.pw2pn(pvect)
  
  lalpha <- matrix(NA,m,nx)
  K <- m+1
  gb <- seq(-gbmax,gbmax,length=K)
  g <- (gb[-1]+gb[-K])*0.5             
  beg <-exp(g/2)   
  gamma <- matrix(0,m,m)
  E=p[4]+p[5]*(g-p[4])
  intlen <- gb[2]-gb[1]
  for (i in 1:m){
    goo <- dnorm(g,E[i],p[6])*intlen
    gamma[i,] <- goo/sum(goo)
  }
  lscale=0
  delta <- dnorm(g,p[4],p[6]/sqrt(1-p[5]^2))*intlen
  for (i in 1:nx){
    if (i==1) foo <- delta*1/beg*pdf_vg((x[1]-p[1]-p[2]*yy[1]-p[3]*beg^2)/beg, p[7])
    if ( i>1) 
      foo <- foo%*%gamma*1/beg*pdf_vg((x[i]-p[1]-p[2]*yy[i]-p[3]*beg^2)/beg, p[7])
    lscale <- lscale+log(sum(foo));
    foo <- foo/sum(foo)
    lalpha[,i] <- log(foo)+lscale 
  }
  return(lalpha)
}
SVMvg.HMM.quantiles=function(x, res.time, m, gbmax, pvect, lalpha, y0){
  nx = length(x)
  yy = c(y0, x[1:(nx-1)])
  K = m+1
  p=svm.pw2pn(pvect)
  gb = seq(-gbmax,gbmax,length=K)
  g = (gb[-1]+gb[-K])*0.5             
  beg = exp(g/2)   
  gamma = matrix(0,m,m)
  E = p[4]+p[5]*(g-p[4])
  intlen = gb[2]-gb[1]
  for(i in 1:m){
    goo = dnorm(g,E[i],p[6])*intlen
    gamma[i,] = goo/sum(goo)
  }
  c=max(lalpha[,res.time-1])   # to avoid numerical underflow
  a=exp(lalpha[,res.time-1]-c) # scalar factor cancels in Pro computation below
  npsr=t(a)%*%(gamma/sum(a))%*%cdf_vg((x[res.time]-p[1]-p[2]*yy[res.time]-p[3]*beg^2)/beg, p[7])
  
  return(npsr)
}
lalphasvmvg=SVMvg.HMM.lalpha(x=ytrain, m=m, gbmax=gmax, pvect=theta_hat, y0=y0)

psressvmvg=rep(NA,length(ytest))
for(k in 1:length(ytest)){ 
  psressvmvg[k]=SVMvg.HMM.quantiles(ytrain, res.time=length(ytest)+k, m=m, gbmax=gmax, 
                                    pvect=theta_hat, lalphasvmvg, y0=y0)
}
qqnorm(qnorm(psressvmvg), main='SVM-VG') 
qqline(qnorm(psressvmvg))
require(tseries)
pvalue=jarque.bera.test(qnorm(psressvmvg))$p.value
if(pvalue>0.05){
  cat('Not reject the hypothesis of normality of the residuals at the 5% nivel.')
}else{
  cat('Reject the hypothesis of normality of the residuals at the 5% nivel.')
}

# now backtesting
# as an example, check how many of the observations in the validation sample are
# smaller than the 1%/5%/10%/50% quantiles of the one-step-ahead forecast 
# distributions, under the fitted model
# should be approximately 0.01*length(obs.val)
VaR=c(sum(psressvmvg<0.01), 0.01*length(ytest))
VaR=rbind(VaR, c(sum(psressvmvg<0.05), 0.05*length(ytest)))
VaR=rbind(VaR, c(sum(psressvmvg<0.1), 0.1*length(ytest)))
VaR=rbind(VaR, c(sum(psressvmvg<0.5), 0.5*length(ytest)))
row.names(VaR)=c('0.01', '0.05', '0.10', '0.50')
colnames(VaR)=c('Observed', 'Expected')
VaR
################################################################################
#save(Results, times, h_hat, DIC, LPS, psressvmvg, 
#     file='~/Documentos/svmHMM/sp500/svmvg.RData')

