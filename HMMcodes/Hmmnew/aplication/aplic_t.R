################################################################################
#rm(list=ls(all=TRUE))
library(parallel)
library(mvtnorm)
Rcpp::sourceCpp("~/HMMnew/aplication/mlogLk_Rcpp.cpp")
Rcpp::sourceCpp("~/HMMnew/aplication/pdf_t.cpp")
### Viterbi
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/HMMcodes/source/svm.viterbi.R')
### Diagnostic
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/HMMcodes/source/ISdiag.R')
################################################################################
num_cores=detectCores(logical=FALSE) 
RcppParallel::setThreadOptions(numThreads=num_cores-1) 
################################################################################
double_eps <- .Machine$double.eps
double_xmax <- .Machine$double.xmax
double_xmin <- .Machine$double.xmin
log_double_xmax <- log(double_xmax)
log_double_xmin <- log(double_xmin)
################################################################################
svmt.sim <- function(mu,phi,sigma,nu,beta,y0,g_dim){
  
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
svm.pn2pw <- function(beta,mu,phi,sigma,nu){
  lbeta1<- beta[1]
  lbeta2<-log((1+beta[2])/(1-beta[2]))
  lbeta3<-beta[3]
  lmu<-mu
  lphi <- log((1+phi)/(1-phi))
  lsigma <- log(sigma)
  # nu > 2
  #lnu = log(nu-2)
  # 2<nu<40
  lnu = log(nu-2)-log(40-nu)
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
  #nu = exp(parvect[7]) + 2
  nu = (40*exp(parvect[7])+2)/(1+exp(parvect[7]))
  return(list(beta=beta,mu=mu,phi=phi,sigma=sigma,nu=nu))
}
## function that will be used to compute 'allprobs' in mllk below
fillallprobs <- function(x,beg,beta,nu,y){
  z = (x-beta[1]-beta[2]*y-beta[3]*beg^2)/beg 
  #w = pdf_vg(y=z, mu=0, sigma=1, nu=nu)  
  w = pdf_t(y=z, df=nu)
  return(w/beg)  
}
svmt.mllk <-function(parvect,y,y0,m,gmax){
  ny <- length(y)
  p <- svm.pw2pn(parvect)
  K = m+1
  b=seq(-gmax,gmax,length=K) 
  bs=(b[-1]+b[-K])*0.5
  E=p$mu+p$phi*(bs-p$mu)
  intlen <- b[2]-b[1]
  sey = exp(bs/2) 
  Gamma=matrix(0,m,m) #06
  for (i in 1:m){
    Gamma[i,]=dnorm(bs, mean=E[i], sd=p$sigma)
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
  allprobs <- outer(xx,sey,"fillallprobs",beta=p$beta,nu=p$nu,yy)
  delta <-dnorm(bs,p$mu,p$sigma/sqrt(1-p$phi^2))*intlen
  foo <- delta*allprobs[1,]
  lscale = mlogLk_Rcpp(allprobs,Gamma,foo,ny) #Rcpp function
  return(-lscale)
}
svmt.prior <-function(parvect){
  lprior = log(dnorm(parvect[1], 0, 10))
  + log(dnorm(parvect[2], 0.5, 10))
  + log(dnorm(parvect[3], 0, 10))
  + log(dnorm(parvect[4], 0, 10))
  + log(dnorm(parvect[5], 4.5, 10))
  + log(dnorm(parvect[6], -1.5, 10)) 
  + log(dnorm(parvect[7], 0, 1))
  return(-lprior)  
}
svmt.posterior <-function(parvect,y,y0,m,gmax){
  return(svmt.mllk(parvect,y,y0,m,gmax)+svmt.prior(parvect))  
}
svmt.map <- function(y, m, parvect0, y0, gmax){
  
  #parvect <- svm.pn2pw(beta=beta0,mu=mu0,phi=phi0,sigma=sigma0,nu=nu0)
  mod <- optim(parvect0,svmt.posterior,y=y,y0=y0,m=m,gmax=gmax,hessian=TRUE)
  mode <- mod$par
  conv = mod$convergence
  return(list(mode=mode,
              lpostsvvg=-mod$value,
              hessian=mod$hessian,
              conv = conv))
}
################################################################################
mu0=0
phi0=0.98
sigma0 = 0.15
beta0=c(0.1 ,0.1 ,-0.05)
nu0=10
y0=0.2
################################################################################
m=200
gmax=5.0
y=log.ret
h=1e3

ytrain=y[1:(length(y)-h)]
ytest=y[(length(y)-h+1):length(y)]

try({
  time = Sys.time()
  ############################################################################
  # Optimize
  eigenvalues = rep(0, 7)
  cont = 0
  parvect0 = svm.pn2pw(beta0, mu0, phi0, sigma0, nu0)
  while( any(eigenvalues <= 0) && cont<5 ){
    
    optim_res = svmt.map(y=y,m=m, parvect0=parvect0, y0=y0, gmax=gmax)
    H1 = signif(solve(optim_res$hessian), 6)
    eigenvalues = eigen(H1, only.values=TRUE)$values
    parvect0 = optim_res$mode
    cont = cont + 1
    
  }
  ############################################################################
  # Test if H1 is positive definite
  if(any(eigenvalues <= 0)){
    stop('Hessian is not positive definite')
  }else{
    k = -optim_res$lpostsvvg
    map = optim_res$mode
    ##########################################################################
    # Weigths Evaluation
    n = 1e3
    X=rmvnorm(n, map, H1)
    #X=rmvt(n, delta=map, sigma=H1, df=df, type='shifted')
    Weigth=array(0,dim=c(n,1))
    for(j in 1:n){
      Weigth[j,1]=exp(k
                      -svmt.posterior(parvect=X[j,],y=y,y0=y0,m=m,gmax=gmax)
                      -dmvnorm(X[j,],mean=map,sigma=H1,log=TRUE)
                      #-dmvt(X[j,],delta=map,sigma=H1,type='shifted',log=TRUE, df=df)
      )
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
    #
    Results=ISdiag(Weigth=Weigth, X=X, nu.lower=2, nu.upper=40)
    times=as.numeric(Sys.time()-time, units='mins')
  }
})
Results
times

h_hat=svm.viterbi(y=y,theta_hat=Results$Results[,1],m=m,gmax=gmax)
plot(abs(y), type='l', col='gray')
lines(exp(0.5*h_hat))

### -log p(y|\theta)
# DIC
theta_hat=apply(X, 2, mean)
D=2*svmt.mllk(parvect=theta_hat,y=y,y0=y0,m=m,gmax=gmax)
D

# Evaluating \bar{D(\theta)}
Dbar=0
for(j in 1:n){
  pv=X[j,]
  Dbar=Dbar+Weigth[j]*svmt.mllk(parvect=pv,y=y,y0=y0,m=m,gmax=gmax)
}
Dbar=2*Dbar
Dbar

pd=Dbar-D
DIC=D+2*pd
DIC

# Log Predictive Score
LPS=0.5*D/length(y)
LPS
################################################################################
SVMt.HMM.lalpha=function(x, m, gbmax, pvect, y0){
  
  nx <- length(x)
  yy = c(y0, x[1:(nx-1)])
  p<-svm.pw2pn(pvect)
  
  lalpha <- matrix(NA,m,nx)
  K <- m+1
  gb <- seq(-gbmax,gbmax,length=K)
  g <- (gb[-1]+gb[-K])*0.5             
  beg <-exp(g/2)   
  gamma <- matrix(0,m,m)
  E=p$mu+p$phi*(g-p$mu)
  intlen <- gb[2]-gb[1]
  for (i in 1:m){
    goo <- dnorm(g,E[i],p$sigma)*intlen
    gamma[i,] <- goo/sum(goo)
  }
  lscale=0
  delta <- dnorm(g,p$mu,p$sigma/sqrt(1-p$phi^2))*intlen
  for (i in 1:nx){
    if (i==1) foo <- delta*1/beg*dt((x[1]-p$beta[1]-p$beta[2]*yy[1]-p$beta[3]*beg^2)/beg, p$nu)
    if ( i>1) 
      foo <- foo%*%gamma*1/beg*dt((x[i]-p$beta[1]-p$beta[2]*yy[i]-p$beta[3]*beg^2)/beg, p$nu)
    lscale <- lscale+log(sum(foo));
    foo <- foo/sum(foo)
    lalpha[,i] <- log(foo)+lscale 
  }
  return(lalpha)
}
SVMt.HMM.quantiles=function(x, res.time, m, gbmax, pvect, lalpha, y0){
  nx = length(x)
  yy = c(y0, x[1:(nx-1)])
  K = m+1
  p=svm.pw2pn(pvect)
  gb = seq(-gbmax,gbmax,length=K)
  g = (gb[-1]+gb[-K])*0.5             
  beg = exp(g/2)   
  gamma = matrix(0,m,m)
  E = p$mu+p$phi*(g-p$mu)
  intlen = gb[2]-gb[1]
  for(i in 1:m){
    goo = dnorm(g,E[i],p$sigma)*intlen
    gamma[i,] = goo/sum(goo)
  }
  c=max(lalpha[,res.time-1])   # to avoid numerical underflow
  a=exp(lalpha[,res.time-1]-c) # scalar factor cancels in Pro computation below
  npsr=t(a)%*%(gamma/sum(a))%*%pt((x[res.time]-p$beta[1]-p$beta[2]*yy[res.time]-p$beta[3]*beg^2)/beg, p$nu)
  
  return(npsr)
}
lalphasvmt=SVMt.HMM.lalpha(x=ytrain, m=m, gbmax=gmax, pvect=theta_hat, y0=y0)


psressvmt=rep(NA,length(ytest))
for(k in 1:length(ytest)){ 
  # note that the function 'SV.HMM.quantiles' could alternatively be written 
  # such that no loop would be required here, but I didn't bother since it's 
  # fast enough as it stands
  psressvmt[k]=SVMt.HMM.quantiles(y, res.time=length(ytest)+k, m=m, gbmax=gmax, 
                                  pvect=theta_hat, lalphasvmt, y0=y0)
}
# check goodness-of-fit of model:
#qqnorm(qnorm(psressvn)) # 'qnorm(psres)', the one-step-ahead forecast 
#pseudo-residuals, under the model fitted to the calibration data - if the model
#fits well, then these will be roughly normally distributed 

# check goodness-of-fit of model:
# # 'qnorm(psres)', the one-step-ahead forecast pseudo-residuals, under the 
# model fitted to the calibration data - if the model fits well, then these will 
# be roughly normally distributed 
qqnorm(qnorm(psressvmt), main='SVM-t') 
qqline(qnorm(psressvmt))
require(tseries)
pvalue=jarque.bera.test(qnorm(psressvmt))$p.value
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
VaR=c(sum(psressvmt<0.01), 0.01*length(ytest))
VaR=rbind(VaR, c(sum(psressvmt<0.05), 0.05*length(ytest)))
VaR=rbind(VaR, c(sum(psressvmt<0.1), 0.1*length(ytest)))
VaR=rbind(VaR, c(sum(psressvmt<0.5), 0.5*length(ytest)))
row.names(VaR)=c('0.01', '0.05', '0.10', '0.50')
colnames(VaR)=c('Observed', 'Expected')
VaR
################################################################################
save(Results, times, h_hat, DIC, LPS, VaR, 
     file='~/Documentos/hmm/ibovespa/svmt.RData')


