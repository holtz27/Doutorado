rm(list=ls(all=TRUE))
library("Rcpp")
library("RcppArmadillo")
library('mvtnorm')
sourceCpp("mlogLk_Rcpp.cpp")
###############################################################################
svm.pn2pw <- function(beta,mu,phi,sigma,nu){
  lbeta1<- beta[1]
  lbeta2<-log((1+beta[2])/(1-beta[2]))
  lbeta3<-beta[3]
  lmu<-mu
  lphi <- log((1+phi)/(1-phi))
  lsigma <- log(sigma)
  # nu > 2
  lnu = log(nu-2)
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
  nu = exp(parvect[7]) + 2
  return(list(beta=beta,mu=mu,phi=phi,sigma=sigma,nu=nu))
}
fillallprobs <- function(x,beg,beta,nu,y){
  return( (1/beg)*dt((x-beta[1]-beta[2]*y-beta[3]*beg^2)/beg,df=nu) )
}
svmt.mllk <-function(parvect,y,y0,m,gmax){
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
  allprobs <- outer(xx,sey,"fillallprobs",beta=p$beta,nu=p$nu,yy)
  delta <-dnorm(bs,p$mu,p$sigma/sqrt(1-p$phi^2))*intlen
  foo <- delta*allprobs[1,]
  lscale = mlogLk_Rcpp(allprobs,Gamma,foo,ny) #Rcpp function
  return(-lscale)
}
svmt.prior <-function(parvect){
  lprior=log(dnorm(parvect[1],0,10))+log(dnorm(parvect[2],0.5,10))
  +log(dnorm(parvect[3],0,10))+log(dnorm(parvect[4],0,10))+log(dnorm(parvect[5],4.5,10))
  +log(dnorm(parvect[6],-1.5,10)) + log(dnorm(parvect[7],0,10))
  return(-lprior)  
}
svmt.posterior <-function(parvect,y,y0,m,gmax){
  return(svmt.mllk(parvect,y,y0,m,gmax)+svmt.prior(parvect))  
}
svmt.map <- function(y,y0,m,beta0,mu0,phi0,sigma0,nu0,gmax){
  time = Sys.time()
  parvect <- svm.pn2pw(beta=beta0,mu=mu0,phi=phi0,sigma=sigma0,nu=nu0)
  mod <- optim(parvect,svmt.posterior,y=y,y0=y0,m=m,gmax=gmax,hessian=T)
  mode <- mod$par
  time = Sys.time()-time
  return(list(mode=mode,lpostsvmt=-mod$value,hessian=mod$hessian,time=time))
}
tsvmeansim <-function(y0,beta,mu,phi,sigma,nu,g_dim){
  y=array(0,dim=g_dim)  
  h=array(0,dim=g_dim)    
  h[1]=rnorm(1,mu,sigma/sqrt(1-phi^2))
  y[1]=beta[1]+beta[2]*y0+beta[3]*exp(h[1])+exp(h[1]/2)*rt(1,df=nu,ncp=0)#rnorm(1,0,1)
  for (j in 2:g_dim){
    h[j]=mu+phi*(h[j-1]-mu)+rnorm(1,0,sigma)
    y[j]=beta[1]+beta[2]*y[j-1]+beta[3]*exp(h[j])+exp(h[j]/2)*rt(1,df=nu,ncp=0)#rnorm(1,0,1)
  }
  return (list(y=y,h=h,y0=y0))
}
quantile <- function(x, weights, probs) {
 
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
###############################################################################
betasim = c(0.2,0.07,-0.18)
mu = 0.1
phi = 0.98
sigma = 0.1
nu = 8
g_dim = 6000
y0 = 0.2
set.seed(8936381)
simsvmt1<-tsvmeansim(y0=y0,beta=betasim,mu=mu,phi=phi,sigma=sigma,nu=nu,g_dim=g_dim)
plot(simsvmt1$y,lty=1,type="l")
hist(simsvmt1$y,breaks = 40)
plot(simsvmt1$h,lty=1,type="l")
###############################################################################
# theta.init
mu0=0.3
phi0=0.96
sigma0=0.2
beta0=c(0.2, 0.1, -0.1)
nu0=10
gmax=4

s1 = 1000
#s2=3000
#s3=6000
######################################################################
m=100
res.svmt200_4_s1 <- svmt.map(simsvmt1$y[1:s1],m=m,
                             beta0=beta0,mu0=mu0,phi0=phi0,
                             sigma0=sigma0,nu0=nu0,y0=simsvmt1$y0,
                             gmax=gmax)
res.svmt200_4_s1
k = -res.svmt200_4_s1$lpostsvmt
#################################################################
H1=signif(solve(res.svmt200_4_s1$hessian),6)
#res.svmt200_4_s1$mode
n = 500
X=rmvnorm(n,res.svmt200_4_s1$mode,H1)
Weigth <- array(0,dim=c(n,1))
for(j in 1:n){
  if(j==1) time = Sys.time()
  Weigth[j,1]=exp(k #5200
                  -svmt.posterior(X[j,],simsvmt1$y[1:s1],
                                  simsvmt1$y0,m=m,gmax=gmax)
                  -dmvnorm(X[j,],res.svmt200_4_s1$mode,sigma=H1,log=T)
                  )
  if(j==n) time = Sys.time()-time
}
time
################################################################
Weigth=Weigth/sum(Weigth)
beta0m=sum(X[,1]*Weigth)
beta1m=sum(((exp(X[,2])-1)/(exp(X[,2])+1))*Weigth)
beta2m=sum(X[,3]*Weigth)
mum=sum(X[,4]*Weigth)
phim=sum(((exp(X[,5])-1)/(exp(X[,5])+1))*Weigth)
sigmam=sum(exp(X[,6])*Weigth)
num=sum((exp(X[,7])+2)*Weigth)
################################################################
cred_int_b0 = quantile(X[,1], Weigth, c(0.025, 0.975))
cred_int_b1 = quantile((exp(X[,2])-1)/(exp(X[,2])+1), 
                                Weigth, c(0.025, 0.975))
cred_int_b2 = quantile(X[,3], Weigth, c(0.025, 0.975))
cred_int_mu = quantile(X[,4], Weigth, c(0.025, 0.975))
cred_int_phi = quantile( (exp(X[,5])-1)/(exp(X[,5])+1), 
                        Weigth, c(0.025, 0.975))
cred_int_sigma = quantile(exp(X[,6]), 
                         Weigth, c(0.025, 0.975))
cred_int_nu = quantile( (exp(X[,7])+2), 
                       Weigth, c(0.025, 0.975))

# Results
results = matrix(c(betasim[1], beta0m, cred_int_b0),nrow = 1)
results = rbind(results,
                c(betasim[2], beta1m, cred_int_b1),
                c(betasim[3], beta2m, cred_int_b2),
                c(mu, mum, cred_int_mu),
                c(phi, phim, cred_int_phi),
                c(sigma, sigmam, cred_int_sigma),
                c(nu, num, cred_int_nu))
row.names(results) = c('b0','b1','b2','mu','phi','sigma','nu')
colnames(results) = c('true','mean','2.5%','97.5%')
round(results, 4)



vbeta0m<-sum(X[,1]^2*Weigth)-beta0m^2
vbeta1m<-sum((((exp(X[,2])-1)/(exp(X[,2])+1)))^2*Weigth)-beta1m^2
vbeta2m<-sum(X[,3]^2*Weigth)-beta2m^2
vmum<-sum(X[,4]^2*Weigth)-mum^2
vphim<-sum((((exp(X[,5])-1)/(exp(X[,5])+1)))^2*Weigth)-phim^2
vsigmam<-sum(exp(2*X[,6])*Weigth)-sigmam^2
vnum = sum((exp(X[,7])+2)^2*Weigth)-num^2
#########################################################
