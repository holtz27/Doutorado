################################################################################
rm(list=ls(all=TRUE))
library("Rcpp")
library("RcppArmadillo")
library(parallel)
library(mvtnorm)
sourceCpp("mlogLk_Rcpp.cpp")
sourceCpp("parallel_pdf_s.cpp")
#sourceCpp("svsl_map_Rcpp.cpp")
################################################################################
num_cores = detectCores() 
RcppParallel::setThreadOptions(numThreads = num_cores - 1) 
################################################################################
ssvm.sim<-function(mu,phi,sigma,nu,beta,y0,g_dim){
  
  y=array(0,dim=g_dim)  
  h=array(0,dim=g_dim)    
  l=array(0,dim=g_dim)
  l=rbeta(g_dim,nu,1,ncp=0)
  
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
  # nu > 1
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
## function that will be used to compute 'allprobs' in mllk below
fillallprobs <- function(x,beg,beta,nu,y){
  return((1/beg)*pdf_s((x-beta[1]-beta[2]*y-beta[3]*beg^2)/beg,0.0,1.0,nu))
}
svms.mllk <-function(parvect,y,y0,m,gmax){
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
svms.prior <-function(parvect){
  lprior=log(dnorm(parvect[1],0,10))+log(dnorm(parvect[2],0.5,10))
  +log(dnorm(parvect[3],0,10))+log(dnorm(parvect[4],0,10))+log(dnorm(parvect[5],4.5,10))
  +log(dnorm(parvect[6],-1.5,10)) + log(dnorm(parvect[7],0,10))
  return(-lprior)  
}
svms.posterior <-function(parvect,y,y0,m,gmax){
  return(svms.mllk(parvect,y,y0,m,gmax)+svms.prior(parvect))  
}
svsl.map <- function(y,m,mu0,phi0,sigma0,nu0,beta0,y0,gmax){
  time = Sys.time()
  parvect <- svm.pn2pw(beta=beta0,mu=mu0,phi=phi0,sigma=sigma0,nu=nu0)
  mod <- optim(parvect,svms.posterior,y=y,y0=y0,m=m,gmax=gmax,hessian=T)
  mode <- mod$par
  time = Sys.time()-time
  return(list(mode=mode,lpostsvs=-mod$value,hessian=mod$hessian,time=time))
}
quantile <- function(x, weights, probs=c(0.025, 0.975)){
  
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
################################################################################
mu=0.1
phi=0.98
sigma=0.1
beta=c(0.2,0.07,-0.18)
nu=2
y0=0.2
g_dim=2000

y = ssvm.sim(mu,phi,sigma,nu,beta,y0,g_dim)
plot(y$y)

################################################################################
mu0=0
phi0=0.99
sigma0 = 0.2
beta0=c(0,0.1,-0.1)
nu0=5
parvect0 = c(beta0, mu0, phi0, sigma0, nu0)
m=50
gmax=2
################################################################################
Weights = function(parvect0, y, y0, n, m, gmax){
  
  res = svsl.map(y, m=m, y0=y0, gmax=gmax,
                 beta0=parvect0[1:3],
                 mu0=parvect0[4],
                 phi0=parvect0[5],
                 sigma0=parvect0[6],
                 nu0=parvect0[7]
                 )
  k = -res$lpostsvs
  H1=solve(res$hessian)
  mode=res$mode
  X=rmvnorm(n,mean=mode,sigma=H1)
  Weigth<-array(0,dim=c(n,1))
  for(j in 1:n){
    Weigth[j,1]=exp(k
                    -svms.posterior(X[j,],y,y0=y0,m=m,gmax=gmax)
                    -log(dmvnorm(X[j,],mean=mode,sigma=H1))
                    )
  }
  return( list('Weigth' = Weigth/sum(Weigth), 'X'=X ) )
}

w = Weights(parvect0 = parvect0, n=500, y=y$y, y0=y0, m=m, gmax=gmax)
w
sourceCpp('Weigths_cpp.cpp')
w2 = Weigth_cpp(parvect_init=parvect0, n=500, y=y$y, y0=y0, m=m, gmax=gmax)
w2
################################################################################
library(microbenchmark)

benchmarks = microbenchmark(
  weigths_R = Weights(parvect0=parvect0, n=500, y=y$y, y0=y0, m=m, gmax=gmax),
  weigths_Rcpp = Weigth_cpp(parvect_init=parvect0, n=500, y=y$y, y0=y0, m=m, gmax=gmax),
  times = 10 
)

print(benchmarks)
boxplot(benchmarks)
################################################################################
Weigth=w2$weigths
X=w2$X
beta0m=sum(X[1,]*Weigth)
beta1m=sum(((exp(X[2,])-1)/(exp(X[2,])+1))*Weigth)
beta2m=sum(X[3,]*Weigth)
mum=sum(X[4,]*Weigth)
phim=sum(((exp(X[5,])-1)/(exp(X[5,])+1))*Weigth)
sigmam=sum(exp(X[6,])*Weigth)
num=sum((exp(X[7,]))*Weigth)
################################################################################
int_b0 = quantile(X[1,], Weigth)
int_b1 = quantile((exp(X[2,])-1)/(exp(X[2,])+1), Weigth)
int_b2 = quantile(X[3,], Weigth)
int_mu = quantile(X[4,], Weigth)
int_phi = quantile( (exp(X[5,])-1)/(exp(X[5,])+1), Weigth)
int_sigma = quantile(exp(X[6,]), Weigth)
int_nu = quantile((exp(X[7,])), Weigth)
################################################################################
# Results
results = matrix(c(beta[1], beta0m, int_b0),nrow = 1)
results = rbind(results,
                c(beta[2], beta1m, int_b1),
                c(beta[3], beta2m, int_b2),
                c(mu, mum, int_mu),
                c(phi, phim, int_phi),
                c(sigma, sigmam, int_sigma),
                c(nu, num, int_nu))
row.names(results) = c('b0','b1','b2','mu','phi','sigma','nu')
colnames(results) = c('true','mean','2.5%','97.5%')
results
