################################################################################
library(parallel)
library(mvtnorm)
Rcpp::sourceCpp("~/Documentos/hmm/mlogLk_Rcpp.cpp")
Rcpp::sourceCpp("~/Documentos/hmm/pdf_s.cpp")
### Viterbi
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/HMMcodes/source/svm.viterbi.R')
### Diagnostic
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/HMMcodes/source/ISdiag.R')
################################################################################
num_cores <- detectCores(logical = FALSE) 
RcppParallel::setThreadOptions(numThreads = num_cores - 1) 
################################################################################
# Definir constantes para evitar overflow e valores extremos
double_eps <- .Machine$double.eps
double_xmax <- .Machine$double.xmax
double_xmin <- .Machine$double.xmin
log_double_xmax <- log(double_xmax)
log_double_xmin <- log(double_xmin)
################################################################################
svms.sim <-function(mu,phi,sigma,nu,beta,y0,g_dim){
  
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
  # nu > 0
  lnu = log(nu)
  parvect <- c(lbeta1,lbeta2,lbeta3,lmu,lphi,lsigma,lnu)
  return(parvect)
}
svm.pw2pn <- function(parvect){
  
  # Inicializar beta
  beta <- array(0, dim = 3)
  beta[1] <- parvect[1]
  if (parvect[2] > log_double_xmax) {
    beta[2] <- 1.0 - double_eps
  } else if (parvect[2] < -log_double_xmax) {
    beta[2] <- -1.0 + double_eps
  } else {
    beta[2] <- (exp(parvect[2]) - 1) / (exp(parvect[2]) + 1)
  }
  beta[3] <- parvect[3]
  
  mu <- parvect[4]
  if (parvect[5] > log_double_xmax) {
    phi <- 1.0 - double_eps
  } else if (parvect[5] < -log_double_xmax) {
    phi <- -1.0 + double_eps
  } else {
    phi <- (exp(parvect[5]) - 1) / (exp(parvect[5]) + 1)
  }
  
  if (parvect[6] > log_double_xmax) {
    sigma <- double_xmax
  } else if (parvect[6] < -log_double_xmax) {
    sigma <- double_xmin
  } else {
    sigma <- exp(parvect[6])
  }
  
  # Verificações para nu
  if (parvect[7] > log_double_xmax) {
    nu <- double_xmax
  } else if (parvect[7] < -log_double_xmax) {
    nu <-  double_eps
  } else {
    nu <- exp(parvect[7]) 
  }
  
  return(list(beta = beta, mu = mu, phi = phi, sigma = sigma, nu = nu))
}
## function that will be used to compute 'allprobs' in mllk below
fillallprobs <- function(x,beg,beta,nu,y){
  z = (x-beta[1]-beta[2]*y-beta[3]*beg^2)/beg
  w = pdf_s(y=z, nu=nu)
  #w = pdf_s(z,0.0,1.0,nu)
  return(w/beg)
}
svms.mllk <-function(parvect,y,y0,m,gmax){
  ny <- length(y)
  p <- svm.pw2pn(parvect)
  K = m+1
  b=seq(-gmax,gmax,length=K) 
  bs=(b[-1]+b[-K])*0.5
  E=p$mu+p$phi*(bs-p$mu)
  intlen <- b[2]-b[1]
  sey= exp(bs/2) 
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
svms.prior <-function(parvect){
  lprior = log(dnorm(parvect[1], 0, 10))
  + log(dnorm(parvect[2], 0.5, 10))
  + log(dnorm(parvect[3], 0, 10))
  + log(dnorm(parvect[4], 0, 10))
  + log(dnorm(parvect[5], 4.5, 10))
  + log(dnorm(parvect[6], -1.5, 10)) 
  + log(dnorm(parvect[7], 0, 10))
  return(-lprior)  
}
svms.posterior <-function(parvect,y,y0,m,gmax){
  return(svms.mllk(parvect,y,y0,m,gmax)+svms.prior(parvect))  
}
svms.map <- function(y, m, parvect0, y0, gmax){
  
  mod <- optim(parvect0, svms.posterior,y=y,y0=y0,m=m,gmax=gmax,hessian=T)
  mode <- mod$par
  conv = mod$convergence
  return(list(mode=mode,
              lpostsvs=-mod$value,
              hessian=mod$hessian,
              conv = conv))
}
################################################################################
mu0=0
phi0=0.98
sigma0 = 0.15
beta0=c(0.1 ,0.1 ,-0.05)
nu0=2
y0=0.2
################################################################################
m=200
gmax=5.0
y=log.ret

try({
  time = Sys.time()
  ############################################################################
  # Optimize
  eigenvalues = rep(0, 7)
  cont = 0
  parvect0 = svm.pn2pw(beta0, mu0, phi0, sigma0, nu0)
  while( any(eigenvalues <= 0) && cont<5 ){
    
    optim_res = svms.map(y=y,m=m, parvect0=parvect0, y0=y0, gmax=gmax)
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
    k = -optim_res$lpostsvs
    map = optim_res$mode
    ##########################################################################
    # Weigths Evaluation
    n = 1e3
    X=rmvnorm(n, map, H1)
    Weigth <- array(0,dim=c(n,1))
    for(j in 1:n){
      Weigth[j,1]=exp(k
                      -svms.posterior(parvect=X[j,],y=y,y0=y0,m=m,gmax=gmax)
                      -dmvnorm(X[j,],mean=map,sigma=H1,log=TRUE)
      )
    }
    s = sum(Weigth)
    if((s != 0) && !is.nan(s)){
      Weigth=Weigth/s  
    }else{
      stop('Error normalize constante weigths!')
    }
    Results = ISdiag(Weigth=Weigth, X=X, nu.lower=0)
    times = as.numeric(Sys.time()-time, units='mins')
  }
})

Results
times

h_hat=svm.viterbi(y=y,theta_hat=Results$Results[,1],m=m,gmax=gmax)
plot(abs(y), type='l', col='gray')
lines(exp(0.5*h_hat))

### -log p(y|\theta)
theta_hat=svm.pn2pw(beta=Results$Results[,1][1:3],
                    mu=Results$Results[,1][4],
                    phi=Results$Results[,1][5],
                    sigma=Results$Results[,1][6],
                    nu=Results$Results[,1][7])
svms.mllk(parvect=theta_hat, y=y, y0=y0, m=m, gmax=gmax)

