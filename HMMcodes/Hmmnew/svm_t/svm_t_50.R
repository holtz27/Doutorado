library('mvtnorm')
library(parallel)
Rcpp::sourceCpp("~/mlogLk_Rcpp.cpp")
Rcpp::sourceCpp("~/pdf_t.cpp")
################################################################################
num_cores <- detectCores(logical = FALSE) 
RcppParallel::setThreadOptions(numThreads = num_cores - 1) 
################################################################################
double_eps <- .Machine$double.eps
double_xmax <- .Machine$double.xmax
double_xmin <- .Machine$double.xmin
log_double_xmax <- log(double_xmax)
log_double_xmin <- log(double_xmin)
################################################################################
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
fillallprobs <- function(x,beg,beta,nu,y){
  return( (1/beg)*pdf_t((x-beta[1]-beta[2]*y-beta[3]*beg^2)/beg,df=nu) )
  #return( (1/beg)*dt((x-beta[1]-beta[2]*y-beta[3]*beg^2)/beg,df=nu) )
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
  lprior=dnorm(parvect[1],0,10, log = TRUE)
  +dnorm(parvect[2],0.5,10, log = TRUE)
  +dnorm(parvect[3],0,10, log = TRUE)
  +dnorm(parvect[4],0,10, log = TRUE)
  +dnorm(parvect[5],4.5,10, log = TRUE)
  +dnorm(parvect[6],-1.5,10, log = TRUE) 
  +dnorm(parvect[7],0, 10, log = TRUE)
  return(-lprior)  
}
svmt.posterior <-function(parvect,y,y0,m,gmax){
  return(svmt.mllk(parvect,y,y0,m,gmax)+svmt.prior(parvect))  
}
svmt.map <- function(y, y0, m, parvect0, gmax){
  
  #parvect <- svm.pn2pw(beta=beta0,mu=mu0,phi=phi0,sigma=sigma0,nu=nu0)
  mod <- optim(parvect0,svmt.posterior,y=y,y0=y0,m=m,gmax=gmax, hessian=TRUE)
  mode <- mod$par
  conv = mod$convergence
  
  return(list(mode=mode,
              lpostsvmt=-mod$value,
              hessian=mod$hessian,
              conv = conv))
}
svmt.sim <-function(y0,beta,mu,phi,sigma,nu,g_dim){
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
### Diagnostic
source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/HMMcodes/source/ISdiag.R')
################################################################################
mu=0.1
phi=0.98
sigma=0.1
beta=c(0.2,0.07,-0.18)
nu=10
theta = c(beta,mu,phi,sigma,nu)
################################################################################
mu0=0.5
phi0=0.99
sigma0 = 0.2
beta0=c(0.5,0.1,-0.1)
nu0=8
y0=0.2
################################################################################
g_dim=1500
m=50
gmax=2.5
################################################################################

N=300
reps=w=sim.y=list()
times = npd = rep(0, N)

for(rep in 1:N){
  
  try({
    time = Sys.time()
    ############################################################################
    # Data simulation
    y = svmt.sim(y0=y0,beta=beta,mu=mu,phi=phi,sigma=sigma,nu=nu,g_dim=g_dim)
    ############################################################################
    # Optimize
    eigenvalues = rep(0, 7)
    cont = 0
    parvect0 = svm.pn2pw(beta0, mu0, phi0, sigma0, nu0)
    while( any(eigenvalues <= 0) && cont<5 ){
      
      optim_res = svmt.map(y=y$y,m=m, parvect0=parvect0, y0=y0, gmax=gmax)
      H1 = signif(solve(optim_res$hessian), 6)
      eigenvalues = eigen(H1, only.values=TRUE)$values
      parvect0 = optim_res$mode
      cont = cont + 1
      
    }
    ############################################################################
    # Test if H1 is positive definite
    if(any(eigenvalues <= 0)){
      npd[rep] = 1
      warning('Hessian is not positive definite')
    }else{
      k = - optim_res$lpostsvmt
      map = optim_res$mode
      ##########################################################################
      # Weigths Evaluation
      n = 1e3
      X=rmvnorm(n, map, H1)
      Weigth <- array(0,dim=c(n,1))
      for(j in 1:n){
        Weigth[j,1]=exp(k
                        -svmt.posterior(parvect=X[j,],y=y$y,y0=y0,m=m,gmax=gmax)
                        -dmvnorm(X[j,],mean=map,sigma=H1,log=TRUE)
        )
      }
      s = sum(Weigth)
      if((s != 0) && !is.nan(s)){
        Weigth=Weigth/s  
      }else{
        warning('Error normalize constante weigths!')
      }
      reps[[rep]] = ISdiag(Weigth=Weigth, X=X, knu=2)
    }
    w[[rep]] = Weigth
    sim.y[[rep]] = y$y
    times[rep] = as.numeric(Sys.time()-time, units='mins')
  })
  
}

save(reps, times, w, sim.y, file = 't_1500_50.RData')
