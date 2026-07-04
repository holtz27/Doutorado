source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R')

# compiling stan models
path='~/tese/SVMHeavy/rstan/svmt.stan'
model_stan = rstan::stan_model(file=path)

###################
## fitting model ##
###################
warmup=5e2
iter=2e2


draws=rstan::sampling(model_stan, 
                      data=list(T=length(log.ret), 
                                y=as.numeric(log.ret),
                                y0=0.2),
                      chains=1,
                      warmup=warmup,
                      iter=warmup+iter,
                      cores=1)
pars=c('mu','phi','s_h','b0','b1','b2','v', 'h','l')
x=rstan::extract(draws, pars=pars)
theta=rbind(x$mu, x$ph, x$s_h, x$b0, x$b1, x$b2, x$v)
names=c('mu','phi','sh','b0','b1','b2','v')

trace_plots(theta, names=names)

print(draws, 
      pars=c('mu','phi','s_h','b0','b1','b2','v'), 
      probs=c(.025,.975))
