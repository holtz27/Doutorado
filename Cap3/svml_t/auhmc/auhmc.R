source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/num_analisys.R')
source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R')
path='~/Documentos/tese/Cap3/svml_t/auhmc/Rcpp_auhmc.cpp'
Rcpp::sourceCpp(path)
# ==============================================================================
# 1. Simulação dos Dados (Modelo SV-t com Leverage)
# ==============================================================================

#set.seed(12345)

b0 <- 0.1; 
b1 <- 0.01; 
b2 <- -0.05
mu <- 0.1; 
phi_true <- 0.98; 
sigma_true <- 0.15
nu_true <- 10; 
rho_true <- -0.45
delta_true <- sigma_true * rho_true          # -0.045
tau_true <- sigma_true * sqrt(1 - rho_true^2) # ~0.143
y0 <- 0;

theta_true=c(b0, b1, b2, mu, phi_true, delta_true, tau_true, nu_true)

N <- 5000
y <- h <- l <- numeric(N)
l[1] <- rgamma(1, shape = 0.5 * nu_true, rate = 0.5 * nu_true)
h[1] <- mu + (sigma_true / sqrt(1 - phi_true^2)) * rnorm(1)
y[1] <- b0 + b1 * y0 + b2 * exp(h[1]) + (exp(0.5 * h[1]) / sqrt(l[1])) * rnorm(1)
for (t in 2:N) {
  mu_prev    <- b0 + b1 * (if (t == 2) y0 else y[t-2]) + b2 * exp(h[t-1])
  sigma_prev <- exp(0.5 * h[t-1]) / sqrt(l[t-1])
  zt_prev    <- (y[t-1] - mu_prev) / sigma_prev
  
  l[t] <- rgamma(1, shape = 0.5 * nu_true, rate = 0.5 * nu_true)
  h[t] <- mu + phi_true * (h[t-1] - mu) + delta_true * zt_prev + tau_true * rnorm(1)
  y[t] <- b0 + b1 * y[t-1] + b2 * exp(h[t]) + (exp(0.5 * h[t]) / sqrt(l[t])) * rnorm(1)
}

par(mfrow=c(1,2))
plot(y, type='l', col='gray')
lines(exp(0.5*h), col='firebrick',lwd=2)
hist(y, breaks=40, freq=FALSE)
par(mfrow=c(1,1))

############################ drawns ##################################
#y = log.ret
#N = length(y)
# N_total = burnin + it * thin
N_total=3e2
burnin=1e2
thin=1
n_draws = (N_total - burnin) / thin
n_draws

drawns=auhmc_svml(y=c(0, y), 
                  n_draws=n_draws, burnin=burnin, thin=thin,
                  nlpfrog=c(10, 10, 10), 
                  epsilon=c(0.01, 0.05, 0.05))

drawns$theta_mean_j_iters
drawns$theta_non_conv_rate

drawns$h_mean_j_iters
drawns$h_non_conv_rate

drawns$l_mean_j_iters
drawns$l_non_conv_rate
################### Numeric Analysis
drawn_theta=drawns$theta
Res=num_analisys(drawn_theta,
                 names=c('b0','b1','b2','mu','phi','delta','tau','v'),
                 digits=4)
#Res
cbind(true=round(theta_true, 4), Res)

################### Plots
# Trace plot
trace_plots(drawn_theta, 
            names=c(expression(beta[0]),
                    expression(beta[1]),
                    expression(beta[2]),
                    expression(mu),
                    expression(phi),
                    expression(delta),
                    expression(tau),
                    expression(nu))
)


drawn_h=drawns$h
min_h=apply(drawn_h, 1, quantile, probs=c(0.025))
max_h=apply(drawn_h, 1, quantile, probs=c(0.975))
par(mfrow=c(2,1))
plot(rowMeans(drawn_h[1:(0.5*N),]), type='l', xlab='', ylab='', main='',
     ylim=c(min(min_h), max(max_h)))
lines(min_h[1:(0.5*N)], lty=2, col='gray50')
lines(max_h[1:(0.5*N)], lty=2, col='gray50')
lines(h[1:(0.5*N)], lty=3, col='red',lwd=2)
plot(rowMeans(drawn_h[(0.5*N+1):N,]), type='l', , xlab='', ylab='', main='',
     ylim=c(min(min_h), max(max_h)))
lines(min_h[(0.5*N+1):N], lty=2, col='gray50')
lines(max_h[(0.5*N+1):N], lty=2, col='gray50')
lines(h[(0.5*N+1):N], lty=3, col='red',lwd=2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
plot(y, type='l', col='gray')
lines(exp(0.5*rowMeans(drawn_h)), col='firebrick',lwd=2)
drawn_l=drawns$l
plot(rowMeans(drawn_l), type='l')
par(mfrow=c(1,1))

