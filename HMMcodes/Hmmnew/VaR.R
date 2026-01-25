LRuc=function(u_t, alpha){

  v <- ifelse(u_t < alpha, 1, 0)
  N <- length(v)
  n1 <- sum(v)
  n0 <- N - n1
  pi_hat <- n1 / N


  L_ratio <- ( (alpha^n1) * (1 - alpha)^n0 )
  L_ratio=L_ratio/( (pi_hat^n1) * (1 - pi_hat)^n0 )

  lr_uc <- -2 * log(L_ratio)
  p_val <- 1 - pchisq(lr_uc, df = 1)

  return(list(stat=lr_uc, p_value=p_val))
}
LRind=function(u_t, alpha){
  v <- ifelse(u_t < alpha, 1, 0)
  N <- length(v)

  n01 <- sum(v[-N] == 0 & v[-1] == 1)
  n00 <- sum(v[-N] == 0 & v[-1] == 0)
  n11 <- sum(v[-N] == 1 & v[-1] == 1)
  n10 <- sum(v[-N] == 1 & v[-1] == 0)

  pi_01 <- n01 / (n00 + n01)
  pi_11 <- n11 / (n10 + n11)
  pi_all <- (n01 + n11) / (n00 + n01 + n10 + n11)

  L_ratio <- ( (1 - pi_all)^(n00 + n10) * pi_all^(n01 + n11) )
  L_ratio = L_ratio/( (1 - pi_01)^n00 * pi_01^n01 * (1 - pi_11)^n10 * pi_11^n11 )
  lr_ind <- -2*log(L_ratio)

  p_val <- 1-pchisq(lr_ind, df=1)

  return(list(stat=lr_ind, p_value=p_val))
}

fpr=svmn_fpr
fpr=svmt_fpr
fpr=svms_fpr
fpr=svmvg_fpr

par(mfrow=c(1, 2))
qqnorm(qnorm(fpr), main='SVM-N')
qqline(qnorm(fpr))
acf(qnorm(fpr))
par(mfrow=c(1, 1))

require(tseries)
pvalue=jarque.bera.test(qnorm(fpr))$p.value
if(pvalue>0.05){
  cat('Not reject the hypothesis of normality of the residuals at the 5% nivel.')
}else{
  cat('Reject the hypothesis of normality of the residuals at the 5% nivel.')
}
pvalue

# now backtesting
# as an example, check how many of the observations in the validation sample are
# smaller than the 1%/5%/10%/50% quantiles of the one-step-ahead forecast
# distributions, under the fitted model
# should be approximately 0.01*length(obs.val)
VaR=c(sum(fpr<0.01), 0.01*length(ytest))
VaR=rbind(VaR, c(sum(fpr<0.05), 0.05*length(ytest)))
#VaR=rbind(VaR, c(sum(fpr<0.1), 0.1*length(ytest)))
#VaR=rbind(VaR, c(sum(fpr<0.5), 0.5*length(ytest)))
row.names(VaR)=c('0.01', '0.05')
colnames(VaR)=c('Observed', 'Expected')
VaR/T

LRuc(u_t=fpr, alpha=0.01)
LRind(u_t=fpr, alpha=0.01)
1-pchisq(LRuc(u_t=fpr, alpha=0.01)$stat+LRind(u_t=fpr, alpha=0.01)$stat, df=2)


LRuc(u_t=fpr, alpha=0.05)
LRind(u_t=fpr, alpha=0.05)
1-pchisq(LRuc(u_t=fpr, alpha=0.05)$stat+LRind(u_t=fpr, alpha=0.05)$stat, df=2)
