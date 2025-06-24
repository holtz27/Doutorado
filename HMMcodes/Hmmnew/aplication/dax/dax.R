dax = quantmod::getSymbols('^GDAXI', 
                            src='yahoo', 
                            from='2002-01-01', to='2025-01-01',
                            auto.assign=FALSE)
dax = na.omit(dax)
dax = data.frame(dax)
dates = as.Date(row.names(dax), '%Y-%m-%d')
dax = dax[, 'GDAXI.Adjusted']

T = length(dax)
log.ret = 100*(log(dax[2:T])-log(dax[1:(T-1)]))
T = length(log.ret)

h=1e3
ytrain=log.ret[1:(length(log.ret)-h)]
ytest=log.ret[(length(log.ret)-h+1):length(log.ret)]

# Plots
library(ggplot2)
df = data.frame(Return=ytrain, Tempo=dates[1:length(ytrain)])

g = ggplot(df) + geom_line(aes(x = Tempo, y = Return))
g = g + scale_x_date(date_breaks = "28 month", date_labels = "%b %Y")
g = g + theme_test() + theme(axis.title.y = element_text(size = 18),
                             axis.text.x = element_text(size = 12),
                             axis.text.y = element_text(size = 18))
g = g + xlab('')
h = ggplot( df, aes(Return) )
h = h + geom_histogram(aes(y = after_stat(density)), bins = 40, color = 'white')
h = h + theme_test() + ylab('')
h = h + theme_test() + theme(axis.title.x = element_text(size = 18),
                             axis.text.x = element_text(size = 18),
                             axis.text.y = element_text(size = 18))
gridExtra::grid.arrange(g, h, nrow = 1, ncol = 2) 
data_summary=matrix(c(T, mean(ytrain),
                      sd(ytrain),
                      min(ytrain),
                      max(ytrain),
                      moments::skewness(ytrain),
                      moments::kurtosis(ytrain)), nrow=1)
colnames(data_summary)=c('T','mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round(data_summary, digits=3)
################################################################################
####### Fitting
################################################################################
# SVM-N
source('~/Documentos/svmHMM/aplic_n.R')
Results
save(Results, times, h_hat, DIC, LPS, file='~/Documentos/svmHMM/dax/n.RData')

# SVM-t
source('~/Documentos/svmHMM/aplic_t.R')
Results
save(Results, times, h_hat, DIC, LPS, file='~/Documentos/svmHMM/dax/t.RData')

# SVM-S
source('~/Documentos/svmHMM/aplic_s.R')
Results
save(Results, times, h_hat, DIC, LPS, file='~/Documentos/svmHMM/dax/s.RData')

# SVM-VG
source('~/Documentos/svmHMM/aplic_vg.R')
Results
save(Results, times, h_hat, DIC, LPS, file='~/Documentos/svmHMM/dax/vg.RData')
###############################################################################
### Analysis Results
################################################################################
# SVM-N
load("~/HMMnew/aplication/dax/n.RData")
Results
times
DIC
LPS

#plot(abs(ytrain), type='l', col='gray')
#lines(exp(0.5*h_hat), lwd=2)
nh_hat=h_hat

qqnorm(qnorm(psressvmn), main='SVM-N', cex.main=2.5, cex.axis=1.8, cex.lab=1.8) 
qqline(qnorm(psressvmn))
require(tseries)
pvalue=jarque.bera.test(qnorm(psressvmn))$p.value
if(pvalue>0.05){
  cat('Not reject the hypothesis of normality of the residuals at the 5% nivel.', '\n')
}else{
  cat('Reject the hypothesis of normality of the residuals at the 5% nivel.', '\n')
}

VaR=c(sum(psressvmn<0.01), 0.01*length(ytest))
VaR=rbind(VaR, c(sum(psressvmn<0.05), 0.05*length(ytest)))
VaR=rbind(VaR, c(sum(psressvmn<0.1), 0.1*length(ytest)))
VaR=rbind(VaR, c(sum(psressvmn<0.5), 0.5*length(ytest)))
row.names(VaR)=c('0.01', '0.05', '0.10', '0.50')
colnames(VaR)=c('Observed', 'Expected')
VaR
###############################################################################
# Unconditional Coverage Test 
n=length(ytest)
x=VaR[1,1]
alpha=0.01
alpha_hat=x/n
LRuc=log(alpha_hat^{x}*(1-alpha_hat)^{n-x})
LRuc=LRuc-log(alpha^{x}*(1-alpha)^{n-x})
LRuc=2*LRuc
pvalue001=pchisq(LRuc, df=1, lower.tail=FALSE)


x=VaR[2,1]
alpha=0.05
alpha_hat=x/n
LRuc=log(alpha_hat^{x}*(1-alpha_hat)^{n-x})
LRuc=LRuc-log(alpha^{x}*(1-alpha)^{n-x})
LRuc=2*LRuc
pvalue005=pchisq(LRuc, df=1, lower.tail=FALSE)

round(pvalue001, 4)
round(pvalue005, 4)
################################################################################
# SVM-t
load("~/HMMnew/aplication/dax/t.RData")
Results
times
DIC
LPS
#plot(abs(ytrain), type='l', col='gray')
#lines(exp(0.5*h_hat), lwd=2)
th_hat=h_hat

qqnorm(qnorm(psressvmt), main='SVM-t', cex.main=2.5, cex.axis=1.8, cex.lab=1.8) 
qqline(qnorm(psressvmt))
require(tseries)
pvalue=jarque.bera.test(qnorm(psressvmt))$p.value
if(pvalue>0.05){
  cat('Not reject the hypothesis of normality of the residuals at the 5% nivel.', '\n')
}else{
  cat('Reject the hypothesis of normality of the residuals at the 5% nivel.', '\n')
}

VaR=c(sum(psressvmn<0.01), 0.01*length(ytest))
VaR=rbind(VaR, c(sum(psressvmt<0.05), 0.05*length(ytest)))
VaR=rbind(VaR, c(sum(psressvmt<0.1), 0.1*length(ytest)))
VaR=rbind(VaR, c(sum(psressvmt<0.5), 0.5*length(ytest)))
row.names(VaR)=c('0.01', '0.05', '0.10', '0.50')
colnames(VaR)=c('Observed', 'Expected')
VaR
###############################################################################
# Unconditional Coverage Test 
n=length(ytest)
x=VaR[1,1]
alpha=0.01
alpha_hat=x/n
LRuc=log(alpha_hat^{x}*(1-alpha_hat)^{n-x})
LRuc=LRuc-log(alpha^{x}*(1-alpha)^{n-x})
LRuc=2*LRuc
pvalue001=pchisq(LRuc, df=1, lower.tail=FALSE)


x=VaR[2,1]
alpha=0.05
alpha_hat=x/n
LRuc=log(alpha_hat^{x}*(1-alpha_hat)^{n-x})
LRuc=LRuc-log(alpha^{x}*(1-alpha)^{n-x})
LRuc=2*LRuc
pvalue005=pchisq(LRuc, df=1, lower.tail=FALSE)

round(pvalue001, 4)
round(pvalue005, 4)
################################################################################
# SVM-S
load("~/HMMnew/aplication/dax/s.RData")
Results
times
DIC
LPS
#plot(abs(ytrain), type='l', col='gray')
#lines(exp(0.5*h_hat), lwd=2)
sh_hat=h_hat

qqnorm(qnorm(psressvms), main='SVM-S', cex.main=2.5, cex.axis=1.8, cex.lab=1.8) 
qqline(qnorm(psressvms))
require(tseries)
pvalue=jarque.bera.test(qnorm(psressvms))$p.value
if(pvalue>0.05){
  cat('Not reject the hypothesis of normality of the residuals at the 5% nivel.', '\n')
}else{
  cat('Reject the hypothesis of normality of the residuals at the 5% nivel.', '\n')
}

VaR=c(sum(psressvmn<0.01), 0.01*length(ytest))
VaR=rbind(VaR, c(sum(psressvms<0.05), 0.05*length(ytest)))
VaR=rbind(VaR, c(sum(psressvms<0.1), 0.1*length(ytest)))
VaR=rbind(VaR, c(sum(psressvms<0.5), 0.5*length(ytest)))
row.names(VaR)=c('0.01', '0.05', '0.10', '0.50')
colnames(VaR)=c('Observed', 'Expected')
VaR
###############################################################################
# Unconditional Coverage Test 
n=length(ytest)
x=VaR[1,1]
alpha=0.01
alpha_hat=x/n
LRuc=log(alpha_hat^{x}*(1-alpha_hat)^{n-x})
LRuc=LRuc-log(alpha^{x}*(1-alpha)^{n-x})
LRuc=2*LRuc
pvalue001=pchisq(LRuc, df=1, lower.tail=FALSE)


x=VaR[2,1]
alpha=0.05
alpha_hat=x/n
LRuc=log(alpha_hat^{x}*(1-alpha_hat)^{n-x})
LRuc=LRuc-log(alpha^{x}*(1-alpha)^{n-x})
LRuc=2*LRuc
pvalue005=pchisq(LRuc, df=1, lower.tail=FALSE)

round(pvalue001, 4)
round(pvalue005, 4)
################################################################################
# SVM-VG
load("~/HMMnew/aplication/dax/vg.RData")
Results
times
DIC
LPS
#plot(abs(ytrain), type='l', col='gray')
#lines(exp(0.5*h_hat), lwd=2)
vgh_hat=h_hat

qqnorm(qnorm(psressvmvg), main='SVM-VG', cex.main=2.5, cex.axis=1.8, cex.lab=1.8) 
qqline(qnorm(psressvmvg))
require(tseries)
pvalue=jarque.bera.test(qnorm(psressvmvg))$p.value
if(pvalue>0.05){
  cat('Not reject the hypothesis of normality of the residuals at the 5% nivel.', '\n')
}else{
  cat('Reject the hypothesis of normality of the residuals at the 5% nivel.', '\n')
}

VaR=c(sum(psressvmn<0.01), 0.01*length(ytest))
VaR=rbind(VaR, c(sum(psressvmvg<0.05), 0.05*length(ytest)))
VaR=rbind(VaR, c(sum(psressvmvg<0.1), 0.1*length(ytest)))
VaR=rbind(VaR, c(sum(psressvmvg<0.5), 0.5*length(ytest)))
row.names(VaR)=c('0.01', '0.05', '0.10', '0.50')
colnames(VaR)=c('Observed', 'Expected')
VaR
###############################################################################
# Unconditional Coverage Test 
n=length(ytest)
x=VaR[1,1]
alpha=0.01
alpha_hat=x/n
LRuc=log(alpha_hat^{x}*(1-alpha_hat)^{n-x})
LRuc=LRuc-log(alpha^{x}*(1-alpha)^{n-x})
LRuc=2*LRuc
pvalue001=pchisq(LRuc, df=1, lower.tail=FALSE)


x=VaR[2,1]
alpha=0.05
alpha_hat=x/n
LRuc=log(alpha_hat^{x}*(1-alpha_hat)^{n-x})
LRuc=LRuc-log(alpha^{x}*(1-alpha)^{n-x})
LRuc=2*LRuc
pvalue005=pchisq(LRuc, df=1, lower.tail=FALSE)

round(pvalue001, 4)
round(pvalue005, 4)


plot(abs(ytrain), type='l', col='gray', 
     main='DAX 30', xlab='', ylab='', cex.main=2.5, cex.axis=1.8, cex.lab=1.8)
lines(exp(0.5*nh_hat), col='black', lwd=2)
lines(exp(0.5*th_hat), col='red', lwd=2)
lines(exp(0.5*sh_hat), col='orange', lwd=2)
lines(exp(0.5*vgh_hat), col='purple', lwd=2)
