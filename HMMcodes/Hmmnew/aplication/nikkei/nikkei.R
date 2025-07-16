nikkei = quantmod::getSymbols('^N225', 
                             src = 'yahoo', 
                             #from='2002-01-01', to='2025-01-01',
                             from='2002-01-01', to='2022-04-06',
                             auto.assign = FALSE)
nikkei = na.omit(nikkei)
nikkei = data.frame(nikkei)
dates = as.Date(row.names(nikkei), '%Y-%m-%d')
nikkei = nikkei[, 'N225.Adjusted']
#View(yen)
T = length(nikkei)
log.ret = 100*(log(nikkei[2:T])-log(nikkei[1:(T-1)]))
T = length(log.ret)
ytrain=log.ret

# Plots
library(ggplot2)
df = data.frame( Return = ytrain, Tempo = dates[1:length(ytrain)] )

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

nikkei = quantmod::getSymbols('^N225', 
                              src = 'yahoo', 
                              #from='2002-01-01', to='2025-01-01',
                              from='2022-04-06', to='2025-01-01',
                              auto.assign = FALSE)
nikkei = na.omit(nikkei)
nikkei = data.frame(nikkei)
dates = as.Date(row.names(nikkei), '%Y-%m-%d')
nikkei = nikkei[, 'N225.Adjusted']
#View(yen)
T = length(nikkei)
log.ret = 100*(log(nikkei[2:T])-log(nikkei[1:(T-1)]))
T = length(log.ret)
ytest=log.ret
################################################################################
####### Fitting
################################################################################
# SVM-N
source('~/Documentos/svmHMM/aplic_n.R')
save(Results, times, h_hat, DIC, LPS, psressvmn,
     file='~/Documentos/svmHMM/nikkei/n.RData')

# SVM-t
source('~/Documentos/svmHMM/aplic_t.R')
save(Results, times, h_hat, DIC, LPS, psressvmt,
     file='~/Documentos/svmHMM/nikkei/t.RData')

# SVM-S
source('~/Documentos/svmHMM/aplic_s.R')
save(Results, times, h_hat, DIC, LPS, psressvms,
     file='~/Documentos/svmHMM/nikkei/s.RData')

# SVM-VG
source('~/Documentos/svmHMM/aplic_vg.R')
save(Results, times, h_hat, DIC, LPS, psressvmvg,
     file='~/Documentos/svmHMM/nikkei/vg.RData')
###############################################################################
### Analysis Results
################################################################################
# SVM-N
load("~/Documentos/svmHMM/nikkei/n.RData")
Results
times
DIC
LPS

plot(abs(ytrain), type='l', col='gray')
lines(exp(0.5*h_hat), lwd=2)
qqnorm(qnorm(psressvmn), main='SVM-N') 
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


################################################################################
# SVM-t
load("~/Documentos/svmHMM/nikkei/t.RData")
Results
times
DIC
LPS
plot(abs(ytrain), type='l', col='gray')
lines(exp(0.5*h_hat), lwd=2)
qqnorm(qnorm(psressvmn), main='SVM-t') 
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

################################################################################
# SVM-S
load("~/Documentos/svmHMM/nikkei/s.RData")
Results
times
DIC
LPS
plot(abs(ytrain), type='l', col='gray')
lines(exp(0.5*h_hat), lwd=2)
qqnorm(qnorm(psressvmn), main='SVM-S') 
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

################################################################################
# SVM-VG
load("~/Documentos/svmHMM/nikkei/vg.RData")
Results
times
DIC
LPS
plot(abs(ytrain), type='l', col='gray')
lines(exp(0.5*h_hat), lwd=2)
qqnorm(qnorm(psressvmn), main='SVM-VG') 
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
