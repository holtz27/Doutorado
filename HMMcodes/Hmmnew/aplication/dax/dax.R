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
load("~/Documentos/svmHMM/dax/n.RData")
Results
times
DIC
LPS
plot(exp(0.5*h_hat), type='l')

load("~/Documentos/svmHMM/dax/t.RData")
Results
times
DIC
LPS
plot(exp(0.5*h_hat), type='l')

load("~/Documentos/svmHMM/dax/s.RData")
Results
times
DIC
LPS
plot(exp(0.5*h_hat), type='l')

load("~/Documentos/svmHMM/dax/vg.RData")
Results
times
DIC
LPS
plot(exp(0.5*h_hat), type='l')

