ibovespa = quantmod::getSymbols('^BVSP', 
                             src = 'yahoo', 
                             from = '1998-03-03', to = '2018-1-01',
                             auto.assign = FALSE)
ibovespa = na.omit(ibovespa)
ibovespa = data.frame(ibovespa)
dates = as.Date(row.names(ibovespa), '%Y-%m-%d')
ibovespa = ibovespa[, 'BVSP.Adjusted']
#View(yen)
T=length(ibovespa)
log.ret=100*(log(ibovespa[2:T])-log(ibovespa[1:(T-1)]))
T=length(log.ret)
# Plots
library(ggplot2)
df = data.frame(Return=log.ret, Tempo=dates[-1])

g = ggplot(df) + geom_line(aes(x = Tempo, y=Return))
g = g + scale_x_date(date_breaks = "28 month", date_labels ="%b %Y")
g = g + theme_test() + theme(axis.title.y = element_text(size=18),
                             axis.text.x = element_text(size=12),
                             axis.text.y = element_text(size=18))
g = g + xlab('')
h = ggplot( df, aes(Return) )
h = h + geom_histogram(aes(y = after_stat(density)), bins = 40, color = 'white')
h = h + theme_test() + ylab('')
h = h + theme_test() + theme(axis.title.x=element_text(size=18),
                             axis.text.x=element_text(size=18),
                             axis.text.y=element_text(size=18))
gridExtra::grid.arrange(g, h, nrow=1, ncol=2) 
data_summary = matrix(c( T, mean(log.ret),
                         sd(log.ret),
                         min(log.ret),
                         max(log.ret),
                         moments::skewness(log.ret),
                         moments::kurtosis(log.ret)), nrow=1)
colnames(data_summary) = c('T','mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round(data_summary, digits=3)

################################################################################
####### Fitting
################################################################################
# SVM-N
source('~/Documentos/hmm/aplic_n.R')
save(Results, times, h_hat, DIC, LPS, file='~/Documentos/hmm/ibovespa/n.RData')

# SVM-t
source('~/Documentos/hmm/aplic_t.R')
save(Results, times, h_hat, DIC, LPS, file='~/Documentos/hmm/ibovespa/t.RData')

# SVM-S
source('~/Documentos/hmm/aplic_s.R')
save(Results, times, h_hat, DIC, LPS, file='~/Documentos/hmm/ibovespa/s.RData')

# SVM-VG
source('~/Documentos/hmm/aplic_vg.R')
save(Results, times, h_hat, DIC, LPS, file='~/Documentos/hmm/ibovespa/vg.RData')

