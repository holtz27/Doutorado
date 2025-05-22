# Obtenha os dados do bitcoin
bitcoin = quantmod::getSymbols('BTC-USD', 
                               src = 'yahoo', 
                               from = '2015-04-07', to = '2020-05-01',
                               auto.assign = FALSE)
bitcoin = na.omit(bitcoin)
bitcoin = data.frame(bitcoin)
dates = as.Date(row.names(bitcoin), '%Y-%m-%d')
bitcoin = bitcoin[,'BTC.USD.Adjusted']
#View(bitcoin)
T = length(bitcoin)
log.ret = 100*(log( bitcoin[2:T])-log(bitcoin[1:(T-1)]))
T = length(log.ret)
# Plots
library(ggplot2)
df = data.frame(Return=log.ret, Tempo=dates[-1])
g = ggplot(df) + geom_line(aes(x=Tempo, y=Return))
g = g + scale_x_date(date_breaks="10 month", date_labels="%b %Y")
g = g + theme_test() + theme(axis.title.y=element_text(size=18),
                             axis.text.x=element_text(size=11),
                             axis.text.y=element_text(size=18))
g = g + xlab('')
h = ggplot(df, aes(Return))
h = h + geom_histogram(aes(y = after_stat(density)), bins = 40, color='white')
h = h + theme_test() + ylab('')
h = h + theme_test() + theme(axis.title.x=element_text(size=18),
                             axis.text.x=element_text(size=12),
                             axis.text.y=element_text(size=18))
gridExtra::grid.arrange(g, h, nrow=1, ncol=2) 
data_summary = matrix(c(mean(log.ret),
                        sd(log.ret),
                        min(log.ret),
                        max(log.ret),
                        moments::skewness(log.ret),
                        moments::kurtosis(log.ret)), nrow=1)
colnames(data_summary)=c('mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round(data_summary, digits=3)



par(mfrow=c(1,2))
y=log.ret-mean(log.ret)
plot(y[1:(length(y)-120)], type='l')
plot(tail(y, 120), type='l')
par(mfrow=c(1,1))

data_summary = matrix(c(mean(y[1:(length(y)-120)]),
                        sd(y[1:(length(y)-120)]),
                        min(y[1:(length(y)-120)]),
                        max(y[1:(length(y)-120)]),
                        moments::skewness(y[1:(length(y)-120)]),
                        moments::kurtosis(y[1:(length(y)-120)])), nrow=1)
round(data_summary, digits=3)
