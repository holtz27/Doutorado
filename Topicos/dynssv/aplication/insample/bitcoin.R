# Obtenha os dados do bitcoin
bitcoin = quantmod::getSymbols('BTC-USD', 
                               src = 'yahoo', 
                               from = '2011-09-13', to = '2017-12-31',
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

g = ggplot(df) + geom_line(aes(x = Tempo, y = Return))
g = g + scale_x_date(date_breaks = "36 month", date_labels = "%b %Y")
g = g + theme_test() + theme(axis.title.y = element_text(size = 18),
                             axis.text.x = element_text(size = 16),
                             axis.text.y = element_text(size = 18))
g = g + xlab('')
h = ggplot(df, aes(Return))
h = h + geom_histogram(aes(y = after_stat(density)), bins = 40, color = 'white')
h = h + theme_test() + ylab('')
h = h + theme_test() + theme(axis.title.x = element_text(size = 18),
                             axis.text.x = element_text(size = 18),
                             axis.text.y = element_text(size = 18))
gridExtra::grid.arrange(g, h, nrow = 1, ncol = 2) 
data_summary = matrix(c( mean(log.ret),
                         sd(log.ret),
                         min(log.ret),
                         max(log.ret),
                         moments::skewness(log.ret),
                         moments::kurtosis(log.ret)), nrow = 1)
colnames(data_summary)=c('mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round(data_summary, digits=3)


