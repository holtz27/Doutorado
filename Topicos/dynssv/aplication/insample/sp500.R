# Obtenha os dados do yen
sp500 = quantmod::getSymbols('^GSPC', 
                           src = 'yahoo', 
                           from = '2005-01-03', to = '2009-12-09',
                           auto.assign = FALSE)
sp500 = na.omit( sp500 )
sp500 = data.frame( sp500 )
dates = as.Date( row.names( sp500 ), '%Y-%m-%d' )
sp500 = sp500[, 'GSPC.Adjusted']
#View(yen)
T = length( sp500 )
log.ret = 100 * ( log( sp500[2:T] ) - log( sp500[1:(T-1)] ) )
T = length( log.ret )
# Plots
library(ggplot2)
df = data.frame( Return = log.ret, Tempo = dates[-1] )

g = ggplot(df) + geom_line(aes(x = Tempo, y = Return))
g = g + scale_x_date(date_breaks = "36 month", date_labels = "%b %Y")
g = g + theme_test() + theme(axis.title.y = element_text(size = 18),
                             axis.text.x = element_text(size = 16),
                             axis.text.y = element_text(size = 18))
g = g + xlab('')
h = ggplot( df, aes(Return) )
h = h + geom_histogram(aes(y = after_stat(density)), bins = 40, color = 'white')
h = h + theme_test() + ylab('')
h = h + theme_test() + theme(axis.title.x = element_text(size = 18),
                             axis.text.x = element_text(size = 18),
                             axis.text.y = element_text(size = 18))
gridExtra::grid.arrange(g, h, nrow = 1, ncol = 2) 
data_summary = matrix(c( mean( log.ret ),
                         sd( log.ret ),
                         min( log.ret ),
                         max( log.ret ),
                         moments::skewness( log.ret ),
                         moments::kurtosis( log.ret ) ), nrow = 1)
colnames( data_summary ) = c( 'mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round( data_summary, digits = 3 )


