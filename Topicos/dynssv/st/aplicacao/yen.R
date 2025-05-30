# Obtenha os dados do yen
yen = quantmod::getSymbols('JPY=X', 
                           src = 'yahoo', 
                           from = '2009-01-01', to = '2019-10-31',
                           auto.assign = FALSE)
yen = na.omit( yen )
yen = data.frame( yen )
dates = as.Date( row.names( yen ), '%Y-%m-%d' )
yen = yen[, 'JPY.X.Adjusted']
#View(yen)
T = length( yen )
log.ret = 100 * ( log( yen[2:T] ) - log( yen[1:(T-1)] ) )
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


