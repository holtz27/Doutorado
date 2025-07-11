# Obtenha os dados do prata
prata = quantmod::getSymbols('SI=F', 
                            src = "yahoo", 
                            from = '2010-01-01', to = '2020-06-04',
                            auto.assign = FALSE)
prata = na.omit( prata )
prata = data.frame( prata )
dates = as.Date( row.names( prata ), "%Y-%m-%d" )
prata = prata[, 'SI.F.Adjusted']
#View(prata)
T = length( prata )
log.ret = 100 * ( log( prata[2:T] ) - log( prata[1:(T-1)] ) )
T = length( log.ret )
# Plots
library(ggplot2)
df = data.frame( Retorno = log.ret, Tempo = dates[-1] )

g = ggplot(df) + geom_line(aes(x = Tempo, y = Retorno))
g = g + scale_x_date(date_breaks = "48 month", date_labels = "%b %Y")
g = g + theme_test() + theme(axis.title.y = element_text(size = 18),
                             axis.text.x = element_text(size = 16),
                             axis.text.y = element_text(size = 18))
g = g + xlab('')
h = ggplot( df, aes(Retorno) )
h = h + geom_histogram(aes(y = after_stat(density)), bins = 40, color = 'white')
h = h + theme_test() + ylab('')
h = h + theme_test() + theme(axis.title.x = element_text(size = 18),
                             axis.text.x = element_text(size = 18),
                             axis.text.y = element_text(size = 18))
gridExtra::grid.arrange(g, h, nrow = 1, ncol = 2) 
data_summary = matrix(c( T, mean( log.ret ),
                         sd( log.ret ),
                         min( log.ret ),
                         max( log.ret ),
                         moments::skewness( log.ret ),
                         moments::kurtosis( log.ret ) ), nrow = 1)
colnames( data_summary ) = c( 'T','mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round( data_summary, digits = 4 )
