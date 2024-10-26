# Reading data
apple = read.csv('https://raw.githubusercontent.com/holtz27/svmsmn/main/aplica%C3%A7%C3%A3o/apple/apple.csv')
apple = apple[, c('Date', 'Adj.Close')]
apple[, 2] = as.numeric( apple[, 2] ) 
apple = na.omit( apple )
dates = as.Date( apple[, 1], "%Y-%m-%d" )
#View(apple)
T = nrow(apple)
log.ret = 100 * ( log( apple[2:T, 2] ) - log( apple[1:(T-1), 2] ) )
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
data_summary = matrix(c( mean( log.ret ),
                         sd( log.ret ),
                         min( log.ret ),
                         max( log.ret ),
                         moments::skewness( log.ret ),
                         moments::kurtosis( log.ret ) ), nrow = 1)
colnames( data_summary ) = c( 'mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round( data_summary, digits = 4 )












































