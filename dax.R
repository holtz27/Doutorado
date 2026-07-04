# -------------------------------------------------------------------------
# Data Collection and Preprocessing - DAX 30
# -------------------------------------------------------------------------

# 1. Download raw data from Yahoo Finance
dax <- quantmod::getSymbols('^GDAXI',
                            src = 'yahoo',
                            from = '2002-01-01', 
                            to = '2022-04-06',
                            auto.assign = FALSE)

# 2. Preprocessing
dax <- na.omit(dax)
dax <- data.frame(dax)
dates <- as.Date(row.names(dax), '%Y-%m-%d')
dax_adj <- dax[, 'GDAXI.Adjusted']

# 3. Calculate daily percentage log-returns
T_total <- length(dax_adj)
log.ret <- 100 * (log(dax_adj[2:T_total]) - log(dax_adj[1:(T_total - 1)]))
T_ret <- length(log.ret)
ytrain <- log.ret

# 4. Descriptive Statistics Summary
data_summary <- matrix(c(T_ret, 
                         mean(log.ret),
                         sd(log.ret),
                         min(log.ret),
                         max(log.ret),
                         moments::skewness(log.ret),
                         moments::kurtosis(log.ret)), nrow = 1)

colnames(data_summary) <- c('T', 'Mean', 'SD', 'Min', 'Max', 'Skewness', 'Kurtosis')
print(round(data_summary, digits = 3))

# 5. Exploratory Data Visualization
library(ggplot2)
df <- data.frame(Return = log.ret, Tempo = dates[-1])

# Time series plot of log-returns
g <- ggplot(df) + 
  geom_line(aes(x = Tempo, y = Return)) +
  scale_x_date(date_breaks = "18 month", date_labels = "%b %Y") +
  theme_test() + 
  theme(axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 18)) +
  xlab('')

# Histogram / Empirical Density plot
h <- ggplot(df, aes(Return)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 40, color = 'white') +
  theme_test() + 
  ylab('') +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 18))

# Arrange and display both plots side by side
gridExtra::grid.arrange(g, h, nrow = 1, ncol = 2)

# 6. Export data for modeling
save(dax_adj, file = 'dax.RData')

