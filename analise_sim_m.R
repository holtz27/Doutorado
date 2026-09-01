library(ggplot2)

################################### CENARIO 1 ##################################
### b=2.5
load("~/tese/HMM/Res_b_2.RData")

M=Time=b0=b1=b2=mu=phi=sigma=nu=matrix(nrow=100,ncol=10)
for(i in 1:100){
  for(j in 1:10){
    status=Results[[i]][[j]]$status
    if(status=="success"){
      M[i,j]=Results[[i]][[j]]$logL
      Time[i,j]=Results[[i]][[j]]$time
      b0[i,j]=Results[[i]][[j]]$summary[1,1]
      b1[i,j]=Results[[i]][[j]]$summary[2,1]
      b2[i,j]=Results[[i]][[j]]$summary[3,1]
      mu[i,j]=Results[[i]][[j]]$summary[4,1]
      phi[i,j]=Results[[i]][[j]]$summary[5,1]
      sigma[i,j]=Results[[i]][[j]]$summary[6,1]
      nu[i,j]=Results[[i]][[j]]$summary[7,1]
    }
  }
}

M=na.omit(M)
Time=na.omit(Time)
b0=na.omit(b0)
b1=na.omit(b1)
b2=na.omit(b2)
mu=na.omit(mu)
phi=na.omit(phi)
sigma=na.omit(sigma)
nu=na.omit(nu)

# Estimates analisys
m=seq(25, 250, by=25)
par(mfrow=c(1,2))
boxplot(phi[,1],phi[,2],phi[,3],phi[,4],phi[,5],
        phi[,6],phi[,7],phi[,8],phi[,9],phi[,10], main='phi', names=m)
abline(h=0.98, col='red', lty=2, lwd=2)
boxplot(sigma[,1],sigma[,2],sigma[,3],sigma[,4],sigma[,5],
        sigma[,6],sigma[,7],sigma[,8],sigma[,9],sigma[,10], main='sigma', names=m)
abline(h=0.1, col='red', lty=2, lwd=2)
par(mfrow=c(1,1))

S25=S50=S75=S100=S125=S150=S175=S200=S225=S250=matrix(0, nrow=7 ,ncol=5)
n=nrow(M)
for(k in 1:n){
  S25  = S25  + Results[[k]][[1]]$summary
  S50  = S50  + Results[[k]][[2]]$summary
  S75  = S75  + Results[[k]][[3]]$summary
  S100 = S100 + Results[[k]][[4]]$summary
  S125 = S125 + Results[[k]][[5]]$summary
  S150 = S150 + Results[[k]][[6]]$summary
  S175 = S175 + Results[[k]][[7]]$summary
  S200 = S200 + Results[[k]][[8]]$summary
  S225 = S225 + Results[[k]][[9]]$summary
  S250 = S250 + Results[[k]][[10]]$summary
}
S25=S25[,c(1,3,5)]/100
S50=S50[,c(1,3,5)]/100
S75=S75[,c(1,3,5)]/100
S100=S100[,c(1,3,5)]/100
S125=S125[,c(1,3,5)]/100
S150=S150[,c(1,3,5)]/100
S175=S175[,c(1,3,5)]/100
S200=S200[,c(1,3,5)]/100
S225=S225[,c(1,3,5)]/100
S250=S250[,c(1,3,5)]/100


# 1. Extração segura dos dados (base R)
list_dfs <- lapply(seq_along(Results), function(i) {
  lapply(seq_along(Results[[i]]), function(j) {
    res <- Results[[i]][[j]]
    if (!is.null(res) && res$status == "success") {
      summary_mat <- res$summary
      data.frame(
        sim = i,
        m = res$m,
        param = rownames(summary_mat),
        mean_val = summary_mat[, "mean"],
        q2.5 = summary_mat[, "2.5%"],
        q97.5 = summary_mat[, "97.5%"]
      )
    }
  })
})

# Empilha todas as listas em um único Data Frame
df_raw <- do.call(rbind, unlist(list_dfs, recursive = FALSE))
# 2. Agregação
df_summary <- aggregate(
  cbind(mean_val, q2.5, q97.5) ~ m + param, 
  data = df_raw, 
  FUN = mean, 
  na.rm = TRUE
)
# Renomeia colunas
colnames(df_summary) <- c("m", "param", "mean_est", "lower_ic", "upper_ic")
ordem_desejada <- c("b0", "b1", "b2", "mu", "phi", "sigma", "nu")
df_summary$param <- factor(df_summary$param, levels = ordem_desejada)
# 3. Data Frame com os valores verdadeiros reordenados
true_values <- data.frame(
  param = factor(ordem_desejada, levels = ordem_desejada),
  true_val = c(0.20, 0.07, -0.18, 0.10, 0.98, 0.10, 10)
)
# 4. Gráfico com ggplot2 (ncol = 3 garante 3 painéis por linha)
ggplot(df_summary, aes(x = factor(m), y = mean_est)) +
  geom_pointrange(aes(ymin = lower_ic, ymax = upper_ic), color = "#1f77b4", size = 0.5) +
  geom_hline(data = true_values, aes(yintercept = true_val), color = "red", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ param, scales = "free_y", ncol = 3) +
  labs(
    x = "Tamanho de m",
    y = "Estimativa (Média ± IC 95%)",
    title = "Convergência dos Parâmetros por tamanho de m"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )



# Eficient analisys
Er=abs(M-M[,ncol(M)])
Er=t(apply(Er, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x - min_val) / (max_val - min_val)
}))
C=Time
C=t(apply(C, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x-min_val)/(max_val-min_val)
}))

A=1-Er
Ef=1-C
# 1. Criação do data.frame com as 3 métricas
df <- data.frame(
  m = rep(m, 3),
  Valor = c(colMeans(Ef), colMeans(A), colMeans((2 * Ef * A) / (Ef + A))),
  Metrica = factor(
    rep(c("Eficiência", "Acurácia", "F-Score"), each = length(m)),
    levels = c("Eficiência", "Acurácia", "F-Score")
  )
)

# 2. Gráfico sem título e sem quadro na legenda
ggplot(df, aes(x = m, y = Valor, color = Metrica, shape = Metrica, linetype = Metrica)) +
  #geom_vline(xintercept = 100, linetype = "dashed", color = "gray60", linewidth = 0.7) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  
  scale_x_continuous(breaks = m) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  
  scale_color_manual(values = c("Eficiência" = "#0072B2", 
                                "Acurácia"   = "#D55E00", 
                                "F-Score"    = "#009E73")) +
  scale_shape_manual(values = c("Eficiência" = 19, "Acurácia" = 19, "F-Score" = 17)) +
  scale_linetype_manual(values = c("Eficiência" = "solid", "Acurácia" = "solid", "F-Score" = "dashed")) +
  
  labs(x = "m", y = NULL, color = NULL, shape = NULL, linetype = NULL) +
  
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14, color = "black"),
    
    legend.text  = element_text(size = 12),
    legend.position = "right",
    legend.background = element_blank(),     
    legend.box.background = element_blank(), 
    panel.grid.minor = element_blank()
  )

################################### CENARIO 1 ##################################
### b=4
load("~/tese/HMM/Res_b_4.RData")
M=Time=b0=b1=b2=mu=phi=sigma=nu=matrix(nrow=100,ncol=10)
for(i in 1:100){
  for(j in 1:10){
    status=Results[[i]][[j]]$status
    if(status=="success"){
      M[i,j]=Results[[i]][[j]]$logL
      Time[i,j]=Results[[i]][[j]]$time
      b0[i,j]=Results[[i]][[j]]$summary[1,1]
      b1[i,j]=Results[[i]][[j]]$summary[2,1]
      b2[i,j]=Results[[i]][[j]]$summary[3,1]
      mu[i,j]=Results[[i]][[j]]$summary[4,1]
      phi[i,j]=Results[[i]][[j]]$summary[5,1]
      sigma[i,j]=Results[[i]][[j]]$summary[6,1]
      nu[i,j]=Results[[i]][[j]]$summary[7,1]
    }
  }
}

M=na.omit(M)
Time=na.omit(Time)
b0=na.omit(b0)
b1=na.omit(b1)
b2=na.omit(b2)
mu=na.omit(mu)
phi=na.omit(phi)
sigma=na.omit(sigma)
nu=na.omit(nu)

# Estimates analisys
m=seq(25, 250, by=25)
par(mfrow=c(1,2))
boxplot(phi[,1],phi[,2],phi[,3],phi[,4],phi[,5],
        phi[,6],phi[,7],phi[,8],phi[,9],phi[,10], main='phi', names=m)
abline(h=0.98, col='red', lty=2, lwd=2)
boxplot(sigma[,1],sigma[,2],sigma[,3],sigma[,4],sigma[,5],
        sigma[,6],sigma[,7],sigma[,8],sigma[,9],sigma[,10], main='sigma', names=m)
abline(h=0.1, col='red', lty=2, lwd=2)
par(mfrow=c(1,1))

S25=S50=S75=S100=S125=S150=S175=S200=S225=S250=matrix(0, nrow=7 ,ncol=5)
n=nrow(M)
for(k in 1:n){
  S25  = S25  + Results[[k]][[1]]$summary
  S50  = S50  + Results[[k]][[2]]$summary
  S75  = S75  + Results[[k]][[3]]$summary
  S100 = S100 + Results[[k]][[4]]$summary
  S125 = S125 + Results[[k]][[5]]$summary
  S150 = S150 + Results[[k]][[6]]$summary
  S175 = S175 + Results[[k]][[7]]$summary
  S200 = S200 + Results[[k]][[8]]$summary
  S225 = S225 + Results[[k]][[9]]$summary
  S250 = S250 + Results[[k]][[10]]$summary
}
S25=S25[,c(1,3,5)]/100
S50=S50[,c(1,3,5)]/100
S75=S75[,c(1,3,5)]/100
S100=S100[,c(1,3,5)]/100
S125=S125[,c(1,3,5)]/100
S150=S150[,c(1,3,5)]/100
S175=S175[,c(1,3,5)]/100
S200=S200[,c(1,3,5)]/100
S225=S225[,c(1,3,5)]/100
S250=S250[,c(1,3,5)]/100


# 1. Extração segura dos dados (base R)
list_dfs <- lapply(seq_along(Results), function(i) {
  lapply(seq_along(Results[[i]]), function(j) {
    res <- Results[[i]][[j]]
    if (!is.null(res) && res$status == "success") {
      summary_mat <- res$summary
      data.frame(
        sim = i,
        m = res$m,
        param = rownames(summary_mat),
        mean_val = summary_mat[, "mean"],
        q2.5 = summary_mat[, "2.5%"],
        q97.5 = summary_mat[, "97.5%"]
      )
    }
  })
})

# Empilha todas as listas em um único Data Frame
df_raw <- do.call(rbind, unlist(list_dfs, recursive = FALSE))
# 2. Agregação
df_summary <- aggregate(
  cbind(mean_val, q2.5, q97.5) ~ m + param, 
  data = df_raw, 
  FUN = mean, 
  na.rm = TRUE
)
# Renomeia colunas
colnames(df_summary) <- c("m", "param", "mean_est", "lower_ic", "upper_ic")
ordem_desejada <- c("b0", "b1", "b2", "mu", "phi", "sigma", "nu")
df_summary$param <- factor(df_summary$param, levels = ordem_desejada)
# 3. Data Frame com os valores verdadeiros reordenados
true_values <- data.frame(
  param = factor(ordem_desejada, levels = ordem_desejada),
  true_val = c(0.20, 0.07, -0.18, 0.10, 0.98, 0.10, 10)
)
# 4. Gráfico com ggplot2 (ncol = 3 garante 3 painéis por linha)
ggplot(df_summary, aes(x = factor(m), y = mean_est)) +
  geom_pointrange(aes(ymin = lower_ic, ymax = upper_ic), color = "#1f77b4", size = 0.5) +
  geom_hline(data = true_values, aes(yintercept = true_val), color = "red", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ param, scales = "free_y", ncol = 3) +
  labs(
    x = "Tamanho de m",
    y = "Estimativa (Média ± IC 95%)",
    title = "Convergência dos Parâmetros por tamanho de m"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )



# Eficient analisys
Er=abs(M-M[,ncol(M)])
Er=t(apply(Er, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x - min_val) / (max_val - min_val)
}))
C=Time
C=t(apply(C, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x-min_val)/(max_val-min_val)
}))

A=1-Er
Ef=1-C
# 1. Criação do data.frame com as 3 métricas
df <- data.frame(
  m = rep(m, 3),
  Valor = c(colMeans(Ef), colMeans(A), colMeans((2 * Ef * A) / (Ef + A))),
  Metrica = factor(
    rep(c("Eficiência", "Acurácia", "F-Score"), each = length(m)),
    levels = c("Eficiência", "Acurácia", "F-Score")
  )
)

# 2. Gráfico sem título e sem quadro na legenda
ggplot(df, aes(x = m, y = Valor, color = Metrica, shape = Metrica, linetype = Metrica)) +
  #geom_vline(xintercept = 100, linetype = "dashed", color = "gray60", linewidth = 0.7) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  
  scale_x_continuous(breaks = m) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  
  scale_color_manual(values = c("Eficiência" = "#0072B2", 
                                "Acurácia"   = "#D55E00", 
                                "F-Score"    = "#009E73")) +
  scale_shape_manual(values = c("Eficiência" = 19, "Acurácia" = 19, "F-Score" = 17)) +
  scale_linetype_manual(values = c("Eficiência" = "solid", "Acurácia" = "solid", "F-Score" = "dashed")) +
  
  labs(x = "m", y = NULL, color = NULL, shape = NULL, linetype = NULL) +
  
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14, color = "black"),
    
    legend.text  = element_text(size = 12),
    legend.position = "right",
    legend.background = element_blank(),     
    legend.box.background = element_blank(), 
    panel.grid.minor = element_blank()
  )


################################### CENARIO 1 ##################################
### b=5
load("~/tese/HMM/Res_b_5.RData")
M=Time=b0=b1=b2=mu=phi=sigma=nu=matrix(nrow=100,ncol=10)
for(i in 1:100){
  for(j in 1:10){
    status=Results[[i]][[j]]$status
    if(status=="success"){
      M[i,j]=Results[[i]][[j]]$logL
      Time[i,j]=Results[[i]][[j]]$time
      b0[i,j]=Results[[i]][[j]]$summary[1,1]
      b1[i,j]=Results[[i]][[j]]$summary[2,1]
      b2[i,j]=Results[[i]][[j]]$summary[3,1]
      mu[i,j]=Results[[i]][[j]]$summary[4,1]
      phi[i,j]=Results[[i]][[j]]$summary[5,1]
      sigma[i,j]=Results[[i]][[j]]$summary[6,1]
      nu[i,j]=Results[[i]][[j]]$summary[7,1]
    }
  }
}

M=na.omit(M)
Time=na.omit(Time)
b0=na.omit(b0)
b1=na.omit(b1)
b2=na.omit(b2)
mu=na.omit(mu)
phi=na.omit(phi)
sigma=na.omit(sigma)
nu=na.omit(nu)

# Estimates analisys
m=seq(25, 250, by=25)
par(mfrow=c(1,2))
boxplot(phi[,1],phi[,2],phi[,3],phi[,4],phi[,5],
        phi[,6],phi[,7],phi[,8],phi[,9],phi[,10], main='phi', names=m)
abline(h=0.98, col='red', lty=2, lwd=2)
boxplot(sigma[,1],sigma[,2],sigma[,3],sigma[,4],sigma[,5],
        sigma[,6],sigma[,7],sigma[,8],sigma[,9],sigma[,10], main='sigma', names=m)
abline(h=0.1, col='red', lty=2, lwd=2)
par(mfrow=c(1,1))

S25=S50=S75=S100=S125=S150=S175=S200=S225=S250=matrix(0, nrow=7 ,ncol=5)
n=nrow(M)
for(k in 1:n){
  S25  = S25  + Results[[k]][[1]]$summary
  S50  = S50  + Results[[k]][[2]]$summary
  S75  = S75  + Results[[k]][[3]]$summary
  S100 = S100 + Results[[k]][[4]]$summary
  S125 = S125 + Results[[k]][[5]]$summary
  S150 = S150 + Results[[k]][[6]]$summary
  S175 = S175 + Results[[k]][[7]]$summary
  S200 = S200 + Results[[k]][[8]]$summary
  S225 = S225 + Results[[k]][[9]]$summary
  S250 = S250 + Results[[k]][[10]]$summary
}
S25=S25[,c(1,3,5)]/100
S50=S50[,c(1,3,5)]/100
S75=S75[,c(1,3,5)]/100
S100=S100[,c(1,3,5)]/100
S125=S125[,c(1,3,5)]/100
S150=S150[,c(1,3,5)]/100
S175=S175[,c(1,3,5)]/100
S200=S200[,c(1,3,5)]/100
S225=S225[,c(1,3,5)]/100
S250=S250[,c(1,3,5)]/100


# 1. Extração segura dos dados (base R)
list_dfs <- lapply(seq_along(Results), function(i) {
  lapply(seq_along(Results[[i]]), function(j) {
    res <- Results[[i]][[j]]
    if (!is.null(res) && res$status == "success") {
      summary_mat <- res$summary
      data.frame(
        sim = i,
        m = res$m,
        param = rownames(summary_mat),
        mean_val = summary_mat[, "mean"],
        q2.5 = summary_mat[, "2.5%"],
        q97.5 = summary_mat[, "97.5%"]
      )
    }
  })
})

# Empilha todas as listas em um único Data Frame
df_raw <- do.call(rbind, unlist(list_dfs, recursive = FALSE))
# 2. Agregação
df_summary <- aggregate(
  cbind(mean_val, q2.5, q97.5) ~ m + param, 
  data = df_raw, 
  FUN = mean, 
  na.rm = TRUE
)
# Renomeia colunas
colnames(df_summary) <- c("m", "param", "mean_est", "lower_ic", "upper_ic")
ordem_desejada <- c("b0", "b1", "b2", "mu", "phi", "sigma", "nu")
df_summary$param <- factor(df_summary$param, levels = ordem_desejada)
# 3. Data Frame com os valores verdadeiros reordenados
true_values <- data.frame(
  param = factor(ordem_desejada, levels = ordem_desejada),
  true_val = c(0.20, 0.07, -0.18, 0.10, 0.98, 0.10, 10)
)
# 4. Gráfico com ggplot2 (ncol = 3 garante 3 painéis por linha)
ggplot(df_summary, aes(x = factor(m), y = mean_est)) +
  geom_pointrange(aes(ymin = lower_ic, ymax = upper_ic), color = "#1f77b4", size = 0.5) +
  geom_hline(data = true_values, aes(yintercept = true_val), color = "red", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ param, scales = "free_y", ncol = 3) +
  labs(
    x = "Tamanho de m",
    y = "Estimativa (Média ± IC 95%)",
    title = "Convergência dos Parâmetros por tamanho de m"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )



# Eficient analisys
Er=abs(M-M[,ncol(M)])
Er=t(apply(Er, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x - min_val) / (max_val - min_val)
}))
C=Time
C=t(apply(C, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x-min_val)/(max_val-min_val)
}))

A=1-Er
Ef=1-C
# 1. Criação do data.frame com as 3 métricas
df <- data.frame(
  m = rep(m, 3),
  Valor = c(colMeans(Ef), colMeans(A), colMeans((2 * Ef * A) / (Ef + A))),
  Metrica = factor(
    rep(c("Eficiência", "Acurácia", "F-Score"), each = length(m)),
    levels = c("Eficiência", "Acurácia", "F-Score")
  )
)

# 2. Gráfico sem título e sem quadro na legenda
ggplot(df, aes(x = m, y = Valor, color = Metrica, shape = Metrica, linetype = Metrica)) +
  #geom_vline(xintercept = 100, linetype = "dashed", color = "gray60", linewidth = 0.7) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  
  scale_x_continuous(breaks = m) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  
  scale_color_manual(values = c("Eficiência" = "#0072B2", 
                                "Acurácia"   = "#D55E00", 
                                "F-Score"    = "#009E73")) +
  scale_shape_manual(values = c("Eficiência" = 19, "Acurácia" = 19, "F-Score" = 17)) +
  scale_linetype_manual(values = c("Eficiência" = "solid", "Acurácia" = "solid", "F-Score" = "dashed")) +
  
  labs(x = "m", y = NULL, color = NULL, shape = NULL, linetype = NULL) +
  
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14, color = "black"),
    
    legend.text  = element_text(size = 12),
    legend.position = "right",
    legend.background = element_blank(),     
    legend.box.background = element_blank(), 
    panel.grid.minor = element_blank()
  )


################################### CENARIO 2 ##################################
### b=2.5
load("~/tese/HMM/2Res_2.5.RData")
M=Time=b0=b1=b2=mu=phi=sigma=nu=matrix(nrow=100,ncol=10)
for(i in 1:100){
  for(j in 1:10){
    status=Results[[i]][[j]]$status
    if(status=="success"){
      M[i,j]=Results[[i]][[j]]$logL
      Time[i,j]=Results[[i]][[j]]$time
      b0[i,j]=Results[[i]][[j]]$summary[1,1]
      b1[i,j]=Results[[i]][[j]]$summary[2,1]
      b2[i,j]=Results[[i]][[j]]$summary[3,1]
      mu[i,j]=Results[[i]][[j]]$summary[4,1]
      phi[i,j]=Results[[i]][[j]]$summary[5,1]
      sigma[i,j]=Results[[i]][[j]]$summary[6,1]
      nu[i,j]=Results[[i]][[j]]$summary[7,1]
    }
  }
}

M=na.omit(M)
Time=na.omit(Time)
b0=na.omit(b0)
b1=na.omit(b1)
b2=na.omit(b2)
mu=na.omit(mu)
phi=na.omit(phi)
sigma=na.omit(sigma)
nu=na.omit(nu)

# Estimates analisys
m=seq(25, 250, by=25)
par(mfrow=c(1,2))
boxplot(phi[,1],phi[,2],phi[,3],phi[,4],phi[,5],
        phi[,6],phi[,7],phi[,8],phi[,9],phi[,10], main='phi', names=m)
abline(h=0.99, col='red', lty=2, lwd=2)
boxplot(sigma[,1],sigma[,2],sigma[,3],sigma[,4],sigma[,5],
        sigma[,6],sigma[,7],sigma[,8],sigma[,9],sigma[,10], main='sigma', names=m)
abline(h=0.15, col='red', lty=2, lwd=2)
par(mfrow=c(1,1))

S25=S50=S75=S100=S125=S150=S175=S200=S225=S250=matrix(0, nrow=7 ,ncol=5)
n=nrow(M)
for(k in 1:n){
  S25  = S25  + Results[[k]][[1]]$summary
  S50  = S50  + Results[[k]][[2]]$summary
  S75  = S75  + Results[[k]][[3]]$summary
  S100 = S100 + Results[[k]][[4]]$summary
  S125 = S125 + Results[[k]][[5]]$summary
  S150 = S150 + Results[[k]][[6]]$summary
  S175 = S175 + Results[[k]][[7]]$summary
  S200 = S200 + Results[[k]][[8]]$summary
  S225 = S225 + Results[[k]][[9]]$summary
  S250 = S250 + Results[[k]][[10]]$summary
}
S25=S25[,c(1,3,5)]/100
S50=S50[,c(1,3,5)]/100
S75=S75[,c(1,3,5)]/100
S100=S100[,c(1,3,5)]/100
S125=S125[,c(1,3,5)]/100
S150=S150[,c(1,3,5)]/100
S175=S175[,c(1,3,5)]/100
S200=S200[,c(1,3,5)]/100
S225=S225[,c(1,3,5)]/100
S250=S250[,c(1,3,5)]/100


# 1. Extração segura dos dados (base R)
list_dfs <- lapply(seq_along(Results), function(i) {
  lapply(seq_along(Results[[i]]), function(j) {
    res <- Results[[i]][[j]]
    if (!is.null(res) && res$status == "success") {
      summary_mat <- res$summary
      data.frame(
        sim = i,
        m = res$m,
        param = rownames(summary_mat),
        mean_val = summary_mat[, "mean"],
        q2.5 = summary_mat[, "2.5%"],
        q97.5 = summary_mat[, "97.5%"]
      )
    }
  })
})

# Empilha todas as listas em um único Data Frame
df_raw <- do.call(rbind, unlist(list_dfs, recursive = FALSE))
# 2. Agregação
df_summary <- aggregate(
  cbind(mean_val, q2.5, q97.5) ~ m + param, 
  data = df_raw, 
  FUN = mean, 
  na.rm = TRUE
)
# Renomeia colunas
colnames(df_summary) <- c("m", "param", "mean_est", "lower_ic", "upper_ic")
ordem_desejada <- c("b0", "b1", "b2", "mu", "phi", "sigma", "nu")
df_summary$param <- factor(df_summary$param, levels = ordem_desejada)
# 3. Data Frame com os valores verdadeiros reordenados
true_values <- data.frame(
  param = factor(ordem_desejada, levels = ordem_desejada),
  true_val = c(0.10, -0.01, -0.05, 0.10, 0.99, 0.15, 5)
)
# 4. Gráfico com ggplot2 (ncol = 3 garante 3 painéis por linha)
ggplot(df_summary, aes(x = factor(m), y = mean_est)) +
  geom_pointrange(aes(ymin = lower_ic, ymax = upper_ic), color = "#1f77b4", size = 0.5) +
  geom_hline(data = true_values, aes(yintercept = true_val), color = "red", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ param, scales = "free_y", ncol = 3) +
  labs(
    x = "Tamanho de m",
    y = "Estimativa (Média ± IC 95%)",
    title = "Convergência dos Parâmetros por tamanho de m"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )



# Eficient analisys
Er=abs(M-M[,ncol(M)])
Er=t(apply(Er, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x - min_val) / (max_val - min_val)
}))
C=Time
C=t(apply(C, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x-min_val)/(max_val-min_val)
}))

A=1-Er
Ef=1-C
# 1. Criação do data.frame com as 3 métricas
df <- data.frame(
  m = rep(m, 3),
  Valor = c(colMeans(Ef), colMeans(A), colMeans((2 * Ef * A) / (Ef + A))),
  Metrica = factor(
    rep(c("Eficiência", "Acurácia", "F-Score"), each = length(m)),
    levels = c("Eficiência", "Acurácia", "F-Score")
  )
)

# 2. Gráfico sem título e sem quadro na legenda
ggplot(df, aes(x = m, y = Valor, color = Metrica, shape = Metrica, linetype = Metrica)) +
  #geom_vline(xintercept = 100, linetype = "dashed", color = "gray60", linewidth = 0.7) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  
  scale_x_continuous(breaks = m) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  
  scale_color_manual(values = c("Eficiência" = "#0072B2", 
                                "Acurácia"   = "#D55E00", 
                                "F-Score"    = "#009E73")) +
  scale_shape_manual(values = c("Eficiência" = 19, "Acurácia" = 19, "F-Score" = 17)) +
  scale_linetype_manual(values = c("Eficiência" = "solid", "Acurácia" = "solid", "F-Score" = "dashed")) +
  
  labs(x = "m", y = NULL, color = NULL, shape = NULL, linetype = NULL) +
  
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14, color = "black"),
    
    legend.text  = element_text(size = 12),
    legend.position = "right",
    legend.background = element_blank(),     
    legend.box.background = element_blank(), 
    panel.grid.minor = element_blank()
  )
################################### CENARIO 2 ##################################
### b=4
load("~/tese/HMM/2Res_4.RData")
M=Time=b0=b1=b2=mu=phi=sigma=nu=matrix(nrow=100,ncol=10)
for(i in 1:100){
  for(j in 1:10){
    status=Results[[i]][[j]]$status
    if(status=="success"){
      M[i,j]=Results[[i]][[j]]$logL
      Time[i,j]=Results[[i]][[j]]$time
      b0[i,j]=Results[[i]][[j]]$summary[1,1]
      b1[i,j]=Results[[i]][[j]]$summary[2,1]
      b2[i,j]=Results[[i]][[j]]$summary[3,1]
      mu[i,j]=Results[[i]][[j]]$summary[4,1]
      phi[i,j]=Results[[i]][[j]]$summary[5,1]
      sigma[i,j]=Results[[i]][[j]]$summary[6,1]
      nu[i,j]=Results[[i]][[j]]$summary[7,1]
    }
  }
}

M=na.omit(M)
Time=na.omit(Time)
b0=na.omit(b0)
b1=na.omit(b1)
b2=na.omit(b2)
mu=na.omit(mu)
phi=na.omit(phi)
sigma=na.omit(sigma)
nu=na.omit(nu)

# Estimates analisys
m=seq(25, 250, by=25)
par(mfrow=c(1,2))
boxplot(phi[,1],phi[,2],phi[,3],phi[,4],phi[,5],
        phi[,6],phi[,7],phi[,8],phi[,9],phi[,10], main='phi', names=m)
abline(h=0.99, col='red', lty=2, lwd=2)
boxplot(sigma[,1],sigma[,2],sigma[,3],sigma[,4],sigma[,5],
        sigma[,6],sigma[,7],sigma[,8],sigma[,9],sigma[,10], main='sigma', names=m)
abline(h=0.15, col='red', lty=2, lwd=2)
par(mfrow=c(1,1))

S25=S50=S75=S100=S125=S150=S175=S200=S225=S250=matrix(0, nrow=7 ,ncol=5)
n=nrow(M)
for(k in 1:n){
  S25  = S25  + Results[[k]][[1]]$summary
  S50  = S50  + Results[[k]][[2]]$summary
  S75  = S75  + Results[[k]][[3]]$summary
  S100 = S100 + Results[[k]][[4]]$summary
  S125 = S125 + Results[[k]][[5]]$summary
  S150 = S150 + Results[[k]][[6]]$summary
  S175 = S175 + Results[[k]][[7]]$summary
  S200 = S200 + Results[[k]][[8]]$summary
  S225 = S225 + Results[[k]][[9]]$summary
  S250 = S250 + Results[[k]][[10]]$summary
}
S25=S25[,c(1,3,5)]/100
S50=S50[,c(1,3,5)]/100
S75=S75[,c(1,3,5)]/100
S100=S100[,c(1,3,5)]/100
S125=S125[,c(1,3,5)]/100
S150=S150[,c(1,3,5)]/100
S175=S175[,c(1,3,5)]/100
S200=S200[,c(1,3,5)]/100
S225=S225[,c(1,3,5)]/100
S250=S250[,c(1,3,5)]/100


# 1. Extração segura dos dados (base R)
list_dfs <- lapply(seq_along(Results), function(i) {
  lapply(seq_along(Results[[i]]), function(j) {
    res <- Results[[i]][[j]]
    if (!is.null(res) && res$status == "success") {
      summary_mat <- res$summary
      data.frame(
        sim = i,
        m = res$m,
        param = rownames(summary_mat),
        mean_val = summary_mat[, "mean"],
        q2.5 = summary_mat[, "2.5%"],
        q97.5 = summary_mat[, "97.5%"]
      )
    }
  })
})

# Empilha todas as listas em um único Data Frame
df_raw <- do.call(rbind, unlist(list_dfs, recursive = FALSE))
# 2. Agregação
df_summary <- aggregate(
  cbind(mean_val, q2.5, q97.5) ~ m + param, 
  data = df_raw, 
  FUN = mean, 
  na.rm = TRUE
)
# Renomeia colunas
colnames(df_summary) <- c("m", "param", "mean_est", "lower_ic", "upper_ic")
ordem_desejada <- c("b0", "b1", "b2", "mu", "phi", "sigma", "nu")
df_summary$param <- factor(df_summary$param, levels = ordem_desejada)
# 3. Data Frame com os valores verdadeiros reordenados
true_values <- data.frame(
  param = factor(ordem_desejada, levels = ordem_desejada),
  true_val = c(0.10, -0.01, -0.05, 0.10, 0.99, 0.15, 5)
)
# 4. Gráfico com ggplot2 (ncol = 3 garante 3 painéis por linha)
ggplot(df_summary, aes(x = factor(m), y = mean_est)) +
  geom_pointrange(aes(ymin = lower_ic, ymax = upper_ic), color = "#1f77b4", size = 0.5) +
  geom_hline(data = true_values, aes(yintercept = true_val), color = "red", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ param, scales = "free_y", ncol = 3) +
  labs(
    x = "Tamanho de m",
    y = "Estimativa (Média ± IC 95%)",
    title = "Convergência dos Parâmetros por tamanho de m"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )



# Eficient analisys
Er=abs(M-M[,ncol(M)])
Er=t(apply(Er, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x - min_val) / (max_val - min_val)
}))
C=Time
C=t(apply(C, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x-min_val)/(max_val-min_val)
}))

A=1-Er
Ef=1-C
# 1. Criação do data.frame com as 3 métricas
df <- data.frame(
  m = rep(m, 3),
  Valor = c(colMeans(Ef), colMeans(A), colMeans((2 * Ef * A) / (Ef + A))),
  Metrica = factor(
    rep(c("Eficiência", "Acurácia", "F-Score"), each = length(m)),
    levels = c("Eficiência", "Acurácia", "F-Score")
  )
)

# 2. Gráfico sem título e sem quadro na legenda
ggplot(df, aes(x = m, y = Valor, color = Metrica, shape = Metrica, linetype = Metrica)) +
  #geom_vline(xintercept = 100, linetype = "dashed", color = "gray60", linewidth = 0.7) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  
  scale_x_continuous(breaks = m) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  
  scale_color_manual(values = c("Eficiência" = "#0072B2", 
                                "Acurácia"   = "#D55E00", 
                                "F-Score"    = "#009E73")) +
  scale_shape_manual(values = c("Eficiência" = 19, "Acurácia" = 19, "F-Score" = 17)) +
  scale_linetype_manual(values = c("Eficiência" = "solid", "Acurácia" = "solid", "F-Score" = "dashed")) +
  
  labs(x = "m", y = NULL, color = NULL, shape = NULL, linetype = NULL) +
  
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14, color = "black"),
    
    legend.text  = element_text(size = 12),
    legend.position = "right",
    legend.background = element_blank(),     
    legend.box.background = element_blank(), 
    panel.grid.minor = element_blank()
  )
################################### CENARIO 2 ##################################
### b=5
load("~/tese/HMM/2Res_5.RData")
M=Time=b0=b1=b2=mu=phi=sigma=nu=matrix(nrow=100,ncol=10)
for(i in 1:100){
  for(j in 1:10){
    status=Results[[i]][[j]]$status
    if(status=="success"){
      M[i,j]=Results[[i]][[j]]$logL
      Time[i,j]=Results[[i]][[j]]$time
      b0[i,j]=Results[[i]][[j]]$summary[1,1]
      b1[i,j]=Results[[i]][[j]]$summary[2,1]
      b2[i,j]=Results[[i]][[j]]$summary[3,1]
      mu[i,j]=Results[[i]][[j]]$summary[4,1]
      phi[i,j]=Results[[i]][[j]]$summary[5,1]
      sigma[i,j]=Results[[i]][[j]]$summary[6,1]
      nu[i,j]=Results[[i]][[j]]$summary[7,1]
    }
  }
}

M=na.omit(M)
Time=na.omit(Time)
b0=na.omit(b0)
b1=na.omit(b1)
b2=na.omit(b2)
mu=na.omit(mu)
phi=na.omit(phi)
sigma=na.omit(sigma)
nu=na.omit(nu)

# Estimates analisys
m=seq(25, 250, by=25)
par(mfrow=c(1,2))
boxplot(phi[,1],phi[,2],phi[,3],phi[,4],phi[,5],
        phi[,6],phi[,7],phi[,8],phi[,9],phi[,10], main='phi', names=m)
abline(h=0.99, col='red', lty=2, lwd=2)
boxplot(sigma[,1],sigma[,2],sigma[,3],sigma[,4],sigma[,5],
        sigma[,6],sigma[,7],sigma[,8],sigma[,9],sigma[,10], main='sigma', names=m)
abline(h=0.15, col='red', lty=2, lwd=2)
par(mfrow=c(1,1))

S25=S50=S75=S100=S125=S150=S175=S200=S225=S250=matrix(0, nrow=7 ,ncol=5)
n=nrow(M)
for(k in 1:n){
  S25  = S25  + Results[[k]][[1]]$summary
  S50  = S50  + Results[[k]][[2]]$summary
  S75  = S75  + Results[[k]][[3]]$summary
  S100 = S100 + Results[[k]][[4]]$summary
  S125 = S125 + Results[[k]][[5]]$summary
  S150 = S150 + Results[[k]][[6]]$summary
  S175 = S175 + Results[[k]][[7]]$summary
  S200 = S200 + Results[[k]][[8]]$summary
  S225 = S225 + Results[[k]][[9]]$summary
  S250 = S250 + Results[[k]][[10]]$summary
}
S25=S25[,c(1,3,5)]/100
S50=S50[,c(1,3,5)]/100
S75=S75[,c(1,3,5)]/100
S100=S100[,c(1,3,5)]/100
S125=S125[,c(1,3,5)]/100
S150=S150[,c(1,3,5)]/100
S175=S175[,c(1,3,5)]/100
S200=S200[,c(1,3,5)]/100
S225=S225[,c(1,3,5)]/100
S250=S250[,c(1,3,5)]/100


# 1. Extração segura dos dados (base R)
list_dfs <- lapply(seq_along(Results), function(i) {
  lapply(seq_along(Results[[i]]), function(j) {
    res <- Results[[i]][[j]]
    if (!is.null(res) && res$status == "success") {
      summary_mat <- res$summary
      data.frame(
        sim = i,
        m = res$m,
        param = rownames(summary_mat),
        mean_val = summary_mat[, "mean"],
        q2.5 = summary_mat[, "2.5%"],
        q97.5 = summary_mat[, "97.5%"]
      )
    }
  })
})

# Empilha todas as listas em um único Data Frame
df_raw <- do.call(rbind, unlist(list_dfs, recursive = FALSE))
# 2. Agregação
df_summary <- aggregate(
  cbind(mean_val, q2.5, q97.5) ~ m + param, 
  data = df_raw, 
  FUN = mean, 
  na.rm = TRUE
)
# Renomeia colunas
colnames(df_summary) <- c("m", "param", "mean_est", "lower_ic", "upper_ic")
ordem_desejada <- c("b0", "b1", "b2", "mu", "phi", "sigma", "nu")
df_summary$param <- factor(df_summary$param, levels = ordem_desejada)
# 3. Data Frame com os valores verdadeiros reordenados
true_values <- data.frame(
  param = factor(ordem_desejada, levels = ordem_desejada),
  true_val = c(0.10, -0.01, -0.05, 0.10, 0.99, 0.15, 5)
)
# 4. Gráfico com ggplot2 (ncol = 3 garante 3 painéis por linha)
ggplot(df_summary, aes(x = factor(m), y = mean_est)) +
  geom_pointrange(aes(ymin = lower_ic, ymax = upper_ic), color = "#1f77b4", size = 0.5) +
  geom_hline(data = true_values, aes(yintercept = true_val), color = "red", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ param, scales = "free_y", ncol = 3) +
  labs(
    x = "Tamanho de m",
    y = "Estimativa (Média ± IC 95%)",
    title = "Convergência dos Parâmetros por tamanho de m"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )



# Eficient analisys
Er=abs(M-M[,ncol(M)])
Er=t(apply(Er, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x - min_val) / (max_val - min_val)
}))
C=Time
C=t(apply(C, 1, function(x){
  min_val <- min(x)
  max_val <- max(x)
  (x-min_val)/(max_val-min_val)
}))

A=1-Er
Ef=1-C
# 1. Criação do data.frame com as 3 métricas
df <- data.frame(
  m = rep(m, 3),
  Valor = c(colMeans(Ef), colMeans(A), colMeans((2 * Ef * A) / (Ef + A))),
  Metrica = factor(
    rep(c("Eficiência", "Acurácia", "F-Score"), each = length(m)),
    levels = c("Eficiência", "Acurácia", "F-Score")
  )
)

# 2. Gráfico sem título e sem quadro na legenda
ggplot(df, aes(x = m, y = Valor, color = Metrica, shape = Metrica, linetype = Metrica)) +
  #geom_vline(xintercept = 100, linetype = "dashed", color = "gray60", linewidth = 0.7) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  
  scale_x_continuous(breaks = m) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  
  scale_color_manual(values = c("Eficiência" = "#0072B2", 
                                "Acurácia"   = "#D55E00", 
                                "F-Score"    = "#009E73")) +
  scale_shape_manual(values = c("Eficiência" = 19, "Acurácia" = 19, "F-Score" = 17)) +
  scale_linetype_manual(values = c("Eficiência" = "solid", "Acurácia" = "solid", "F-Score" = "dashed")) +
  
  labs(x = "m", y = NULL, color = NULL, shape = NULL, linetype = NULL) +
  
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14, color = "black"),
    
    legend.text  = element_text(size = 12),
    legend.position = "right",
    legend.background = element_blank(),     
    legend.box.background = element_blank(), 
    panel.grid.minor = element_blank()
  )