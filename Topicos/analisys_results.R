source('https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/source/res_sim.R')
#load("~/Doutorado/Pesquisa/Papers/Topicos/DS-SV/Code/st/simulation/simulacao_sa_0.001.RData")
#load("~/Doutorado/Pesquisa/Papers/Topicos/DS-SV/Code/st/simulation/simulacao_sa_0.05.RData")
load("~/Doutorado/Pesquisa/Papers/Topicos/DS-SV/Code/st/simulation/simulacao_sa_0.1.RData")

s1 = s2 = s3 = list()
N = length( result )
for(i in 1:N){
  s1[[ i ]] = result[[ i ]]$summary1[ c(1:5, 8), ]
  s2[[ i ]] = result[[ i ]]$summary2[ c(1:5, 8), ]
  s3[[ i ]] = result[[ i ]]$summary3[ c(1:5, 8), ]
} 

# Parameter estimation
res_sim( s1, theta_vdd, med.abs = TRUE )$metricas
res_sim( s1, theta_vdd, med.abs = TRUE )$indx

res_sim( s2, theta_vdd, med.abs = TRUE )$metricas
res_sim( s2, theta_vdd, med.abs = TRUE )$indx

res_sim( s3, theta_vdd, med.abs = TRUE )$metricas
res_sim( s3, theta_vdd, med.abs = TRUE )$indx

# resumo
res_sim( s1, theta_vdd, med.abs = TRUE )$resumo
res_sim( s2, theta_vdd, med.abs = TRUE )$resumo
res_sim( s3, theta_vdd, med.abs = TRUE )$resumo
