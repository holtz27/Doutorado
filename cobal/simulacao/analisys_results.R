source( 'https://raw.githubusercontent.com/holtz27/Doutorado/refs/heads/main/source/res_sim.R' )
load("~/cobal/simulation/sim_xi_0.RData")

s1 = s2 = list()
b1 = b2 = eqm1 = eqm2 = matrix(0, nrow = 5, ncol=1)
N = length( result )
for(i in 1:N){
  s1[[ i ]] = result[[ i ]]$summary1[ c(1:5, 7), ]
  s2[[ i ]] = result[[ i ]]$summary2[ c(1:5, 7), ]
  b1 = b1 + result[[ i ]]$err1
  b2 = b2 + result[[ i ]]$err2
  eqm1 = eqm1 + result[[ i ]]$err1^2
  eqm2 = eqm2 + result[[ i ]]$err2^2
} 

# Parameter estimation
res_sim( s1, theta_vdd )$metricas
res_sim( s2, theta_vdd )$metricas

# Volatilitie Prediction
res.h = cbind( b1 / N, sqrt(eqm1 / N), b2 / N, sqrt(eqm2 / N) )
colnames(res.h) = c('bias1', 'reqm1', 'bias2', 'reqm2')
round( res.h, 3 )



