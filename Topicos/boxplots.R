load("~/Doutorado/Pesquisa/Papers/Topicos/DS-SV/Code/st/simulation/simulacao_sa_0.001.RData")
N = length(result)
err1.1 = err2.1 = err3.1 = err4.1 = numeric(N)
for(i in 1:N){
  err1.1[i] = result[[ i ]]$summary1[ 'ls_a', 1 ]
  err2.1[i] = result[[ i ]]$summary2[ 'ls_a', 1 ]
  err3.1[i] = result[[ i ]]$summary3[ 'ls_a', 1 ]
  #err4.1[i] = result[[ i ]]$summary4[ 's_a', 1 ]
}

load("~/Doutorado/Pesquisa/Papers/Topicos/DS-SV/Code/st/simulation/simulacao_sa_0.05.RData")
err1.2 = err2.2 = err3.2 = err4.2 = numeric(N)
for(i in 1:N){
  err1.2[i] = result[[ i ]]$summary1[ 'ls_a', 1 ]
  err2.2[i] = result[[ i ]]$summary2[ 'ls_a', 1 ]
  err3.2[i] = result[[ i ]]$summary3[ 'ls_a', 1 ]
  #err4.2[i] = result[[ i ]]$summary4[ 's_a', 1 ]
}

load("~/Doutorado/Pesquisa/Papers/Topicos/DS-SV/Code/st/simulation/simulacao_sa_0.1.RData")
err1.3 = err2.3 = err3.3 = err4.3 = numeric(N)
for(i in 1:N){
  err1.3[i] = result[[ i ]]$summary1[ 'ls_a', 1 ]
  err2.3[i] = result[[ i ]]$summary2[ 'ls_a', 1 ]
  err3.3[i] = result[[ i ]]$summary3[ 'ls_a', 1 ]
  #err4.3[i] = result[[ i ]]$summary3[ 's_a', 1 ]
}

par(mfrow = c(1, 3))
boxplot(err1.1, err2.1, err3.1, #err4.1,
        names = c('IG', 'PCP', 'EXP'), 
        main = expression(sigma[a] == 0.001), pch = 19, cex = 1.5,
        cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.2)
abline(h=log(0.001))
boxplot(err1.2, err2.2, err3.2, #err4.2,
        names = c('IG', 'PCP', 'EXP'),
        main = expression(sigma[a] == 0.05), pch = 19, cex = 1.5,
        cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.2)
abline(h=log(0.05))
boxplot(err1.3, err2.3, err3.3, #err4.3, 
        names = c('IG', 'PCP', 'EXP'),
        main = expression(sigma[a] == 0.1), pch = 19, cex = 1.5,
        cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.2)
abline(h=log(0.1))
par(mfrow = c(1, 1))

