res_sim = function(s , theta_vdd, med.abs = TRUE, digits = 4, names){
  
  T = length(s)
  conv = T
  rows = dim(s[[1]])[1]
  if(med.abs){
    theta = matrix(1, nrow = rows, ncol = 1)
  }else{
    theta = theta_vdd
    theta[theta == 0] = 1
  }
  
  theta_hat = err = amp = matrix(0, nrow=rows, ncol=T) 
  x1 = prob.cob = 0
  errors = NULL
  
  for(i in 1:T){
    
    if(sum(is.nan(s[[i]]) + is.infinite(s[[i]])) != 0){
      errors = cbind( errors, i ) 
      i = i + 1
    }else{ if(sum(abs( s[[ i ]][ , 'CD']) > 1.96 ) > 0){
      conv = conv-1
    }else{
      x1 = x1 + s[[i]]
      theta_hat[,i] = matrix(s[[ i ]][, 1], ncol = 1)
      err[,i] = (theta_hat[,i] - theta_vdd)/theta
    }
    }
  } 
  ##############################################################################
  indx = NULL
  for(t in 1:rows){
    x = boxplot(err[t,], plot = FALSE)
    out = x$out
    if(length(out)!=0) indx = cbind(indx, which(err[t,]==max(out)))
  }
  indx = unique(as.vector(indx))
  err = err[,-indx]
  ##############################################################################
  # Prob cob
  L = 1:T
  L = L[-indx]
  for(i in L){
    l1 = as.numeric(s[[i]][,3] < theta_vdd)
    l2 = as.numeric(s[[i]][,5] > theta_vdd)
    amp[,i] = t(s[[i]][,5]-s[[i]][,3])
    prob.piv = matrix(round(0.5*(l1 + l2),0), ncol=1)
    prob.cob = prob.cob + prob.piv
  }
  
  Data = cbind(theta_vdd,
               matrix(apply(theta_hat, MARGIN = 1, mean), ncol = 1),
               matrix(apply(err, MARGIN = 1, mean ), ncol = 1),
               matrix(apply(err^2, MARGIN = 1, mean ), ncol = 1),
               prob.cob/length(L), 
               apply(amp, 1, mean) )
  row.names(Data) = names #c('mu_h', 'phi_h', 's_h', 'phi_a', 's_a', 'v')
  colnames(Data) = c('Theta', 'Theta_hat', 'vies', 'reqm', 'prob.cob', 'amplitude')
  
  return(list(resumo = round(x1/conv, digits), 
              metricas = round(Data, digits),
              conv = conv,
              errors = errors,
              indx = indx)) 
}
