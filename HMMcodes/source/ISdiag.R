quant = function(x, weights, probs=c(0.025, 0.5, 0.975)){
  
  ix = order(x)
  x_sorted = x[ix]
  weights_sorted = weights[ix]
  weights_sorted = weights_sorted / sum(weights_sorted)
  cumulative_weights = cumsum(weights_sorted)
  quants = numeric(length(probs))
  for(i in seq_along(probs)){
    index = which(cumulative_weights >= probs[i])[1]
    quants[i] = x_sorted[index]
  }
  
  return(quants)
}

ISdiag = function(Weigth, X, nu.lower=0, nu.upper=Inf, alpha){
  
  if(nu.lower==Inf){
    # svm-n
    p=function(parvect){
      
      beta=array(0,dim=3)
      beta[1]= parvect[1]
      beta[2]=(exp(parvect[2])-1)/(exp(parvect[2])+1)
      beta[3]=parvect[3]
      mu=parvect[4]
      phi = (exp(parvect[5])-1)/(exp(parvect[5])+1)
      sigma = exp(parvect[6])
      
      return(c(beta, mu, phi, sigma))
    }
  }else{
    if(nu.upper==Inf){
      p=function(parvect, nu.lower){
        
        beta=array(0,dim=3)
        beta[1]= parvect[1]
        beta[2]=(exp(parvect[2])-1)/(exp(parvect[2])+1)
        beta[3]=parvect[3]
        mu=parvect[4]
        phi = (exp(parvect[5])-1)/(exp(parvect[5])+1)
        sigma = exp(parvect[6])
        nu=exp(parvect[7])+nu.lower
        
        return(c(beta, mu, phi, sigma, nu))
      }
    }else{
      p=function(parvect, nu.lower, nu.upper, alpha){
        
        beta=array(0,dim=3)
        beta[1]= parvect[1]
        beta[2]=(exp(parvect[2])-1)/(exp(parvect[2])+1)
        beta[3]=parvect[3]
        mu=parvect[4]
        phi = (exp(parvect[5])-1)/(exp(parvect[5])+1)
        sigma = exp(parvect[6])
        
        #nu = (40*exp(parvect[7])+2)/(1+exp(parvect[7]))
        #alpha=0.1
        nu=0.5*((nu.upper-nu.lower)*tanh(0.5*alpha*parvect[7])+(nu.upper+nu.lower))
        
        return(c(beta, mu, phi, sigma, nu))
      }
    }
    
  }
  
  ### mean #####################################################################
  if(nu.lower==Inf){
    Thetas=t(apply(X=X, MARGIN=1, FUN=p))
  }else{
    if(nu.upper==Inf){
      Thetas=t(apply(X=X, MARGIN=1, FUN=p, nu.lower=nu.lower))
    }else{
      Thetas=t(apply(X=X, MARGIN=1, FUN=p, nu.lower=nu.lower, nu.upper=nu.upper))  
    }
  }
  wThetas = apply(X=Thetas, MARGIN = 2, FUN='*', Weigth)
  theta_hat = apply(X=wThetas, MARGIN = 2, sum)
  ### var ######################################################################
  Vars = t(apply(X=Thetas, 1, '-', theta_hat))
  wVars = apply(X=Vars, MARGIN = 2, FUN='*', Weigth)
  wVars = wVars^2
  Vars_hat = apply(X=wVars, MARGIN = 2, sum)
  ### quantiles ################################################################
  Quants = t(apply(Thetas, 2, quant, Weigth))
  Results = cbind(theta_hat, Vars_hat, Quants)
  colnames(Results) = c('mean', 'var', '2.5%', '50%', '97.5%')
  if(nu.lower==Inf){
    row.names(Results) = c('b0','b1','b2','mu','phi','sigma')
  }else{
    row.names(Results) = c('b0','b1','b2','mu','phi','sigma','nu')
  }
  Results = round(Results, 4)#signif(Results, 3)
  
  return(list(Results=Results, ess = diagis::ess(Weigth)))
}
