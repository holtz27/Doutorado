quantile = function(x, weights, probs=c(0.025, 0.5, 0.975)){
  
  ix <- order(x)
  x_sorted <- x[ix]
  weights_sorted <- weights[ix]
  weights_sorted <- weights_sorted / sum(weights_sorted)
  cumulative_weights <- cumsum(weights_sorted)
  quantiles <- numeric(length(probs))
  for (i in seq_along(probs)) {
    index <- which(cumulative_weights >= probs[i])[1]
    quantiles[i] <- x_sorted[index]
  }
  
  return(quantiles)
}
ISdiag = function(Weigth, X, knu=2){
  
  #if(sum(Weigth) != 1){
  #  warning('Weigth is unnormalized')
  #  Weigth = Weigth/sum(Weigth)
  #}
  
  p <- function(parvect, knu=knu){
    
    # Inicializar beta
    beta <- array(0, dim = 3)
    beta[1] <- parvect[1]
    if (parvect[2] > log_double_xmax) {
      beta[2] <- 1.0 - double_eps
    } else if (parvect[2] < -log_double_xmax) {
      beta[2] <- -1.0 + double_eps
    } else {
      beta[2] <- (exp(parvect[2]) - 1) / (exp(parvect[2]) + 1)
    }
    beta[3] <- parvect[3]
    
    mu <- parvect[4]
    if (parvect[5] > log_double_xmax) {
      phi <- 1.0 - double_eps
    } else if (parvect[5] < -log_double_xmax) {
      phi <- -1.0 + double_eps
    } else {
      phi <- (exp(parvect[5]) - 1) / (exp(parvect[5]) + 1)
    }
    
    if (parvect[6] > log_double_xmax) {
      sigma <- double_xmax
    } else if (parvect[6] < -log_double_xmax) {
      sigma <- double_xmin
    } else {
      sigma <- exp(parvect[6])
    }
    
    # Verificações para nu
    if (parvect[7] > log_double_xmax) {
      nu <- double_xmax
    } else if (parvect[7] < -log_double_xmax) {
      nu <-  double_eps
    } else {
      nu <- exp(parvect[7]) + knu
    }
    
    return(c(beta, mu, phi, sigma, nu))
  }
  ### mean #####################################################################
  Thetas = t(apply(X=X, MARGIN=1, FUN=p, knu=knu))
  wThetas = apply(X=Thetas, MARGIN = 2, FUN='*', Weigth)
  theta_hat = apply(X=wThetas, MARGIN = 2, sum)
  ### var ######################################################################
  Vars = t(apply(X=Thetas, 1, '-', theta_hat))
  wVars = apply(X=Vars, MARGIN = 2, FUN='*', Weigth)
  wVars = wVars^2
  Vars_hat = apply(X=wVars, MARGIN = 2, sum)
  ### quantiles ################################################################
  Quants = t(apply(Thetas, 2, quantile, Weigth))
  
  Results = cbind(theta_hat, Vars_hat, Quants)
  colnames(Results) = c('mean', 'var', '2.5%', '50%', '97.5%')
  row.names(Results) = c('b0','b1','b2','mu','phi','sigma','nu')
  Results = signif(Results, 3)
  
  return(list(Results=Results, ess = diagis::ess(Weigth)))
}
