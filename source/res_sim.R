res_sim = function( s , theta_vdd, med.abs = TRUE, digits = 4 ){
  
  T = length( s )
  conv = T
  rows = dim( s[[1]] )[ 1 ]
  if( med.abs ){
    theta = matrix(1, nrow = rows, ncol = 1)
  }else{
    theta = theta_vdd
    theta[ theta == 0 ] = 1
  }
  
  err = matrix(0, nrow = rows, ncol = T) 
  x1 = prob.cob = 0
  errors = NULL
  
  for( i in 1:T ){
    
    if( sum( is.nan( s[[ i ]] ) + is.infinite( s[[ i ]] ) ) != 0 ){
      errors = cbind( errors, i ) 
      i = i + 1
    }else{ if( sum( abs( s[[ i ]][ , 'CD'] ) > 1.96 ) > 0 ){
      conv = conv - 1
    }else{
      x1 = x1 + s[[ i ]]
      err[ ,i ] = ( matrix(s[[ i ]][, 1], ncol = 1) - theta_vdd ) / theta
      # prob cob
      l1 = as.numeric( s[[ i ]][, 3] < theta_vdd )
      l2 = as.numeric( s[[ i ]][, 5] > theta_vdd )
      prob.piv = matrix( round( 0.5 * (l1 + l2), 0 ), ncol = 1)
      prob.cob = prob.cob + prob.piv
    }
    }
  } 
  ##############################################################################
  indx = NULL
  for(t in 1:rows){
    x = boxplot(err[t, ], plot = FALSE)
    out = x$out
    if(length(out)!=0) indx = cbind(indx, which(err[t, ]==max(out)))
  }
  indx = unique(as.vector(indx))
  err = err[, -indx]
  ##############################################################################
  Data = cbind( matrix( apply( err, MARGIN = 1, mean ), ncol = 1 ),
                matrix( apply( err^2, MARGIN = 1, mean ), ncol = 1 ),
                prob.cob / conv )
  #row.names( Data ) = c('mu_h', 'phi_h', 's_h', 'mu_a', 'phi_a', 'ls_a', 'v')
  colnames( Data ) = c('vies', 'reqm', 'prob.cob')
  
  return( list( resumo = round(x1 / conv, digits), 
                metricas = round(Data, digits),
                conv = conv,
                errors = errors,
                indx = indx) 
  ) 
}
