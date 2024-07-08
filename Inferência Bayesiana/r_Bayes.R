r.Bayes = function( a , theta ){
  
  N = dim( a )[1]
  T = dim( a )[2]
  if( length( theta ) != N ) print( 'Erro!' )
  L = rho = 0
  # Amostras L(theta, delta)
  for( n in 1:N ){
    for( m in 1:N ){
      a.piv = a[ n, ]
      L = L + ( theta[ m ] - sum( ( a.piv[2:T] - a.piv[1:(T-1)] )**2 ) )**2
    }
    rho = rho + L / N
  }
 return( rho / N ) 
}
