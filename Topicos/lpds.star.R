lpds.star = function( h_T, a_T, theta_hat, yobs, kernel='epanechnikov'){
  
  # theta = (mu, phi_h, s_h, phi_a, s_a, v)
  theta = theta_hat
  M = dim( theta )[ 2 ]
  
  # h_new = mu + phi * ( h_T - mu ) + sigma_h * eps_h
  eps_h = rnorm( M, sd = sqrt(theta[3, ]) )
  h_new = theta[1, ] + theta[2, ] * ( h_T - theta[1, ] ) + eps_h
  # a_new = phi_a * a_T + eps_a
  eps_a = rnorm( M, sd = sqrt(theta[5, ]) )
  a_new = theta[4, ] * a_T + eps_a
  delta_new = a_new / sqrt(1 + a_new^2)
  k1 = sqrt(0.5 * theta[6,]) * gamma(0.5*(theta[6,]-1))/gamma(0.5*theta[6,])
  k2 = theta[6,] / (theta[6,] - 2)
  omega_new = 1/sqrt(k2 - 2/pi * delta_new^2 * k1^2)
  gamma_new = -sqrt(2/pi) * omega_new * delta_new * k1
  u_new = rgamma(M, 0.5*theta[6,], 0.5*theta[6,])
  w_new = qnorm( 0.5 * (runif(M) + 1) )
  mu_new = (gamma_new + omega_new*delta_new*w_new/sqrt(u_new))*exp(0.5*h_T)
  s_new = w_new * (1 - delta_new^2)*exp(0.5*h_T)/sqrt(u_new)
  # y_new = mu + a_T * exp{ h_T } + eps_y
  y_new = s_new * rnorm(M, mean = mu_new)
  
  # Calcula a densidade usando a função density do R
  dens = density(y_new, kernel = kernel)
  # Aproxima a densidade no ponto yobs
  dens_yobs = approx(dens$x, dens$y, xout = yobs)$y
  
  return( log(dens_yobs) ) 
}
