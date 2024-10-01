data {
  int<lower=0> T;
  vector[T] y;
  real<lower=0> lambda1;
  real<lower=0> lambda2;
}
parameters{
  real mu;                      
  real<lower=0, upper=1> phiT_h;
  real<lower=0, upper=1> phiT_a;
  real<lower=0> s_h;
  real<lower=0> s_a;
  vector[T] h_std;
  vector[T] a_std;
  real a1;
  vector<lower=0>[T] W;
  vector<lower=-1, upper=1>[T] delta;
}
transformed parameters{
  real<lower=-1,upper=1> phi_h;
  real<lower=-1,upper=1> phi_a;
  real<lower=0> s2_h;
  real<lower=0> s2_a;
  
  vector[T] a = a_std * s_a;
  vector[T] h = h_std * s_h;  
  vector[T] mu_t;
  vector[T] sigma_t;
  vector[T] omega = 1 ./ sqrt(1 - 2 * square( delta ) / pi());
  vector[T] mean_sn = - sqrt(2 / pi()) * delta .* omega;
  
  phi_h = ( 2 * phiT_h - 1 );
  phi_a = ( 2 * phiT_a - 1 );
  s2_h = pow( s_h, 2 );
  s2_a = pow( s_a, 2 );
  
  h[1] /= sqrt(1 - phi_h * phi_h);
  a[1] += a1;
  h += mu;
  
  //omega[1] = 1 / sqrt(1 - 2 * square( delta[1] ) / pi());
  //mean_sn[1] = - sqrt(2 / pi()) * delta[1] * omega[1];
  mu_t[1] = mean_sn[1] + omega[1] * delta[1] * W[1] * exp(0.5 * h[1]) ;
  sigma_t[1] = omega[1] * sqrt(1 - square ( delta[1] )) * exp(0.5 * h[1]) ;
  
  for (t in 2:T){
    h[t] += phi_h * (h[t - 1] - mu);
    a[t] += phi_a * a[t-1];
    //omega[t] = 1 / sqrt(1 - 2 * square( delta[t] ) / pi());
    //mean_sn[t] = - sqrt(2 / pi()) * delta[t] * omega[t];
    mu_t[t] = mean_sn[t] + omega[t] * delta[t] * W[t] * exp(0.5 * h[t]);
    sigma_t[t] = omega[t] * sqrt(1 - square ( delta[t] )) * exp(0.5 * h[t]) ;
  }
}
model {
  // Prioris h
  mu ~ normal( 0, sqrt(10) );
  phiT_h ~ beta( 20, 1.5 );
  s2_h ~ inv_gamma( 2.5, 0.025 );
  
  // Prioris a
  a1 ~ normal( 0, 1 );
  s2_a ~ weibull( 0.5, pow(lambda2, -2 ) );
  //s2_a ~ inv_gamma( 4.5, 0.065 );
  
  // model
  h_std ~ std_normal();
  a_std ~ std_normal();
  y ~ normal( mu_t, sigma_t );
  target += - 0.5 * square( W );
  target += - lambda1 * sqrt(1 - phi_a) - 0.5 * log(1 - phi_a);

}
