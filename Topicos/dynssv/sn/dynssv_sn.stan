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
}
transformed parameters{
  real<lower=-1,upper=1> phi_h;
  real<lower=-1,upper=1> phi_a;
  real<lower=0> s2_h;
  real<lower=0> s2_a;
  vector[T] a = a_std * s_a;
  vector[T] h = h_std * s_h;  
  phi_h = ( 2 * phiT_h - 1 );
  phi_a = ( 2 * phiT_a - 1 );
  s2_h = pow( s_h, 2 );
  s2_a = pow( s_a, 2 );
  
  h[1] /= sqrt(1 - phi_h * phi_h);
  a[1] += a1;
  h += mu;
  for (t in 2:T){
    h[t] += phi_h * (h[t - 1] - mu);
    a[t] += phi_a * a[t-1];
  } 
  
  vector<lower=-1, upper=1>[T] delta = a ./ sqrt(1 + square(a) );
  vector[T] omega = 1 ./ sqrt(1 - 2 * square( delta ) / pi());
  vector[T] mean_sn = - sqrt(2 / pi()) * delta .* omega;
  vector[T] mu_t = mean_sn + omega .* delta .* W .* exp(0.5 * h);
  vector[T] sigma_t = omega .* sqrt(1 - square ( delta )) .* exp(0.5 * h);
  
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
