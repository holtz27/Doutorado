data {
  int<lower=0> T;
  vector[T] y;
  real<lower=0> lambda1;
  real<lower=0> lambda2;
}
parameters{
  real mu_h;
  real mu_a;
  real<lower=0, upper=1> phiT_h;
  real<lower=0, upper=1> phiT_a;
  real<lower=0> s_h;
  real<lower=0> s_a;
  vector[T] h_std;
  vector[T] a_std;
  vector<lower=0>[T] W;
  vector<lower=0>[T] U;
  
  real<lower=2> v;
  real<lower=0> k;
}
transformed parameters{
  real<lower=-1,upper=1> phi_h;
  real<lower=-1,upper=1> phi_a;
  real<lower=0> s2_h;
  real<lower=0> s2_a;
  vector[T] a = a_std * s_a;
  vector[T] h = h_std * s_h;  
  real ls_a = log( s_a );
  phi_h = ( 2 * phiT_h - 1 );
  phi_a = ( 2 * phiT_a - 1 );
  s2_h = pow( s_h, 2 );
  s2_a = pow( s_a, 2 );
  
  h[1] /= sqrt(1 - phi_h * phi_h);
  a[1] /= sqrt(1 - phi_a * phi_a);
  h += mu_h;
  a += mu_a;
  for (t in 2:T){
    h[t] += phi_h * (h[t - 1] - mu_h);
    a[t] += phi_a * (a[t - 1] - mu_a);
  }
  
  vector<lower=-1, upper=1>[T] delta = a ./ sqrt(1 + square(a) );
  real k1 = sqrt(0.5 * v) * tgamma(0.5*(v-1)) / tgamma(0.5 * v);
  real k2 = v / (v-2);
  vector[T] omega = 1 ./ sqrt(k2 - 2 * square( delta * k1 ) / pi()); 
  vector[T] mean_st = - sqrt(2 / pi()) * k1 * delta .* omega;
  vector[T] mu_t = mean_st + omega .* delta .* W .* exp(0.5 * h) ./ sqrt( U );
  vector[T] sigma_t = omega .* sqrt(1 - square ( delta )) .* exp(0.5 * h) ./ sqrt( U );
  
}
model {
  // Prioris h
  mu_h ~ normal( 0, sqrt(10) );
  phiT_h ~ beta( 20, 1.5 );
  s2_h ~ inv_gamma( 2.5, 0.025 );
  
  // Prioris a
  mu_a ~ normal( 0, sqrt(10) );
  
  //priori1
  s2_a ~ weibull( 0.5, pow(lambda2, -2 ) );
  
  //prior2
  //s2_a ~ inv_gamma( 2.5, 0.025 );
  
  //prior3
  //k ~ gamma(0.1, 0.1);
  //s_a ~ exponential( k );
  
  target += - lambda1 * sqrt(1 - phi_a) - 0.5 * log(1 - phi_a);
  
  //tails
  v ~ gamma(12, .8);
  
  // model
  U ~ gamma(0.5 * v, 0.5 * v);
  h_std ~ std_normal();
  a_std ~ std_normal();
  y ~ normal( mu_t, sigma_t );
  target += - 0.5 * square( W );
  
}
