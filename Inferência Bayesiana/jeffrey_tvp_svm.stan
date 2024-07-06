data {
  int<lower=0> T;
  vector[T] y;
}
parameters{
  real tau;
  real mu;                      
  real<lower=0,upper=1> phiT;  
  real<lower=0> s_h;             
  real<lower=0> s_a;
  vector[T] h_std; 
  real a1;
  vector[T] a_std;
}
transformed parameters{
 real<lower=-1,upper=1> phi;
 real<lower=0> s2_h;
 real<lower=0> s2_a;
 real log_s2_a;
 vector[T] mu_t;
 phi = ( 2 * phiT - 1 );
 s2_h = pow( s_h, 2 );
 s2_a = pow(s_a, 2);
 log_s2_a = log( s2_a );
 
 vector[T] h = h_std * s_h;  // now h ~ normal(0, sigma)
 vector[T] a = a_std * s_a;
  h[1] /= sqrt(1 - phi * phi);  // rescale h[1]
  h += mu;
  a[1] = a1;
  mu_t[1] = tau + a[1] * exp( h[1] );
  for (t in 2:T){
    h[t] += phi * (h[t - 1] - mu);
    a[t] += a[t - 1];
    mu_t[t] = tau + a[t] * exp( h[t] );
  }
}
model {
  // Prioris
  tau ~ normal( 0, sqrt(10) );
  mu ~ normal( 0, sqrt(10) );
  phiT ~ beta( 20, 1.5 );
  s2_h ~ inv_gamma( 2.5, 0.025 );
  a1 ~ normal( -1.0, 10 );
  // Jeffrey
  target += -log( s2_a );
  
  // model
  h_std ~ std_normal();
  a_std ~ std_normal();
  y ~ normal( mu_t, exp( 0.5 * h ) );
  
}
