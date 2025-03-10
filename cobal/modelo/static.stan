data {
  int<lower=0> T;
  vector[T] y;
}
parameters{
  real b;
  real mu;                      
  real<lower=0,upper=1> phiT;  
  real<lower=0> s_h;             
  vector[T] h_std; 
  real a;
}
transformed parameters{
 real<lower=-1,upper=1> phi;
 real<lower=0> s2_h;
 vector[T] mu_t;
 phi = ( 2 * phiT - 1 );
 s2_h = pow( s_h, 2 );
 
 vector[T] h = h_std * s_h;  // now h ~ normal(0, sigma)
 h[1] /= sqrt(1 - phi * phi);  // rescale h[1]
 h += mu;
 mu_t[1] = b + a * exp( h[1] );
 for (t in 2:T){
   h[t] += phi * (h[t - 1] - mu);
   mu_t[t] = b + a * exp( h[t] );
 }
}
model {
  // Prioris
  b ~ normal( 0, sqrt(10) );
  mu ~ normal( 0, sqrt(10) );
  phiT ~ beta( 20, 1.5 );
  s2_h ~ inv_gamma( 2.5, 0.025 );
  a ~ normal( 0, sqrt(10) );
  
  // model
  h_std ~ std_normal();
  y ~ normal( mu_t, exp( 0.5 * h ) );
  
}
