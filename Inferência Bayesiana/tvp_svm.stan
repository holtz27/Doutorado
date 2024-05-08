data {
  int<lower=0> T;
  vector[T] y;
  real a1;
  real<lower=0> lambda;
}
parameters{
  real tau;
  real mu;                      
  real<lower=0,upper=1> phiT;  
  real<lower=0> s_h;             
  real<lower=0> s_a;
  vector[T] h_std; // std log volatility time t               
  vector[T] a_std;
}
transformed parameters{
 real<lower=-1,upper=1> phi;
 real<lower=0> s2_h;
 real<lower=0> s2_a;
 vector[T] mu_t;
 phi = ( 2*phiT - 1 );
 s2_h = pow( s_h, 2 );
 s2_a = pow( s_a, 2 );
 
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
  // PC priori
  target += - 0.5 * log( s2_a ) - lambda * sqrt( s2_a );
  
  // Verossimilhança
  h_std ~ std_normal();
  a_std ~ std_normal();
  y ~ normal( mu_t, exp( 0.5 * h ) );
  
}

