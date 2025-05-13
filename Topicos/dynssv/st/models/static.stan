data {
  int<lower=0> T;
  vector[T] y;
}
parameters{
  real mu;
  real<lower=0, upper=1> phiT_h;
  real<lower=0> s_h;
  vector[T] h_std;
  vector<lower=0>[T] W;
  vector<lower=0>[T] U;
  real a;
  real<lower=2> v;
  real<lower=0> k;
}
transformed parameters{
  real<lower=-1,upper=1> phi_h;
  real<lower=0> s2_h;
  vector[T] h = h_std*s_h;
  phi_h = (2*phiT_h - 1);
  s2_h = pow(s_h, 2);
  
  h[1] /= sqrt(1 - phi_h*phi_h);
  h += mu;
  for (t in 2:T){
    h[t] += phi_h*(h[t-1]-mu);
  }

  real<lower=-1, upper=1> delta = a / sqrt(1 + square(a));
  real k1 = sqrt(0.5 * v) * tgamma(0.5*(v-1)) / tgamma(0.5 * v);
  real k2 = v / (v-2);
  real omega = 1 / sqrt(k2 - 2 * square( delta * k1 ) / pi()); 
  real gamma = -sqrt(2 / pi()) * k1 * delta * omega;
  vector[T] mu_t = gamma + omega * delta * W .* exp(0.5 * h) ./ sqrt( U );
  vector<lower=0>[T] sigma_t = omega * sqrt(1 - square(delta)) * exp(0.5 * h) ./ sqrt( U );
  
}
model {
  // Prioris h
  mu ~ normal( 0, sqrt(10) );
  phiT_h ~ beta( 20, 1.5 );
  s2_h ~ inv_gamma( 2.5, 0.025 );
  
  // Prioris a
  k ~ gamma(0.1, 0.1);
  a ~ double_exponential(0, 1/k);

  //tails
  v ~ gamma(12, .8);
  
  // model
  U ~ gamma(0.5 * v, 0.5 * v);
  h_std ~ std_normal();
  y ~ normal( mu_t, sigma_t );
  target += - 0.5 * square( W );
  
}
