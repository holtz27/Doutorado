data{
  int<lower=0> T;
  vector[T] y;
  real<lower=0> lambda;
}
parameters{
  real mu;
  real<lower=0, upper=1> phiT_h;
  real<lower=0> ka;
  real a1;
  real<lower=0> s_h;
  real<lower=0> s_a;
  vector[T] h_std;
  vector[T] a_std;
  vector<lower=0>[T] W;
  vector<lower=0,upper=1>[T] U;
  real<lower=1> v;
}
transformed parameters{
  real<lower=-1,upper=1> phi_h;
  real<lower=0> s2_h;
  vector[T] a = a_std*s_a;
  vector[T] h = h_std*s_h;
  phi_h = (2*phiT_h - 1);
  s2_h = pow(s_h, 2);
  
  h[1] /= sqrt(1-phi_h*phi_h);
  h += mu;
  a[1] += a1;
  for (t in 2:T){
    h[t] += phi_h * (h[t-1]-mu);
    a[t] += a[t-1];
  }

  vector<lower=-1, upper=1>[T] delta = a./sqrt(1+square(a));
  real k1 = 2*v/(2*v-1);
  real k2 = v/(v-1);
  vector[T] omega = 1 ./sqrt(k2 - 2*square(delta*k1)/pi()); 
  vector[T] mean_st = -sqrt(2/pi())*k1*delta.*omega;
  vector[T] mu_t = mean_st + omega.*delta.*W.*exp(0.5*h)./sqrt(U);
  vector[T] sigma_t = omega.*sqrt(1-square(delta)).*exp(0.5*h)./sqrt(U);
}
model{
  // Prioris h
  mu ~ normal(0, sqrt(10));
  phiT_h ~ beta(20, 1.5);
  s2_h ~ inv_gamma(2.5, 0.025);
  
  // Prioris a
  ka ~ gamma(0.1, 0.1);
  a1 ~ double_exponential(0, 1/ka);
  s_a ~ exponential(lambda);
  
  //tails
  v ~ gamma(0.08, 0.04);
  
  // model
  U ~ beta(v, 1.0);
  h_std ~ std_normal();
  a_std ~ std_normal();
  y ~ normal(mu_t, sigma_t);
  target += -0.5*square(W);
}
