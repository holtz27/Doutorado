data{
  int<lower=0> T;
  vector[T] y;
  real y0;
}
parameters{
  real b0;
  real<lower=-1,upper=1>b1;
  real b2;
  real mu;                      
  real<lower=-1,upper=1> phi;  
  real<lower=0> s_h;             
  vector[T] h_std;
  vector<lower=0,upper=1>[T] l;
  real<lower=0> v;
}
transformed parameters{
 real<lower=0,upper=1> b1T=0.5*(b1+1);
 real<lower=0,upper=1> phiT=0.5*(phi+1);
 real<lower=0> s_hT=pow(s_h,2);
 
 vector[T] mu_t;
 vector[T] h = h_std*s_h;  
 h[1] /= sqrt(1-phi*phi);  
 h += mu;
 mu_t[1]=b0+b1*y0+b2*exp(h[1]);
 for(t in 2:T){
   h[t] += phi*(h[t-1]-mu);
   mu_t[t] = b0+b1*y[t-1]+b2*exp(h[t]);
 }
 vector<lower=0>[T] sigma_t=exp(0.5*h)./sqrt(l);
}
model {
  // Prioris
  b0 ~ normal(0, sqrt(10));
  b1T ~ beta(5.0, 1.5);
  b2 ~ normal(0, sqrt(10));
  mu ~ normal(0, sqrt(10));
  phiT ~ normal(4.5, 10);
  s_hT ~ inv_gamma(2.5, 0.025);
  v ~ gamma(0.08,0.04);
  
  // model
  l ~ beta(v, 1);
  h_std ~ std_normal();
  y ~ normal(mu_t, sigma_t);
}
