data {
  int<lower=0> T;
  vector[T] a;
  real<lower=0> lambda;
}
parameters{
  real<lower=0> xi;
}
transformed parameters{
 real delta = exp( xi );
}
model{
  
  // PC priori
  xi ~ weibull( 0.5, pow(lambda, -2 ) );
  
  for(t in 2:T) a[t] ~ normal( a[t-1], sqrt(xi) );
  
}
