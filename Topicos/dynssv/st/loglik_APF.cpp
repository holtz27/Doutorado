// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

mat moms(vec a, vec h, double v, int N){
  
  mat M(N, 2);
  
  double invsqrt_U = sqrt(0.5*v)*tgamma(0.5*(v-1))/tgamma(0.5*v);
  double W = sqrt(2/M_PI);
  
  vec delta_t = a / sqrt(1 + pow(a, 2));
  double k1 = sqrt(0.5*v)*tgamma(0.5*(v-1))/tgamma(0.5*v);
  double k2 = v/(v-2);
  vec omega_t = 1/sqrt(k2-(2/M_PI) * pow(k1*delta_t, 2) );
  vec gamma_t = -sqrt(2/M_PI) * k1 * delta_t % omega_t;
  vec mu_t = gamma_t + W * omega_t % delta_t % exp(0.5*h); 
  mu_t *= invsqrt_U;
  vec sigma_t = omega_t % sqrt(1-pow(delta_t, 2)) % exp(0.5*h);
  sigma_t *= invsqrt_U;
  
  M.col(0) = mu_t;
  M.col(1) = sigma_t;
  
  return M;
}


// [[Rcpp::export]]
List loglik_APF(vec y, vec theta, int N){
  
  Function sample("sample.int");
  
  //arma::arma_rng::set_seed(42);
  
  int T = y.n_elem;
  double mu = theta[0];
  double phi = theta[1];
  double sh = theta[2];
  double sa = theta[3];
  double v = theta[4]; 
  
  mat h_par(N, T);
  mat a_par(N, T);
  mat w(N, T);
  
  vec h_par0 = randn( N, distr_param(mu, sh/sqrt(1-phi*phi)) );
  vec a_par0(N, fill::zeros);
  vec w0(N, fill::ones);
  w0 *= 1.0 / N;
  
  vec hbar(N), abar(N), pk(N);
  
  hbar = mu + phi*(h_par0 - mu);
  abar = a_par0;
  mat Ms = moms(abar, hbar, v, N);
  
  pk = normpdf(y(0), Ms.col(0), Ms.col(1)) % w0;
  pk /= sum(pk);
  
  IntegerVector indx = sample(N, N, true, pk);
  indx = as<arma::uvec>(indx);
  
  //new states
  vec X(N);
  for(int i=0; i < N; i++){
    X(i) = mu + phi*( hbar(indx(i)-1)-mu );
  }
  h_par.col(0) = X + sh * randn(N);
  
  for(int i=0; i < N; i++){
    X(i) = abar(indx(i)-1);
  }
  a_par.col(0) = X+sa*randn(N);
  
  Ms = moms(a_par.col(0), h_par.col(0), v, N);
  vec dens1 = log_normpdf(y(0), Ms.col(0), Ms.col(1));
    
  Ms = moms(abar, hbar, v, N);
  vec dens2 = log_normpdf(y(0), Ms.col(0), Ms.col(1));
    
  w.col(0) = exp(dens1-dens2);
  w.col(0) = w.col(0)/sum(w.col(0));

  double loglik = log( sum( exp(dens1) % w.col(0) ) );
  
  for(int t=1; t<T; t++){
    
    hbar = mu + phi*(h_par.col(t-1) - mu);
    abar = a_par.col(t-1);
    Ms = moms(abar, hbar, v, N);
    pk = normpdf(y(t), Ms.col(0), Ms.col(1)) % w.col(t-1);
    pk /= sum(pk);
    indx = sample(N, N, true, pk);
    indx = as<arma::uvec>(indx);
    
    //new states
    for(int i=0; i < N; i++){
      X(i) = mu + phi*( hbar(indx(i)-1)-mu );
    }
    h_par.col(t) = X+sh*randn(N);
    
    for(int i=0; i < N; i++){
      X(i) = abar(indx(i)-1);
    }
    a_par.col(t) = X+sa*randn(N);
    
    Ms = moms(a_par.col(t), h_par.col(t), v, N);
    dens1 = log_normpdf(y(t), Ms.col(0), Ms.col(1));
    Ms = moms(abar, hbar, v, N);
    dens2 = log_normpdf(y(t), Ms.col(0), Ms.col(1));
    w.col(t) = exp(dens1-dens2);
    w.col(t) = w.col(t)/sum(w.col(t));
    
    loglik += log( sum( exp(dens1) % w.col(t) ) );
    
  }
  
  return List::create(Named("loglik") = loglik,
                      Named("pred_h") = h_par,
                      Named("pred_a") = a_par,
                      Named("w") = w);
}