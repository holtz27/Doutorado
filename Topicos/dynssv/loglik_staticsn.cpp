// [[Rcpp::depends(RcppArmadillo, RcppGSL)]]

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace Rcpp;
using namespace arma;

IntegerVector sample_int(int size, vec prob){
  int n = prob.n_elem;
  
  if (n == 0 || size < 0) {
    Rcpp::stop("Invalid size or empty probability vector");
  }
  
  double sum_prob = arma::sum(prob);
  if (sum_prob <= 0.0) {
    Rcpp::stop("Sum of probabilities must be positive");
  }
  
  // Inicializa o gerador de números aleatórios do GSL
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, time(NULL));  // semente (pode substituir por valor fixo para reprodutibilidade)
  
  // Cria a estrutura de aliasing da GSL
  gsl_ran_discrete_t* g = gsl_ran_discrete_preproc(n, prob.memptr());
  
  // Gera as amostras
  Rcpp::IntegerVector result(size);
  for (int i = 0; i < size; ++i) {
    result[i] = gsl_ran_discrete(r, g); 
  }
  
  // Libera memória
  gsl_ran_discrete_free(g);
  gsl_rng_free(r);
  
  return result;
}

mat moms(double a, vec h, int N){
  
  mat M(N, 2);
  
  double W = sqrt(2/M_PI);
  
  double delta_t = a / sqrt(1 + pow(a, 2));
  double k1 = 1;
  double k2 = 1;
  double omega_t = 1/sqrt(k2-(2/M_PI) * pow(k1*delta_t, 2) );
  double gamma_t = -sqrt(2/M_PI) * k1 * delta_t * omega_t;
  vec mu_t = gamma_t + W * omega_t * delta_t * exp(0.5*h); 
  vec sigma_t = omega_t * sqrt(1-pow(delta_t, 2)) * exp(0.5*h);
  
  M.col(0) = mu_t;
  M.col(1) = sigma_t;
  
  return M;
}

// [[Rcpp::export]]
List loglik_staticsn(vec y, vec theta, int N){
  
  //Function sample("sample.int");
  
  int T = y.n_elem;
  double mu = theta[0];
  double phi = theta[1];
  double sh = theta[2];
  double a1 = theta[3];
  
  mat h_par(N, T);
  mat w(N, T);
  
  vec h_par0 = randn( N, distr_param(mu, sh/sqrt(1-phi*phi)) );
  vec w0(N, fill::value(1.0/N));
  
  vec hbar(N, fill::zeros), abar(N, fill::zeros), pk(N, fill::zeros);
  
  hbar = mu + phi*(h_par0 - mu);
  mat Ms = moms(a1, hbar, N);
  
  pk = log_normpdf(y(0), Ms.col(0), Ms.col(1)) + log(w0);
  //double max_pk = pk.max();
  vec prob = exp(pk);
  //prob = prob / sum(prob);
  
  //IntegerVector indx = sample(N, N, true, prob);
  IntegerVector indx = sample_int(N, prob);
  indx = as<arma::uvec>(indx);
  
  //new states
  vec X(N);
  for(int i=0; i < N; i++){
    X(i) = mu + phi*( hbar(indx(i))-mu );
  }
  h_par.col(0) = X + sh * randn(N);
  
  Ms = moms(a1, h_par.col(0), N);
  vec dens1 = log_normpdf(y(0), Ms.col(0), Ms.col(1));
  
  Ms = moms(a1, hbar, N);
  vec dens2 = log_normpdf(y(0), Ms.col(0), Ms.col(1));
  
  w.col(0) = exp(dens1-dens2);
  double sum_w = sum(w.col(0));
  if( sum_w > 0){
    w.col(0) = w.col(0)/sum_w;  
  }else{
    stop("Error normalize weight");
  }
  
  double loglik = log( sum( exp(dens1) % w.col(0) ) );
  
  for(int t=1; t<T; t++){
    
    hbar = mu + phi*(h_par.col(t-1) - mu);
    Ms = moms(a1, hbar, N);
    
    
    pk = log_normpdf(y(t), Ms.col(0), Ms.col(1)) + log(w.col(t-1));
    prob = exp(pk);
    prob = prob / sum(prob);
    
    indx = sample_int(N, prob);
    indx = as<arma::uvec>(indx);
    
    //new states
    for(int i=0; i < N; i++){
      X(i) = mu + phi*( hbar(indx(i))-mu );
    }
    h_par.col(t) = X+sh*randn(N);
    
  
    Ms = moms(a1, h_par.col(t), N);
    dens1 = log_normpdf(y(t), Ms.col(0), Ms.col(1));
    Ms = moms(a1, hbar, N);
    dens2 = log_normpdf(y(t), Ms.col(0), Ms.col(1));
    w.col(t) = exp(dens1-dens2);
    
    sum_w = sum(w.col(t));
    if( sum_w > 0){
      w.col(t) = w.col(t)/sum_w;  
    }else{
      stop("Error normalize weight");
    }
    
    // Resample 
    indx= sample_int(N, w.col(t));
    indx = as<arma::uvec>(indx);
    for(int i=0; i < N; i++){
      h_par(i, t) = h_par(indx(i), t);
    }
    // Reset pesos
    w.col(t).fill(1.0 / N);
    
    loglik += log( sum( exp(dens1) % w.col(t) ) );
    
  }
  
  return List::create(Named("loglik") = loglik,
                      Named("pred_h") = h_par,
                      Named("w") = w);
}