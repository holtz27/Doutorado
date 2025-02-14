// [[Rcpp::depends(RcppArmadillo, RcppGSL)]]

#include <RcppArmadillo.h>
#include <gsl/gsl_sf_gamma.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pdf_s2(NumericVector y, double nu){
  int n = y.size();
  NumericVector result(n);
  
  double a = nu / sqrt(2 * M_PI);
  double s = nu + 0.5;
  double z, x;
  double tol = 1e-8;
  
  for(int i = 0; i < n; i++){
    
    x = 0.5 * y[i] * y[i];
    if(fabs(x) < tol){
      z = a / s;
    }else{
      z = a*pow(1.0 / x, s)*gsl_sf_gamma_inc_P(s, x);
      z *= gsl_sf_gamma(s);
      //z = a*pow(1.0 / x, s)*(gsl_sf_gamma(s) - gsl_sf_gamma_inc(s, x));
    }
    result[i] = z;
    
  }
  
  return result;
}
