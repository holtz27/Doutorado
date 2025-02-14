// [[Rcpp::depends(RcppArmadillo, RcppGSL)]]

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>

using namespace Rcpp;
using namespace arma;

struct my_vg_params { double x; double mu; double sigma; double nu; };
double my_function_vg(double u, void *p) {
  struct my_vg_params *params = (struct my_vg_params *)p;
  double x = (params->x);
  double mu = (params->mu);
  double sigma = (params->sigma);
  double nu = (params->nu);
  double f = (1.0 / sqrt(M_2PI)) * (sqrt(u) / sigma) * 
    exp((-u * 0.5 / (sigma * sigma)) * (x - mu) * (x - mu)) * 
    exp((nu / 2) * log(nu / 2)) * exp(-0.5 * nu / u) * 
    exp(-gsl_sf_lngamma(0.5*nu)) *
    //exp(-lgamma(0.5 * nu)) * 
    exp(-(0.5 * nu + 1) * log(u));
  return f;
}

// [[Rcpp::export]]
NumericVector pdf_vg(NumericVector y, double mu, double sigma, double nu){
  
  int nn = y.size();
  NumericVector res(nn);
  
  if(nu < 100){
    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(1000);
    double lower_limit = 0.0;
    double abs_error = 1.0e-4;
    double rel_error = 1.0e-4;
    double result, error;
    
    gsl_function G;
    G.function = &my_function_vg;
    
    for(int i = 0; i < nn; i++){
      
      struct my_vg_params params = { y[i], mu, sigma, nu };
      G.params = &params;
      gsl_integration_qagiu(&G, lower_limit, abs_error, rel_error, 1000, 
                            work_ptr, &result, &error);
      res[i] = result;
      
    }
    gsl_integration_workspace_free(work_ptr);
    
  }else{
    res = dnorm(y);
  }
  
 return res;
}
