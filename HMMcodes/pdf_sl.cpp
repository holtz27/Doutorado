// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include<Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_integration.h>

// set up struct that contains an Rcpp Function, this will be used in the function provided to gsl
struct my_f_params {Rcpp::Function G;};
gsl_integration_workspace *work_ptr =
  gsl_integration_workspace_alloc (1000);

/*********************************************************************************/
struct my_sl_params {double x;double mu; double sigma; double nu;};
/*********************************************************************************/
double my_function_sl (double u, void *p){
  struct my_sl_params *params=(struct my_sl_params *)p;
  double x=(params->x);
  double mu = (params->mu);
  double sigma = (params->sigma);
  double nu=(params->nu);
  double f=(1.0/sqrt(M_2PI))*(sqrt(u)/sigma)*exp((-u*0.5/(sigma*sigma))*(x-mu)*(x-mu))*nu*exp((nu-1)*log(u));
  return  f;
}
/*********************************************************************************/

// [[Rcpp::export]]
Rcpp::NumericVector pdf_sl(Rcpp::NumericVector  y,double mu, double sigma, double nu){
  double lower_limit=0.0;
  double upper_limit =1.0;
  double abs_error = 1.0e-8;
  double rel_error = 1.0e-8;
  double result;
  double error;
  
  int nn=y.size();
  int kk=0;
  Rcpp::NumericVector res(nn);
  
  gsl_function G;
  
  for (kk=0;kk<nn;kk++){
    
    struct my_sl_params params={y[kk],mu,sigma,nu};
    G.function =&my_function_sl;
    G.params=&params;
    gsl_integration_qag(&G, lower_limit, upper_limit, abs_error, rel_error, 1000, GSL_INTEG_GAUSS21, work_ptr, &result, &error);
    res[kk]=result;
    
  }
  return res;
}
/*********************************************************************************/

