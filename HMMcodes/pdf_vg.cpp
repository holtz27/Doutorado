// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include<Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_integration.h>

// set up struct that contains an Rcpp Function, this will be used in the function provided to gsl
struct my_f_params {Rcpp::Function G;};
gsl_integration_workspace *work_ptr =
  gsl_integration_workspace_alloc (1000);

/******************************************************************/
double lgamma(double x){
  double x0,x2,xp,gl,gl0;
  int n,k;
  static double a[] = {
    8.333333333333333e-02,
    -2.777777777777778e-03,
    7.936507936507937e-04,
    -5.952380952380952e-04,
    8.417508417508418e-04,
    -1.917526917526918e-03,
    6.410256410256410e-03,
    -2.955065359477124e-02,
    1.796443723688307e-01,
    -1.39243221690590};
  
  x0 = x;
  if (x <= 0.0) return 1e308;
  else if ((x == 1.0) || (x == 2.0)) return 0.0;
  else if (x <= 7.0) {
    n = (int)(7-x);
    x0 = x+n;
  }
  x2 = 1.0/(x0*x0);
  xp = 2.0*M_PI;
  gl0 = a[9];
  for (k=8;k>=0;k--) {
    gl0 = gl0*x2 + a[k];
  }
  gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
  if (x <= 7.0) {
    for (k=1;k<=n;k++) {
      gl -= log(x0-1.0);
      x0 -= 1.0;
    }
  }
  return gl;
}
struct my_vg_params {double x;double mu; double sigma; double nu;};
double my_function_vg (double u, void *p){
  struct my_vg_params *params=(struct my_vg_params *)p;
  double x=(params->x);
  double mu = (params->mu);
  double sigma = (params->sigma);
  double nu=(params->nu);
  double f=(1.0/sqrt(M_2PI))*(sqrt(u)/sigma)*exp((-u*0.5/(sigma*sigma))
                                                   *(x-mu)*(x-mu))*exp((nu/2)*log(nu/2))*exp(-0.5*nu/u)*
                                                     exp(-lgamma(0.5*nu))
    *exp(-(0.5*nu+1)*log(u));
  return  f;
}
/*********************************************************************************/
// [[Rcpp::export]]
Rcpp::NumericVector pdf_vg(Rcpp::NumericVector  y,double mu, double sigma,double nu){
  double lower_limit=0.0;
  
  double abs_error = 1.0e-8;
  double rel_error = 1.0e-8;
  double result;
  double error;
  
  int nn=y.size();
  int kk=0;
  Rcpp::NumericVector res(nn);
  
  gsl_function G;
  
  for (kk=0;kk<nn;kk++){
    
    struct my_vg_params params={y[kk],mu,sigma,nu};
    G.function =&my_function_vg;
    G.params=&params;
    gsl_integration_qagiu(&G,lower_limit,abs_error,rel_error,1000,work_ptr, &result, &error);
    res[kk]=result;
    
  }
  return res;
}

