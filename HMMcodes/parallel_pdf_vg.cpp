// [[Rcpp::depends(RcppGSL, RcppParallel)]]
#include <Rcpp.h>
#include <RcppGSL.h>
#include <RcppParallel.h>
#include <gsl/gsl_integration.h>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>


struct my_vg_params {
  double x, mu, sigma, nu;
};
double my_function_vg (double u, void *p){
  struct my_vg_params *params=(struct my_vg_params *)p;
  double x=(params->x);
  double mu = (params->mu);
  double sigma = (params->sigma);
  double nu=(params->nu);
  double f=(1.0/sqrt(M_2PI))*(sqrt(u)/sigma)*exp((-u*0.5/(sigma*sigma))*
            (x-mu)*(x-mu))*exp((nu/2)*log(nu/2))*exp(-0.5*nu/u)*
            exp(-gsl_sf_lngamma(0.5 * nu))* 
            exp(-(0.5*nu+1)*log(u));
  return  f;
}
struct PDFWorker : public RcppParallel::Worker {
  // Entrada
  const Rcpp::NumericVector& y;
  double mu, sigma, nu;
  
  // Saída
  RcppParallel::RVector<double> res;
  
  // Construtor
  PDFWorker(const Rcpp::NumericVector& y_, double mu_, double sigma_, double nu_, Rcpp::NumericVector& res_)
    : y(y_), mu(mu_), sigma(sigma_), nu(nu_), res(res_) {}
  
  // Função que processa cada índice no intervalo
  void operator()(std::size_t begin, std::size_t end){ //override {
    gsl_integration_workspace* work_ptr = gsl_integration_workspace_alloc(1000);
    double lower_limit = 0.0;
    double abs_error = 1.0e-4;
    double rel_error = 1.0e-4;
    double result, error;
    
    gsl_function G;
    G.function = &my_function_vg;
    
    for(std::size_t i = begin; i < end; ++i){
      my_vg_params params = { y[i], mu, sigma, nu };
      G.params = &params;
      
      gsl_integration_qagiu(&G,lower_limit,abs_error,rel_error,1000,work_ptr, &result, &error);
      res[i] = result;
    }
    gsl_integration_workspace_free(work_ptr);
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector pdf_vg(const Rcpp::NumericVector& y, double mu, double sigma, double nu){
  
  Rcpp::NumericVector res(y.length());
  
  PDFWorker worker(y, mu, sigma, nu, res);
  RcppParallel::parallelFor(0, y.size(), worker);
  
  return res;
}