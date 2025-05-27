// [[Rcpp::depends(RcppArmadillo, RcppGSL, RcppParallel)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <gsl/gsl_sf_gamma.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace RcppParallel;

struct PdfS2Worker : public Worker {
  const RVector<double> y;
  const double nu;
  RVector<double> result;
  
  PdfS2Worker(const NumericVector& y, double nu, NumericVector& result)
    : y(y), nu(nu), result(result) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    double a = nu / sqrt(2 * M_PI);
    double s = nu + 0.5;
    double tol = 1e-8;
    double z, x;
    
    for (std::size_t i = begin; i < end; i++) {
      
      if(nu < 100){
        x = 0.5 * y[i] * y[i];
        if (fabs(x) < tol) {
          z = a / s;
        } else {
          z = a*pow(1.0 / x, s)*gsl_sf_gamma_inc_P(s, x);
          z *= gsl_sf_gamma(s);
        }
        
        if(std::isnan(z)){
          z = 0.0;
        }
        
        if(std::isinf(z)){
          stop("sdens must be finite");
        }
        
      }else{
        z = R::dnorm(y[i], 0.0, 1.0, false);
      }
      
      result[i] = z;
    }
  }
};

// [[Rcpp::export]]
NumericVector pdf_s(NumericVector y, double nu) {
  NumericVector result(y.size());
  PdfS2Worker worker(y, nu, result);
  parallelFor(0, y.size(), worker);
  return result;
}
