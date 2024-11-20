// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

class VG: public Func{
private:
  double x;
  double nu;
public:
  VG(double x_, double nu_) : x(x_), nu(nu_) {}
  
  double operator()(const double& u) const
  {
    return sqrt(0.5*u/M_PI) * exp(-0.5*u*x*x-0.5*nu/u + 0.5*nu*log(0.5*nu)-lgamma(0.5*nu)-(1+0.5*nu)*log(u));
  }
};

// [[Rcpp::export]]
Rcpp::List pdf_vg(Rcpp::NumericVector y, double mu, double sigma, double nu){
  
  double lower = 0, upper = R_PosInf, err_est;
  int err_code;
  
  int n = y.length();
  Rcpp::NumericVector res(n);
  Rcpp::NumericVector err_vec(n);
  
  for(int i = 0; i < n; i++){
    VG f(y[i], nu);
    res[i] = integrate(f, lower, upper, err_est, err_code);
    err_vec[i] = err_est;
    
    if(err_code != 0) {
      Rcpp::Rcout << "Erro na integração para y[" << i << "]: código " << err_code << "\n";
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("result") = res,      
    Rcpp::Named("error_estimate") = err_vec 
  );
}