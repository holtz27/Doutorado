// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppNumerical.h>
#include <RcppParallel.h>
using namespace Numer;
using namespace RcppParallel;

// Classe VG para o integrando
class VG : public Func {
private:
  double x;
  double nu;
public:
  VG(double x_, double nu_) : x(x_), nu(nu_) {}
  
  double operator()(const double& u) const {
    return sqrt(0.5 * u / M_PI) * exp(
        -0.5 * u * x * x 
    - 0.5 * nu / u 
    + 0.5 * nu * log(0.5 * nu) 
    - lgamma(0.5 * nu) 
    - (1 + 0.5 * nu) * log(u)
    );
  }
};

// Classe paralela para calcular as integrais
struct VGWorker : public Worker {
  // Input
  const Rcpp::NumericVector& y;
  const double nu;
  const double lower;
  const double upper;
  
  // Output
  Rcpp::NumericVector& result;
  Rcpp::NumericVector& error_estimates;
  
  // Constructor
  VGWorker(const Rcpp::NumericVector& y_, double nu_, 
           double lower_, double upper_,
           Rcpp::NumericVector& result_, Rcpp::NumericVector& error_estimates_)
    : y(y_), nu(nu_), lower(lower_), upper(upper_), 
      result(result_), error_estimates(error_estimates_) {}
  
  // Operador para processamento paralelo
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      double err_est = 0;
      int err_code = 0;
      VG f(y[i], nu);
      
      result[i] = integrate(f, lower, upper, err_est, err_code);
      error_estimates[i] = err_est;
    }
  }
};

// [[Rcpp::export]]
Rcpp::List pdf_vg(Rcpp::NumericVector y, double mu, double sigma, double nu) {
  const double lower = 0, upper = R_PosInf;
  int n = y.length();
  
  Rcpp::NumericVector res(n);
  Rcpp::NumericVector err_vec(n);
  
  // Criar a instância da classe paralela
  VGWorker worker(y, nu, lower, upper, res, err_vec);
  
  // Processar em paralelo
  parallelFor(0, n, worker);
  
  return Rcpp::List::create(
    Rcpp::Named("result") = res,
    Rcpp::Named("error_estimate") = err_vec
  );
}
