// [[Rcpp::depends(RcppArmadillo, RcppGSL, RcppParallel)]]

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct PdfVG2Worker : public Worker {
  const RVector<double> y;
  const double nu;
  RVector<double> vgdens;
  
  PdfVG2Worker(const NumericVector& y, double nu, NumericVector& vgdens)
    : y(y), nu(nu), vgdens(vgdens) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      double yi = y[i];
      if (nu < 100) {
        double c = std::sqrt(2.0 / M_PI) / gsl_sf_gamma(0.5 * nu);
        c *= std::pow(0.5 * nu, 0.5 * nu);
        c *= std::pow(nu, 0.25 * (1 - nu));
        
        double tol = 1e-5;
        
        if (std::abs(yi) < tol) {
          vgdens[i] = 0.5 * std::sqrt(nu / M_PI) * gsl_sf_gamma(0.5 * (nu - 1)) / gsl_sf_gamma(0.5 * nu);
        } else {
          double bessel_val = gsl_sf_bessel_Knu(0.5 * (nu - 1), std::abs(yi) * std::sqrt(nu));
          vgdens[i] = c * std::pow(std::abs(yi), 0.5 * (nu - 1)) * bessel_val;
        }
        
        if (std::isnan(vgdens[i])) {
          //Rcout << "warnings: y = " << yi << " nu = " << nu << std::endl;
          //warning("vgdens result was NaN");
          vgdens[i] = 0.0;
        }
        
        if (std::isinf(vgdens[i])) {
          //Rcout << "y = " << yi << " nu = " << nu << std::endl;
          stop("vgdens must be finite!");
        }
      } else {
        vgdens[i] = R::dnorm(yi, 0.0, 1.0, false);
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector pdf_vg(NumericVector y, double nu) {
  NumericVector vgdens(y.size());
  PdfVG2Worker worker(y, nu, vgdens);
  parallelFor(0, y.size(), worker);
  return vgdens;
}
