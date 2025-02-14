// [[Rcpp::depends(RcppArmadillo, RcppGSL)]]

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pdf_vg2(NumericVector y, double nu){
  
  int n = y.size();
  Rcpp::NumericVector vgdens(n);
  
  for(int i = 0; i < n; i++){
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
        Rcout << "warnings: y = " << yi << " nu = " << nu << std::endl;
        warning("vgdens result was NaN");
        vgdens[i] = 0.0;
      }
      
      if (std::isinf(vgdens[i])) {
        Rcout << "y = " << yi << " nu = " << nu << std::endl;
        stop("vgdens must be finite!");
      }
    } else {
      vgdens[i] = R::dnorm(yi, 0.0, 1.0, false);
    }
  }
  
  return vgdens;
}
