// [[Rcpp::depends(RcppArmadillo, RcppGSL, RcppParallel)]]

#include <RcppArmadillo.h>
#include <gsl/gsl_randist.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct PdfVG2Worker : public Worker{
  const RVector<double> y;
  RVector<double> ndens;
  
  PdfVG2Worker(const NumericVector& y, NumericVector& ndens)
    : y(y), ndens(ndens) {}
  
  void operator()(std::size_t begin, std::size_t end){
    for(std::size_t i = begin; i < end; i++){
      double yi = y[i];
      ndens[i]=gsl_ran_gaussian_pdf(yi, 1.0);
    }
  }
};

// [[Rcpp::export]]
NumericVector pdf_n(NumericVector y) {
  NumericVector ndens(y.size());
  PdfVG2Worker worker(y, ndens);
  parallelFor(0, y.size(), worker);
  return ndens;
}

