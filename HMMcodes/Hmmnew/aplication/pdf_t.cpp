// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace RcppParallel;

struct PdfS2Worker : public Worker {
  const RVector<double> y;
  const double df;
  RVector<double> result;
  
  PdfS2Worker(const NumericVector& y, double df, NumericVector& result)
    : y(y), df(df), result(result) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    
    for (std::size_t i = begin; i < end; i++){
      
      result[i] = R::dt(y[i], df, 0);
      
    }
  }
};


// [[Rcpp::export]]
NumericVector pdf_t(NumericVector y, double df) {
  NumericVector result(y.size());
  PdfS2Worker worker(y, df, result);
  parallelFor(0, y.size(), worker);
  return result;
}
