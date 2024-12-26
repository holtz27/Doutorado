// [[Rcpp::depends(RcppArmadillo, RcppGSL, RcppParallel)]]

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <RcppParallel.h>
#include <gsl/gsl_integration.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

vec svm_pw2pn(vec parvect){
  
  vec p(parvect.size());
  
  //(b0, b1, b2)
  p[0] = parvect[0];
  p[1] = (exp(parvect[1])-1)/(exp(parvect[1])+1);
  p[2] = parvect[2];
  // (mu, phi, sigma)
  p[3] = parvect[3];
  p[4] = (exp(parvect[4])-1)/(exp(parvect[4])+1);
  p[5] = exp(parvect[5]);
  // nu
  p[6] = exp(parvect[6]);
  
  return p;
}

double mlogLk_Rcpp(mat allprobs, mat egamma, rowvec foo, int n){
  double lscale=0;
  double sumfoo;
  
  for(int i=1;i<n;i++){		
    foo=foo*egamma%allprobs.row(i);
    sumfoo=sum(foo);
    lscale=lscale+log(sumfoo);
    foo=foo/sumfoo;
  }
  
  return lscale;
}

// Estrutura para os parâmetros da função
struct my_sl_params {
  double x, mu, sigma, nu;
};
// Função que será integrada
double my_function_sl(double u, void* p) {
  struct my_sl_params* params = (struct my_sl_params*)p;
  double x = params->x;
  double mu = params->mu;
  double sigma = params->sigma;
  double nu = params->nu;
  return (1.0 / sqrt(M_2PI)) * (sqrt(u) / sigma) * exp((-u * 0.5 / (sigma * sigma)) * (x - mu) * (x - mu)) * nu * exp((nu - 1) * log(u));
}
// Classe paralela para calcular os resultados
struct PDFWorker : public RcppParallel::Worker{
  // Entrada
  const NumericVector& y;
  double mu, sigma, nu;
  
  // Saída
  RcppParallel::RVector<double> res;
  
  // Construtor
  PDFWorker(const NumericVector& y_, double mu_, double sigma_, double nu_, NumericVector& res_)
    : y(y_), mu(mu_), sigma(sigma_), nu(nu_), res(res_) {}
  
  // Função que processa cada índice no intervalo
  void operator()(std::size_t begin, std::size_t end) override {
    gsl_integration_workspace* work_ptr = gsl_integration_workspace_alloc(1000);
    double lower_limit = 0.0;
    double upper_limit = 1.0;
    double abs_error = 1.0e-8;
    double rel_error = 1.0e-8;
    double result, error;
    
    gsl_function G;
    G.function = &my_function_sl;
    
    for (std::size_t i = begin; i < end; ++i) {
      my_sl_params params = { y[i], mu, sigma, nu };
      G.params = &params;
      
      gsl_integration_qag(&G, lower_limit, upper_limit, abs_error, rel_error, 1000,
                          GSL_INTEG_GAUSS21, work_ptr, &result, &error);
      
      res[i] = result;
    }
    
    gsl_integration_workspace_free(work_ptr);
  }
};
Rcpp::NumericVector pdf_s(const NumericVector& y, double mu, double sigma, double nu){
  Rcpp::NumericVector res(y.size());
  
  // Instancia a classe paralela
  PDFWorker worker(y, mu, sigma, nu, res);
  
  // Executa o loop em paralelo
  RcppParallel::parallelFor(0, y.size(), worker);
  
  return res;
}

mat OUTER(vec y, vec h, vec y0, vec beta, double nu){
  
  int n = y.size(), m = h.size();
  vec z = (y-beta[0]-beta[1]*y0-beta[2]*h[0]*h[0])/h[0];
  
  for(int t=1; t<m; t++){
    z = join_vert(z, (y-beta[0]-beta[1]*y0-beta[2]*h[t]*h[t])/h[t]);
  }
  
  NumericVector q(z.begin(), z.end());
  q = pdf_s(q, 0.0, 1.0, nu);
  
  mat M(q);
  M.reshape(n, m);
  M.each_row() /= h.t();
  
  return M;
}

// [[Rcpp::export]]
double svms_mllk2(vec parvect, vec y, double y0, int m, double gmax){
  int ny = y.size();
  vec p = svm_pw2pn(parvect);
  
  Function seq("seq");
  NumericVector b = seq(-gmax, gmax, 2*gmax/m); 
  vec b2 = b;
  vec bs = 0.5*(b2(span(1, b2.size()-1))+b2(span(0, b2.size()-2)));
  
  vec E = p[3]+p[4]*(bs-p[3]);
  vec sey = exp(0.5*bs);
  
  double intlen = b[2]-b[1];
  mat Gamma(m, m, fill::zeros);
  for(int i=0; i<m; i++){
    Gamma.row(i) = arma::normpdf(bs, E[i], p[5]).t() * intlen;
    //Gamma.row(i) /= sum(Gamma.row(i));
  }
  Gamma.each_col() /= sum(Gamma, 1);
  
  //vec xx = y;
  vec yy(ny);
  yy[0] = y0;
  yy(span(1, ny-1)) = y(span(0,ny-2));
  mat allprobs = OUTER(y, sey, yy, p.subvec(0, 2), p[6]);
  vec delta = arma::normpdf(bs, p[3], p[5]/sqrt(1-p[4]*p[4])) * intlen;  
  rowvec foo = delta.t() % allprobs.row(0);

  return -mlogLk_Rcpp(allprobs, Gamma, foo, ny);
}
