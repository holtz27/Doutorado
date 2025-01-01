// [[Rcpp::depends(RcppArmadillo, RcppGSL, RcppParallel)]]

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
//-----------------------------------------------------------------------------
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

struct my_sl_params {double x; double mu; double sigma; double nu;};
double my_function_sl(double u, void* p) {
  struct my_sl_params* params = (struct my_sl_params*)p;
  double x = (params->x);
  double mu = (params->mu);
  double sigma = (params->sigma);
  double nu = (params->nu);
  double f = (1.0 / sqrt(M_2PI)) * (sqrt(u) / sigma) * 
    exp((-u * 0.5 / (sigma * sigma)) * (x - mu) * (x - mu)) * 
    nu * exp((nu - 1) * log(u));
  return f;
}
double pdf_s2(double y, double mu, double sigma, double nu) {
  gsl_integration_workspace* work_ptr = gsl_integration_workspace_alloc(1000);
  double lower_limit = 0.0;
  double upper_limit = 1.0;
  double abs_error = 1.0e-8;
  double rel_error = 1.0e-8;
  double result;
  double error;
  
  gsl_function G;
  struct my_sl_params params = {y, mu, sigma, nu};
  G.function = &my_function_sl;
  G.params = &params;
  
  gsl_integration_qag(&G, lower_limit, upper_limit, abs_error, rel_error, 
                      1000, GSL_INTEG_GAUSS21, work_ptr, &result, &error);
  gsl_integration_workspace_free(work_ptr);
  
  return result;
}
struct Outer2Worker : public Worker {
  const vec& y;
  const vec& sey;
  const vec& beta;
  double nu;
  const vec& ylag;
  mat& M;
  
  Outer2Worker(const vec& y, const vec& sey, const vec& beta, double nu, const vec& ylag, mat& M)
    : y(y), sey(sey), beta(beta), nu(nu), ylag(ylag), M(M) {}
  
  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t index = begin; index < end; index++) {
      std::size_t i = index / sey.n_elem;
      std::size_t j = index % sey.n_elem;
      double z = (y[i] - beta[0] - beta[1] * ylag[i] - beta[2] * sey[j] * sey[j]) / sey[j];
      M(i, j) = pdf_s2(z, 0.0, 1.0, nu) / sey[j];
    }
  }
};


double logprior(vec parvect){
  
  double logp = log_normpdf(parvect[0], 0.0, 10.0);
  logp += log_normpdf(parvect[1], 0.5, 10.0);
  logp += log_normpdf(parvect[2], 0.0, 10.0);
  logp += log_normpdf(parvect[3], 0.0, 10.0);
  logp += log_normpdf(parvect[4], 4.5, 10.0);
  logp += log_normpdf(parvect[5], -1.5, 10.0);
  logp += log_normpdf(parvect[6], 0.0, 10.0);
  
  return -logp;
}
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
  }
  Gamma.each_col() /= sum(Gamma, 1);
  
  vec ylag(ny);
  ylag[0] = y0;
  ylag(span(1, ny-1)) = y(span(0, ny-2));
  
  mat allprobs(ny, m, fill::zeros);
  Outer2Worker worker(y, sey, p.subvec(0, 2), p[6], ylag, allprobs);
  parallelFor(0, ny * m, worker);
  
  
  //mat allprobs = OUTER(y, sey, yy, p.subvec(0, 2), p[6]);
  
  vec delta = arma::normpdf(bs, p[3], p[5]/sqrt(1-p[4]*p[4])) * intlen;  
  rowvec foo = delta.t() % allprobs.row(0);
  
  double result = -mlogLk_Rcpp(allprobs, Gamma, foo, ny); 
  
  return result;
}
double posterior(vec parvect, vec y, double y0, int m, double gmax){
  return svms_mllk2(parvect, y, y0, m, gmax) + logprior(parvect);  
}

List svsl_map_Rcpp(const vec parvect_init, vec y, double y0, int m, double gmax){
  
  Environment stats("package:stats"); 
  Function optim = stats["optim"];
  
  List opt_results = optim(_["par"]    = parvect_init,
                           _["fn"]     = InternalFunction(&posterior),
                           _["method"] = "BFGS",
                           _["hessian"] = true,
                           _["y"] = y,
                           _["y0"] = y0,
                           _["m"] = m,
                           _["gmax"] = gmax);
  return opt_results;
}

double MVDnorm(vec x, vec mu, mat L, double log_det){
  int dim = mu.n_elem;
  
  // Verificar compatibilidade das dimensões
  //if (sigma.n_rows != dim || sigma.n_cols != dim) {
  //  stop("A matriz de covariância deve ser quadrada e corresponder ao tamanho do vetor de médias.");
  //}
  //if (x.n_elem != dim) {
  //  stop("O vetor de entrada deve ter a mesma dimensão que o vetor de médias.");
  //}
  
  // Decomposição de Cholesky
  //arma::mat L = arma::chol(sigma, "lower");
  
  // Calcular o determinante usando a decomposição de Cholesky
  //double log_det_sigma = 2 * sum(log(L.diag()));
  
  // Calcular o vetor z = L^(-1) * (x - mu)
  arma::vec diff = x - mu;
  arma::vec z = arma::solve(L, diff, arma::solve_opts::fast);
  
  // Calcular a densidade
  double quad_form = arma::dot(z, z);
  double log_density = -0.5 * (dim * std::log(2 * M_PI) + log_det + quad_form);
  
  return log_density;
}
mat MVRnorm(int n, vec mu, mat sigma){
  int dim = mu.n_rows;
  
  // Verificar compatibilidade das dimensões
  if (sigma.n_rows != dim || sigma.n_cols != dim) {
    stop("A matriz de covariância deve ser quadrada e corresponder ao tamanho do vetor de médias.");
  }
  
  // Configurar gerador
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default; // Usa o gerador padrão
  gsl_rng *r = gsl_rng_alloc(T);
  
  // Converter mu para gsl_vector
  gsl_vector *gsl_mu = gsl_vector_alloc(dim);
  for (int i = 0; i < dim; ++i) {
    gsl_vector_set(gsl_mu, i, mu[i]);
  }
  
  // Converter sigma para gsl_matrix
  gsl_matrix *gsl_sigma = gsl_matrix_alloc(dim, dim);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      gsl_matrix_set(gsl_sigma, i, j, sigma(i, j));
    }
  }
  if (gsl_linalg_cholesky_decomp(gsl_sigma) != GSL_SUCCESS) {
    gsl_matrix_free(gsl_sigma);
    gsl_rng_free(r);
    stop("Erro ao realizar a decomposição de Cholesky.");
  }
  
  // Matriz de resultados
  mat result(dim, n);
  gsl_vector *gsl_sample = gsl_vector_alloc(dim);
  
  // Gerar amostras
  for (int i = 0; i < n; ++i) {
    if (gsl_ran_multivariate_gaussian(r, gsl_mu, gsl_sigma, gsl_sample) != GSL_SUCCESS){
      gsl_rng_free(r);
      gsl_vector_free(gsl_mu);
      gsl_vector_free(gsl_sample);
      gsl_matrix_free(gsl_sigma);
      stop("Erro ao gerar amostra.");
    }
    
    for (int j = 0; j < dim; ++j) {
      result(j, i) = gsl_vector_get(gsl_sample, j);
    }
  }
  
  // Liberar memória
  gsl_rng_free(r);
  gsl_vector_free(gsl_mu);
  gsl_vector_free(gsl_sample);
  gsl_matrix_free(gsl_sigma);
  
  return result;
}

// [[Rcpp::export]]
List Weigth_cpp2(const vec parvect_init, vec y, double y0, int m, double gmax, int n){
  
  List opt_results = svsl_map_Rcpp(parvect_init, y, y0, m, gmax);
  
  double k = opt_results["value"];
  vec map = opt_results["par"];
  mat H1 = opt_results["hessian"];
  H1 = H1.i(); 
  mat L = arma::chol(H1, "lower");
  double log_det = 2 * sum(log(L.diag()));
  
  mat X2 = MVRnorm(n, map, H1);
  
  vec weigths(n);
  for(int i=0; i < n; i++){
    weigths(i) = exp(k -posterior(X2.col(i), y, y0, m, gmax)
                       -MVDnorm(X2.col(i), map, L, log_det));
  }
  if(sum(weigths) == 0) stop("Instabilidade numérica");
  weigths /= sum(weigths);
  
  //return weigths;
  return List::create(Named("weigths") = weigths , _["X"] = X2);
}