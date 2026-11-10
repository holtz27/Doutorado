// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <chrono>
#include "rmlmc.h"

using namespace Rcpp;
using namespace arma;

// Função utilitária para exibir a barra de progresso customizada
inline void print_progress(int current, int total, int bar_width = 30) {
  int step = total / 100;
  if (step == 0) step = 1;
  
  if (current % step == 0 || current == total) {
    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(bar_width * progress);
    
    Rcpp::Rcout << "\r[";
    for (int i = 0; i < bar_width; ++i) {
      if (i < pos) Rcpp::Rcout << "=";
      else if (i == pos) Rcpp::Rcout << ">";
      else Rcpp::Rcout << " ";
    }
    Rcpp::Rcout << "] " << int(progress * 100.0) << "% (" 
                << current << "/" << total << ")";
    Rcpp::Rcout.flush();
  }
}

// [[Rcpp::export]]
List rmlmc_svml(const vec& y,
                int n_draws, int burnin, int thin, vec nlpfrog, 
                vec epsilon, const vec alpha, bool show_progress = true)
{
  int T = y.n_elem - 1; // y = (y_0, y_1, ..., y_T)
  
  // Matrizes para armazenar as amostras da posterior (draws)
  mat draws_theta(8, n_draws);
  mat draws_h(T, n_draws);
  mat draws_l_star(T, n_draws);
  
  // Contadores de aceitação do RMHMC
  int acc_theta = 0, acc_h = 0, acc_l = 0;
  int div_theta = 0, div_h = 0, div_l = 0;
  
  // 1. Inicialização automatizada de theta_cur
  vec theta_cur(8);
  theta_cur[0] = 0.1;                        // b0
  theta_cur[1] = std::atanh(0.0);            // b1* = atanh(0) = 0
  theta_cur[2] = 0.0;                        // b2
  theta_cur[3] = 0.1;                        // mu
  theta_cur[4] = std::atanh(0.98);           // phi* = atanh(0.98)
  theta_cur[5] = -0.06;                      // delta
  theta_cur[6] = std::log(0.14);             // tau* = log(0.14)
  theta_cur[7] = -1.3;                       // nu*
  
  // 2. Inicialização automatizada de l_star_cur (log-lambda)
  vec l_star_cur = randg<vec>(T, arma::distr_param(5.0, 1.0/5.0));
  l_star_cur = log(l_star_cur);
  
  // 3. Inicialização automatizada de h_cur com base nos dados y, theta e l_star
  const vec y_obs  = y.subvec(1, T);      
  const vec y_lag1 = y.subvec(0, T - 1);  
  vec mu_t = theta_cur[0] + std::tanh(theta_cur[1]) * y_lag1;
  vec res2 = square(y_obs - mu_t) + 1e-4;
  vec h_cur = l_star_cur + log(res2) + 1.2704;
  
  // Total de iterações do laço Gibbs (burnin + n_draws * thin)
  int total_iters = burnin + n_draws * thin;
  
  // Medição do tempo de execução do loop Gibbs
  auto start_time = std::chrono::high_resolution_clock::now();
  
  int saved_draws = 0; // Ponteiro para a coluna salva
  
  for (int iter = 0; iter < total_iters; ++iter) 
  {
    if (show_progress) {
      print_progress(iter + 1, total_iters, 30);
    }
    
    // 1. Amostragem de theta | h, l*, y
    theta_cur = RMLMC::updatetheta(theta_cur, y, h_cur, exp(l_star_cur), T, 
                                   alpha[0], 
                                   nlpfrog[0], 
                                   epsilon[0], 
                                   acc_theta, div_theta);
    
    // 2. Amostragem de h | theta, l*, y
    h_cur = RMLMC::updateh(h_cur, exp(l_star_cur), theta_cur, y, T, 
                           alpha[1],
                           nlpfrog[1], 
                           epsilon[1], 
                           acc_h, div_h);
    
    // 3. Amostragem de l* | theta, h, y (no espaço não-restrito R^T)
    l_star_cur = RMLMC::updatel(l_star_cur, h_cur, theta_cur, y, T, 
                                alpha[2],
                                nlpfrog[2], 
                                epsilon[2], 
                                acc_l, div_l);
    
    // Salva a amostra se estiver pós burn-in e no passo múltiplo de 'thin'
    if (iter >= burnin && ((iter - burnin) % thin == 0) && saved_draws < n_draws) {
      draws_theta.col(saved_draws)  = theta_cur;
      draws_h.col(saved_draws)      = h_cur;
      draws_l_star.col(saved_draws) = l_star_cur;
      saved_draws++;
    }
  }
  
  if (show_progress) {
    Rcpp::Rcout << std::endl; // Quebra linha ao finalizar a barra
  }
  
  auto end_time = std::chrono::high_resolution_clock::now();
  double elapsed_sec = std::chrono::duration<double>(end_time - start_time).count();
  
  // Cálculo das taxas de aceitação finais
  double acc_theta_rate = static_cast<double>(acc_theta) / total_iters;
  double acc_h_rate     = static_cast<double>(acc_h) / total_iters;
  double acc_l_rate     = static_cast<double>(acc_l) / total_iters;
  
  double div_theta_rate = static_cast<double>(div_theta);  // / total_iters;
  double div_h_rate     = static_cast<double>(div_h) / total_iters;
  double div_l_rate     = static_cast<double>(div_l) / total_iters;
  
  double elapsed_min    = elapsed_sec / 60.0;
  
  // Impressão das estatísticas diretamente no console do R
  if (show_progress) {
    Rcpp::Rcout << "  theta      h      l   time(min)" << std::endl;
    Rcpp::Rcout << "   " << std::fixed << std::setprecision(2)
                << acc_theta_rate << "   " 
                << acc_h_rate << "   " 
                << acc_l_rate << "   " 
                << elapsed_min << std::endl;
  }
  
  // 4. Mapeamento de theta* para a escala original theta
  draws_theta.row(1) = tanh(draws_theta.row(1));                    // beta1
  draws_theta.row(4) = tanh(draws_theta.row(4));                    // phi
  draws_theta.row(6) = exp(draws_theta.row(6));                     // tau
  draws_theta.row(7) = 2.0 + 38.0 / (1.0 + exp(-draws_theta.row(7))); // nu
  
  return List::create(
    Named("theta")        = draws_theta,
    Named("h")            = draws_h,
    Named("l")            = exp(draws_l_star), 
    Named("acc_theta")    = acc_theta_rate,
    Named("acc_h")        = acc_h_rate,
    Named("acc_l")        = acc_l_rate,
    Named("div_theta")    = div_theta_rate,
    Named("div_h")        = div_h_rate,
    Named("div_l")        = div_l_rate,
    Named("elapsed_sec")  = elapsed_sec
  );
}