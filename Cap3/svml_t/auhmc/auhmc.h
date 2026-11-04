#ifndef AUHMC_H_INCLUDED
#define AUHMC_H_INCLUDED

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Estatisticas que dependem apenas de (y, h, lambda). Somas em t = 1..T-1
// (indices 0..T-2), mais os termos t = T.
struct SvmStats {
  int    T;
  double h1;                 // h_1
  // termos t = T
  double s2_T;               // s_T^2
  double ylag_T;             // y_{T-1}
  double lam_T;              // lambda_T
  double lamexp_T;           // lambda_T * exp(h_T)
  // somas em t = 1..T-1
  double S_s2, S_s2y, S_s2y2;          // sum s^2, s^2 y, s^2 y^2
  double S_lam, S_lamy, S_lamexp;      // sum lam, lam y, lam e^h
  double S_s, S_sy, S_sexp;            // sum s, s y, s e^h
  double S_sh, S_syh, S_sexph;         // sum s h, s y h, s e^h h
  double S_h, S_h2;                    // sum h, sum h^2
  // soma em t = 1..T
  double sum_logl_minus_l;             // sum (log lambda - lambda)
};

class AUHMC {

private:
  // Método privado auxiliar para preenchimento de estatísticas estáticas
  static SvmStats make_stats(const arma::vec& y, const arma::vec& h, 
                             const arma::vec& l, const int T);
  static inline vec solve_L(const vec& ldiag, const vec& msubdiag, const vec& b, 
                            const int T);
  static inline vec solve_Lt(const vec& ldiag, const vec& msubdiag, const vec& b, 
                             const int T);
  static inline vec solve_M(const vec& ldiag, const vec& msubdiag, const vec& b, 
                            const int T);
  static void chol_tridiagonal(const vec& b, const vec& cc, vec& l_diag, 
                               vec& m_subdiag, const int T);
  
  static void apply_chol_tri(const vec& ldiag, const vec& msubdiag, const vec& x,
                             vec& vP, const int T);
  
public:
  AUHMC() = default;
  
  static double jitter_epsilon(const double eps_base, const double Umin, 
                               const double Umax);
  
  // #############################################################################
  // ############################ theta ##########################################
  // #############################################################################
  static double lpost_theta(const vec& theta, const vec& y, const vec& h, 
                                   const vec& l, const int T);
  
  static vec grad_theta(const vec& theta, const vec& y, const vec& h, const vec& l,
                               const int T);
  
  static void I_theta(const vec& theta, const SvmStats& st, 
                             mat& I);
  
  static vec updatetheta(const vec& theta_cur, const vec& y, const vec& h_cur, 
                         const vec& l_cur, const int T, const int nlpfrog, 
                         const double epsilon, int& acc, 
                         int& total_j_iters, int& non_conv_count);
  
  // #############################################################################
  // ############################ h ##############################################
  // #############################################################################
  // 1. Densidade posterior condicional log p(h | theta, y, lambda) e seu gradiente
  static double lpost_h(const vec& h, const vec& theta, const vec& y, const vec& l, 
                        const int T);
  
  static vec grad_h(const vec& h, const vec& theta, const vec& y, const vec& l, 
                    const int T);
  
  // 2. Matriz de Informacao de Fisher tridiagonal I(h_T)
  static bool I_h(const vec& h, const vec& theta, const vec& l, const int T,
                  vec& ldiag_I, vec& msubdiag_I, vec& ldiag, vec& msubdiag);
  
  static vec updateh(const vec& theta_cur, const vec& y, const vec& h_cur, 
                            const vec& l_cur, const int T, const int nlpfrog, 
                            const double epsilon, int& acc, 
                            int& total_j_iters, int& non_conv_count);
  
  // #############################################################################
  // ############################ l ##############################################
  // #############################################################################
  static double lpost_l(const vec& l_star, const vec& h, const vec& theta, 
                               const vec& y, const int T);
  
  static vec grad_l(const vec& l_star, const vec& h, const vec& theta, const vec& y, 
                           const int T);
  
  static vec I_l(const vec& l_star, const vec& theta, const int T);
  
  static vec updatel(const vec& l_star_cur, const vec& h_cur, const vec& theta_cur,
                            const vec& y, const int T, const int nlpfrog,
                            const double epsilon, int& acc,
                            int& total_j_iters, int& non_conv_count);

};

#endif // AUHMC_H_INCLUDED
