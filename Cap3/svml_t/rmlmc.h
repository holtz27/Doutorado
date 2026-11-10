#ifndef RMLMC_H_INCLUDED
#define RMLMC_H_INCLUDED

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Estatisticas que dependem apenas de (y, h, lambda). Somas em t = 1..T-1
// (indices 0..T-2), mais os termos t = T.
struct SvmStats {
  int    T;
  double h1;
  
  // termos t = T
  double s2_T, ylag_T, lam_T, lamexp_T;
  
  // somas em t = 1..T-1
  double S_s2, S_s2y, S_s2y2;
  double S_lam, S_lamy, S_lamexp;
  double S_s, S_sy, S_sexp;
  double S_sh, S_syh, S_sexph;
  double S_h, S_h2;
  
  double sum_logl_minus_l;
  
  // vetores completos (novos), todos de comprimento T
  vec s;      // s_t = sqrt(lambda_t) * exp(-h_t/2),  t = 1..T
  vec yt;     // y_1..y_T
  vec ylag;   // y_0..y_{T-1}
  vec eh;     // exp(h_t)
  vec h;      // h_1..h_T
};

class RMLMC {

private:
  
  static SvmStats make_stats(const vec& y, const vec& h, const vec& l, const int T);
    
public:
  RMLMC() = default;
  
  // #############################################################################
  // ############################ theta ##########################################
  // #############################################################################
  static double lpost_theta(const vec& theta, const vec& y, const vec& h, 
                     const vec& l, const int T);
  
  static vec grad_theta(const vec& theta, const vec& y, const vec& h, const vec& l,
                    const int T);
  
  static void H_theta(const vec& theta, const SvmStats& st, mat& H, const int T);
  
  static void sqrtG_theta(const vec& theta, const vec& y, const vec& h, 
                                 const vec& l,
                                 const int T, const double alpha, mat& sqrtG);
  
  static double Etheta(const vec& theta, const vec& y, const vec& h, 
                       const vec& l, const int T, const vec& v, const double alpha);
    
  static vec updatetheta(const vec& theta_cur, const vec& y, const vec& h_cur, 
                         const vec& l_cur, const int T, const double alpha,
                         const int nlpfrog, const double epsilon, int& acc, 
                         int& div);
  
  // #############################################################################
  // ############################ h ##############################################
  // #############################################################################
  static double lpost_h(const vec& h, const vec& theta, const vec& y, const vec& l, 
                        const int T);
    
  static  vec grad_h(const vec& h, const vec& theta, const vec& y, const vec& l, 
                      const int T);
    
  static  void H_h(const vec& h_star, const vec& l_cur, const vec& theta, 
                    const vec& y, const int T, vec& diagH, vec& subH);
    
  static  vec mult_H_tridiag(const vec& diagH, const vec& subH, const vec& x, 
                               const int T);
    
  static vec sample_v_h(const vec& h_cur, const vec& l_cur, const vec& theta_cur, 
                          const vec& y, const int T, const double alpha);
    
  static double Eh(const vec& h_star, const vec& l_cur, const vec& theta, const vec& y, 
                     const int T, const vec& v, const double alpha);
    
  static vec updateh(const vec& h_cur, const vec& l_cur, const vec& theta_cur,
                       const vec& y, const int T, const double alpha,
                       const int nlpfrog, const double epsilon,
                       int& acc, int& div);
  
  // #############################################################################
  // ############################ l ##############################################
  // #############################################################################
  static double lpost_l(const vec& l_star, const vec& h, const vec& theta, 
                               const vec& y, const int T);
  
  static vec grad_l(const vec& l_star, const vec& h, const vec& theta, const vec& y, 
                           const int T);
  
  static vec Hdiag_l(const vec& l_star, const vec& h, const vec& theta,
                            const vec& y, const int T);
  
  static vec sample_v(const vec& l_star, const vec& h, const vec& theta, 
                             const vec& y, const int T, const double alpha);
  
  static double El(const vec& l_star, const vec& h, const vec& theta, const vec& y, 
                          const int T, const vec& v, const double alpha);
  
  static vec updatel(const vec& l_cur, const vec& h_cur, const vec& theta_cur,
                            const vec& y, const int T, const double alpha,
                            const int nlpfrog, const double epsilon,
                            int& acc, int& div);

  
};

#endif // RMLMC_H_INCLUDED