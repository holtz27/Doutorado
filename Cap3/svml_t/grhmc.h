#ifndef GRHMC_H_INCLUDED
#define GRHMC_H_INCLUDED

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

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

class GRHMC {

private:
  static SvmStats make_stats(const vec& y, const vec& h, const vec& l, 
                                    const int T);
  
  static void chol_tridiagonal(const vec& b, const vec& cc, vec& l_diag, 
                                     vec& m_subdiag, const int T);
  
  static void apply_chol_tri(const vec& ldiag, const vec& msubdiag, const vec& x,
                                   vec& vP, const int T);
  static double sample_r(const double m, const double c, const int d);
  
  static vec sample_pI(const double m, const double c, const int d);
  
public:
  
  GRHMC() = default;
  
  // #############################################################################
  // ############################ theta ##########################################
  // #############################################################################
  static double lpost_theta(const vec& theta, const vec& y, const vec& h, 
                            const vec& l, const int T);
  static vec grad_theta(const vec& theta, const vec& y, const vec& h, const vec& l,
                        const int T, const mat& Lchol, const mat& Iparb1, 
                        const mat& Iparmu, const mat& Iparphi, const mat& Ipardelta, 
                        const mat& Ipartau, const mat& Iparnu);
  
  static vec nuHtheta(const mat& Lchol, const mat& Iparb1, const mat& Iparmu, 
                      const mat& Iparphi, const mat& Ipardelta, const mat& Ipartau, 
                      const mat& Iparnu, const vec& p, const double m, 
                      const double c); 
  
  static vec gradHmomtheta(const mat& Lchol, const vec& p, const double m, 
                           const double c);
  
  static void I_theta(const vec& theta, const SvmStats& st, 
                      mat& I, mat& Iparb1, mat& Iparmu, mat& Iparphi, 
                      mat& Ipardelta, mat& Ipartau, mat& Iparnu);
  
  static void I_theta_ref(const vec& theta, const vec& y, const vec& h, const vec& l, 
                          const int T, mat& I, mat& Iparb1, mat& Iparmu, mat& Iparphi, 
                          mat& Ipardelta, mat& Ipartau, mat& Iparnu);
  
  static vec updatetheta(const vec& theta_cur, const vec& y, const vec& h_cur, 
                         const vec& l_cur, const int T, const int nfix,
                         const int nlpfrog, const double epsilon, const double m, 
                         const double c, int& acc);
  
  // #############################################################################
  // ############################ h ##############################################
  // #############################################################################
  static double lpost_h(const vec& h, const vec& theta, const vec& y, 
                        const vec& l, const int T);
  
  static vec grad_h(const vec& h, const vec& theta, const vec& y, const vec& l, 
                    const int T, const vec& ldiag, const vec& msubdiag);
  
  static vec nuHh(const vec& h, const vec& theta, const vec& l, const int T,
                  const vec& ldiag, const vec& msubdiag, const vec& p, const double m, 
                  const double c);
  
  static vec gradHmomh(const vec& ldiag, const vec& msubdiag, const vec& p, 
                       const int T, const double m, 
                       const double c);
  
  static bool I_h(const vec& h, const vec& theta, const vec& l, const int T,
                  vec& ldiag_I, vec& msubdiag_I, vec& ldiag, vec& msubdiag);
  
  static vec updateh(const vec& h_cur, const vec& theta_cur, const vec& y, 
                     const vec& l_cur, const int T, const int nfix, 
                     const int nlpfrog, const double epsilon, const double m, 
                     const double c, int& acc);
 
  
  // #############################################################################
  // ############################ l ##############################################
  // #############################################################################
  static double lpost_l(const vec& l_star, const vec& h, const vec& theta, 
                        const vec& y, const int T);
  
  static vec grad_l(const vec& l_star, const vec& h, const vec& theta, 
                    const vec& y, 
                    const int T);
  
  static vec nuHl(const vec& l_star, const vec& theta, const vec& p, 
                  const int T, 
                  const double m, const double c);
  
  static vec gradHmoml(const vec& l_star, const vec& theta, const vec& p, 
                       const int T, 
                       const double m, const double c);
  
  static vec I_l(const vec& l_star, const vec& theta, const int T);
  
  static vec updatel(const vec& l_star_cur, const vec& h_cur, const vec& theta_cur, 
                     const vec& y, const int T, const int nfix, 
                     const int nlpfrog, const double epsilon,
                     const double m, const double c, int& acc);
 
  
};

#endif // GRHMC_H_INCLUDED