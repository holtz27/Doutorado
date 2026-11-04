#include "auhmc.h"
#include <cmath>

// Parameters priori
const double mu_b0 = 0.0, s_b0 = 3.2;
const double a_b1 = 5.0, b_b1 = 1.5;
const double mu_b2 = 0.0, s_b2 = 3.2;
const double mu_0 = 0.0, s_0 = 3.2;
const double a_phi = 20.0, b_phi = 1.5;
const double mu_delta = 0.0, s_delta = 0.1;
const double a_tau = 2.5, b_tau = 0.02;      
const double mu_nu = -1.0, s_nu = 1.0;

// Variâncias simples (usadas no gradiente)
const double var_b0    = s_b0 * s_b0;
const double var_b2    = s_b2 * s_b2;
const double var_mu    = s_0 * s_0;
const double var_delta = s_delta * s_delta;
const double var_nu    = s_nu * s_nu;

// Variâncias duplas (2 * sigma^2, usadas na log-posterior)
const double var2_b0    = 2.0 * var_b0;
const double var2_b2    = 2.0 * var_b2;
const double var2_mu    = 2.0 * var_mu;
const double var2_delta = 2.0 * var_delta;
const double var2_nu    = 2.0 * var_nu;


SvmStats AUHMC::make_stats(const vec& y, const vec& h, const vec& l, const int T)
{
  SvmStats st;
  st.T  = T;
  st.h1 = h[0];
  
  const vec ylag = y.subvec(0, T - 1);       // y_0, ..., y_{T-1}
  const vec hh   = h.subvec(0, T - 1);
  const vec ll   = l.subvec(0, T - 1);
  const vec eh   = exp(hh);
  const vec s    = sqrt(ll) % exp(-0.5 * hh);
  const vec s2   = square(s);
  
  // termos t = T
  st.s2_T     = s2[T - 1];
  st.ylag_T   = ylag[T - 1];
  st.lam_T    = ll[T - 1];
  st.lamexp_T = ll[T - 1] * eh[T - 1];
  
  // somas em t = 1..T-1 (indices 0..T-2)
  const vec sa  = s.subvec(0, T - 2);
  const vec s2a = s2.subvec(0, T - 2);
  const vec la  = ll.subvec(0, T - 2);
  const vec ya  = ylag.subvec(0, T - 2);
  const vec eha = eh.subvec(0, T - 2);
  const vec ha  = hh.subvec(0, T - 2);
  
  st.S_s2     = accu(s2a);
  st.S_s2y    = dot(s2a, ya);
  st.S_s2y2   = dot(s2a, square(ya));
  st.S_lam    = accu(la);
  st.S_lamy   = dot(la, ya);
  st.S_lamexp = dot(la, eha);
  st.S_s      = accu(sa);
  st.S_sy     = dot(sa, ya);
  st.S_sexp   = dot(sa, eha);
  st.S_sh     = dot(sa, ha);
  st.S_syh    = accu(sa % ya % ha);
  st.S_sexph  = accu(sa % eha % ha);
  st.S_h      = accu(ha);
  st.S_h2     = dot(ha, ha);
  
  st.sum_logl_minus_l = accu(log(ll) - ll);
  return st;
}

void AUHMC::chol_tridiagonal(const vec& b, const vec& cc, vec& l_diag, 
                             vec& m_subdiag, const int T) 
{
  
  l_diag.zeros(T);
  m_subdiag.zeros(T - 1);
  l_diag[0]=sqrt(b[0]);
  
  for(int j=0; j<T-1 ;j++)
  {
    m_subdiag[j]=cc[j]/l_diag[j];
    l_diag[j+1]=sqrt(b[j+1]-m_subdiag[j]*m_subdiag[j]);
  }
  
}

void AUHMC::apply_chol_tri(const vec& ldiag, const vec& msubdiag, const vec& x,
                           vec& vP, const int T)
{
  
  vP.set_size(T);
  vP[0] = ldiag[0] * x[0];
  for(int i = 1; i < T; ++i) 
  {
    vP[i] = msubdiag[i-1] * x[i-1] + ldiag[i] * x[i];
  }
  
}

// Resolve L * x = b em O(T), onde L é bidiagonal inferior
inline vec AUHMC::solve_L(const vec& ldiag, const vec& msubdiag, const vec& b, 
                          const int T) 
{
  vec x(T);
  x[0] = b[0] / ldiag[0];
  for (int t = 1; t < T; ++t) {
    x[t] = (b[t] - msubdiag[t - 1] * x[t - 1]) / ldiag[t];
  }
  return x;
}

// Resolve L^T * x = b em O(T), onde L^T é bidiagonal superior
inline vec AUHMC::solve_Lt(const vec& ldiag, const vec& msubdiag, const vec& b, 
                           const int T) 
{
  vec x(T);
  x[T - 1] = b[T - 1] / ldiag[T - 1];
  for (int t = T - 2; t >= 0; --t) {
    x[t] = (b[t] - msubdiag[t] * x[t + 1]) / ldiag[t];
  }
  return x;
}

// Resolve M * x = b => (L * L^T) * x = b em O(T)
inline vec AUHMC::solve_M(const vec& ldiag, const vec& msubdiag, const vec& b, 
                          const int T) 
{
  vec y = solve_L(ldiag, msubdiag, b, T);
  return solve_Lt(ldiag, msubdiag, y, T);
}

double AUHMC::jitter_epsilon(const double eps_base, const double Umin, 
                      const double Umax)
{
  return eps_base * (Umin + (Umax-Umin) * R::unif_rand());
}

// #############################################################################
// ############################ theta ##########################################
// #############################################################################
double AUHMC::lpost_theta(const vec& theta, const vec& y, const vec& h, 
                          const vec& l, const int T)
{
  const double b0    = theta[0];
  const double b1    = std::tanh(theta[1]); 
  const double b2    = theta[2];
  const double mu    = theta[3];
  const double phi   = std::tanh(theta[4]);
  const double delta = theta[5];
  const double tau   = std::exp(theta[6]);
  const double logistic_nu = 1.0/(1+std::exp(-theta[7]));
  const double nu          = 40.0*logistic_nu + 2.0*(1.0-logistic_nu);
  
  // Variáveis temporais (e_t, s_t, z_t, A_t)
  const vec y_obs  = y.subvec(1, T);       // y_1, ..., y_T
  const vec y_lag1 = y.subvec(0, T - 1);   // y_0, ..., y_{T-1}
  vec e_t = y_obs - b0 - b1 * y_lag1 - b2 * exp(h.subvec(0, T - 1));
  
  vec s_t = sqrt(l.subvec(0, T - 1)) % exp(-0.5 * h.subvec(0, T - 1));
  vec zt  = e_t % s_t;
  
  // A_t para t = 1, ..., T-1 (índices 0 a T-2)
  vec At = h.subvec(1, T - 1) - mu - phi * (h.subvec(0, T - 2) - mu) - delta * zt.subvec(0, T - 2);
  
  // Componentes da verossimilhança L(theta)
  const double phi2    = phi * phi;
  const double tau2    = tau * tau;
  const double diff_h1 = h[0] - mu;
  
  double term_z  = -0.5 * sum(square(zt));
  double term_h1 = 0.5 * std::log(1.0 - phi2) - ((1.0 - phi2) / (2.0 * tau2)) * (diff_h1 * diff_h1);
  double term_A  = -T * std::log(tau) - (1.0 / (2.0 * tau2)) * sum(square(At));
  
  // Componente t-Student (nu e lambda_t)
  const double half_nu = nu / 2.0;
  double term_nu_l = T * half_nu * std::log(half_nu) - T * std::lgamma(half_nu) +
    (half_nu - 1.0) * sum(log(l.subvec(0, T - 1))) -
    half_nu * sum(l.subvec(0, T - 1));
  
  // Priori dos parâmetros
  const double diff_b0    = b0 - mu_b0;
  const double diff_b2    = b2 - mu_b2;
  const double diff_mu    = mu - mu_0;
  const double diff_delta = delta - mu_delta;
  const double diff_nu    = theta[7] - mu_nu;
  
  double term_priors = - (diff_b0 * diff_b0) / var2_b0
  - (diff_b2 * diff_b2) / var2_b2
  - (diff_mu * diff_mu) / var2_mu
  - (diff_delta * diff_delta) / var2_delta
  - (diff_nu * diff_nu) / var2_nu
  + (a_b1 - 1.0) * std::log(1.0 + b1) + (b_b1 - 1.0) * std::log(1.0 - b1)
  + (a_phi - 1.0) * std::log(1.0 + phi) + (b_phi - 1.0) * std::log(1.0 - phi)
  - (2.0 * a_tau + 1.0) * std::log(tau) - b_tau / tau2;
  
  double L_theta = term_z + term_h1 + term_A + term_nu_l + term_priors;
  
  // Log-Jacobiano (theta[6] = tau* = log(tau))
  double log_J = std::log(1.0 - b1 * b1) + std::log(1.0 - phi2) + theta[6];
  
  return L_theta + log_J;  
}

vec AUHMC::grad_theta(const vec& theta, const vec& y, const vec& h, const vec& l,
                      const int T)
{
  const double b0          = theta[0];
  const double theta1      = theta[1];
  const double b1          = std::tanh(theta1); 
  const double b2          = theta[2];
  const double mu          = theta[3];
  const double theta4      = theta[4];
  const double phi         = std::tanh(theta4);
  const double delta       = theta[5];
  const double tau         = std::exp(theta[6]);
  const double theta7      = theta[7];
  const double logistic_nu = 1.0 / (1 + std::exp(-theta7));
  const double nu          = 40.0 * logistic_nu + 2.0 * (1.0 - logistic_nu);
  
  const vec y_obs  = y.subvec(1, T);       // y_1, ..., y_T
  const vec y_lag1 = y.subvec(0, T - 1);   // y_0, ..., y_{T-1}
  vec exp_h  = exp(h.subvec(0, T - 1));
  vec e_t    = y_obs - b0 - b1 * y_lag1 - b2 * exp_h;
  
  vec s_t = sqrt(l.subvec(0, T - 1)) % exp(-0.5 * h.subvec(0, T - 1));
  vec zt  = e_t % s_t;
  
  vec At = h.subvec(1, T - 1) - mu - phi * (h.subvec(0, T - 2) - mu) - delta * zt.subvec(0, T - 2);
  
  const double phi2       = phi * phi;
  const double tau2       = tau * tau;
  const double diff_h1    = h[0] - mu;
  const double diff_h1_sq = diff_h1 * diff_h1;
  const double d_b1       = 1.0 - phi2;
  
  vec zs_t = zt % s_t;                         
  vec s_A  = s_t.subvec(0, T - 2) % At;         
  
  vec grad(8);
  
  // dL / d b0
  double sum_zs = sum(zs_t);
  double sum_sA = sum(s_A);
  grad[0] = sum_zs - (delta / tau2) * sum_sA - (b0 - mu_b0) / var_b0;
  
  // dL / d b1*
  double sum_sA_ylag = sum(s_A % y_lag1.subvec(0, T - 2));
  double sum_zs_ylag = sum(zs_t % y_lag1);
  double term_chain1 = 1.0 - b1 * b1;
  grad[1] = term_chain1 * (sum_zs_ylag - (delta / tau2) * sum_sA_ylag)
    + (a_b1 - 1.0) * (1.0 - b1) - (b_b1 - 1.0) * (1.0 + b1) - 2.0 * b1;
  
  // dL / d b2
  double sum_zs_exph = sum(zs_t % exp_h);
  double sum_sA_exph = sum(s_A % exp_h.subvec(0, T - 2));
  grad[2] = sum_zs_exph - (delta / tau2) * sum_sA_exph - (b2 - mu_b2) / var_b2;
  
  // dL / d mu
  double sum_At = sum(At);
  grad[3] = (d_b1 / tau2) * diff_h1 + ((1.0 - phi) / tau2) * sum_At - (mu - mu_0) / var_mu;
  
  // dL / d phi*
  double sum_At_hdiff = sum(At % (h.subvec(0, T - 2) - mu));
  grad[4] = - phi + (phi * d_b1 * diff_h1_sq) / tau2 + (d_b1 / tau2) * sum_At_hdiff
  + (a_phi - 1.0) * (1.0 - phi) - (b_phi - 1.0) * (1.0 + phi) - 2.0 * phi;
  
  // dL / d delta
  double sum_At_zt = sum(At % zt.subvec(0, T - 2));
  grad[5] = (1.0 / tau2) * sum_At_zt - (delta - mu_delta) / var_delta;
  
  // dL / d tau*
  double sum_At2 = sum(square(At));
  grad[6] = (d_b1 * diff_h1_sq) / tau2 - T + (1.0 / tau2) * sum_At2
  - (2.0 * a_tau + 1.0) + (2.0 * b_tau) / tau2 + 1.0;
  
  // dL / d nu*
  const double half_nu = nu / 2.0;
  double sum_log_l_minus_l = sum(log(l.subvec(0, T - 1)) - l.subvec(0, T - 1));
  double digamma_half_nu   = R::digamma(half_nu);
  
  double bracket_nu = (T / 2.0) * (1.0 + std::log(half_nu)) - (T / 2.0) * digamma_half_nu
  + 0.5 * sum_log_l_minus_l;
  
  double jacobian_nu = ((40.0 - nu) * (nu - 2.0)) / 38.0;
  grad[7] = jacobian_nu * bracket_nu - (theta7 - mu_nu) / var_nu;
  
  return grad;
}

void AUHMC::I_theta(const vec& theta, const SvmStats& st, mat& I)
{
  I.zeros(8, 8);
  const double T = static_cast<double>(st.T);
  
  // 1. Unpack dos parametros
  const double b1          = std::tanh(theta[1]);
  const double mu          = theta[3];
  const double phi         = std::tanh(theta[4]);
  const double delta       = theta[5];
  const double tau         = std::exp(theta[6]);
  const double logistic_nu = 1.0 / (1.0 + std::exp(-theta[7]));
  const double nu          = 40.0 * logistic_nu + 2.0 * (1.0 - logistic_nu);
  
  // 2. Termos auxiliares
  const double phi2        = phi * phi;
  const double tau2        = tau * tau;
  const double delta2      = delta * delta;
  const double b1_2        = b1 * b1;
  const double d_b1        = 1.0 - phi2;
  const double term_b1_sq  = (1.0 - b1_2) * (1.0 - b1_2);
  const double ratio_tau_d = (tau2 + delta2) / tau2;
  
  const double diff_h1     = st.h1 - mu;
  const double diff_h1_sq  = diff_h1 * diff_h1;
  
  // Termos t = T e somas que nao dependem de theta
  const double s2T       = st.s2_T;
  const double yT        = st.ylag_T;
  const double lamT      = st.lam_T;
  const double lamexpT   = st.lamexp_T;
  
  const double sum_s2           = st.S_s2;
  const double sum_s2_ylag      = st.S_s2y;
  const double sum_s2_ylag2     = st.S_s2y2;
  const double sum_lambda       = st.S_lam;
  const double sum_lambda_ylag  = st.S_lamy;
  const double sum_lambda_exph  = st.S_lamexp;
  const double sum_s            = st.S_s;
  const double sum_s_ylag       = st.S_sy;
  const double sum_s_exph       = st.S_sexp;
  
  // Somas que dependem de mu (expansao polinomial), t = 1..T-1
  const double sum_s_hdiff       = st.S_sh    - mu * st.S_s;
  const double sum_s_ylag_hdiff  = st.S_syh   - mu * st.S_sy;
  const double sum_s_exph_hdiff  = st.S_sexph - mu * st.S_sexp;
  const double sum_hdiff         = st.S_h     - mu * (T - 1.0);
  const double sum_hdiff2        = st.S_h2 - 2.0 * mu * st.S_h + mu * mu * (T - 1.0);
  
  // =============================================================
  // 1. Matriz de informacao I(theta*)
  // =============================================================
  
  // Bloco (b0, .)
  I(0, 0) = s2T + ratio_tau_d * sum_s2 + (1.0 / var_b0);
  I(0, 1) = (1.0 - b1_2) * (s2T * yT + ratio_tau_d * sum_s2_ylag);
  I(0, 2) = lamT + ratio_tau_d * sum_lambda;
  I(0, 3) = - (1.0 - phi) * delta / tau2 * sum_s;
  I(0, 4) = - d_b1 * delta / tau2 * sum_s_hdiff;
  
  // Bloco (b1*, .)
  I(1, 0) = I(0, 1);
  I(1, 1) = term_b1_sq * s2T * (yT * yT)
    + term_b1_sq * ratio_tau_d * sum_s2_ylag2
  + (a_b1 + b_b1) * (1.0 - b1_2);
  I(1, 2) = (1.0 - b1_2) * (yT * lamT + ratio_tau_d * sum_lambda_ylag);
  I(1, 3) = - (1.0 - b1_2) * delta * (1.0 - phi) / tau2 * sum_s_ylag;
  I(1, 4) = - (1.0 - b1_2) * delta * d_b1 / tau2 * sum_s_ylag_hdiff;
  
  // Bloco (b2, .)
  I(2, 0) = I(0, 2);
  I(2, 1) = I(1, 2);
  I(2, 2) = lamexpT + ratio_tau_d * sum_lambda_exph + (1.0 / var_b2);
  I(2, 3) = - (1.0 - phi) * delta / tau2 * sum_s_exph;
  I(2, 4) = - d_b1 * delta / tau2 * sum_s_exph_hdiff;
  
  // Bloco (mu, .)
  I(3, 0) = I(0, 3);
  I(3, 1) = I(1, 3);
  I(3, 2) = I(2, 3);
  I(3, 3) = (d_b1 / tau2) + ((1.0 - phi) * (1.0 - phi) * (T - 1.0)) / tau2 + (1.0 / var_mu);
  I(3, 4) = (2.0 * phi * d_b1 * diff_h1) / tau2 + (d_b1 * (1.0 - phi) / tau2) * sum_hdiff;
  I(3, 6) = (2.0 * d_b1 * diff_h1) / tau2;
  
  // Bloco (phi*, .)
  I(4, 0) = I(0, 4);
  I(4, 1) = I(1, 4);
  I(4, 2) = I(2, 4);
  I(4, 3) = I(3, 4);
  I(4, 4) = - (diff_h1_sq * d_b1 * d_b1) / tau2
  + (2.0 * phi2 * d_b1 * diff_h1_sq) / tau2
  + (d_b1 * d_b1 / tau2) * sum_hdiff2
  + (1.0 - phi2) * (a_phi + b_phi + 1.0);
  I(4, 6) = (2.0 * phi * d_b1 * diff_h1_sq) / tau2;
  
  // Bloco (delta, .)
  I(5, 5) = ((T - 1.0) / tau2) + (1.0 / var_delta);
  
  // Bloco (tau*, .)
  I(6, 3) = I(3, 6);
  I(6, 4) = I(4, 6);
  I(6, 6) = (2.0 / tau2) * (d_b1 * diff_h1_sq + (T - 1.0) * tau2 + 2.0 * b_tau);
  
  // Bloco (nu*, .)  (versao -H com lambda fixo)
  const double half_nu      = nu / 2.0;
  const double digamma_half_nu  = R::digamma(half_nu);
  const double trigamma_half_nu = R::trigamma(half_nu);
  
  const double d_nu         = (T / 2.0) * (1.0 + std::log(half_nu)) - (T / 2.0) * digamma_half_nu
  + 0.5 * st.sum_logl_minus_l;
  const double d_prime_nu   = (T / (2.0 * nu)) - (T / 4.0) * trigamma_half_nu;
  const double c_nu         = ((40.0 - nu) * (nu - 2.0)) / 38.0;
  const double c_prime_nu   = (21.0 - nu) / 19.0;
  
  I(7, 7) = - (c_nu * (c_prime_nu * d_nu + c_nu * d_prime_nu) - (1.0 / var_nu));
  
  //I = 0.5 * (I + I.t());
  //I.diag() += 1e-4;
  
  
}

vec AUHMC::updatetheta(const vec& theta_cur, const vec& y, const vec& h_cur, 
                       const vec& l_cur, const int T, const int nlpfrog, 
                       const double epsilon, int& acc, 
                       int& total_j_iters, int& non_conv_count)
{
  
  // 1. Checagens de tamanho de vetor
  if (y.n_elem != static_cast<uword>(T + 1))
    Rcpp::stop("y deve ter T+1 elementos: (y0, y1, ..., yT).");
  if (l_cur.n_elem < static_cast<uword>(T) || h_cur.n_elem < static_cast<uword>(T))
    Rcpp::stop("l_star e h devem ter pelo menos T elementos.");
  
  const int d = theta_cur.n_elem; 
  SvmStats st_cur = make_stats(y, h_cur, l_cur, T);
  
  const int max_fixed_point_iters = 1e2;
  const double tol = 1e-5;
  double norm1 = 0.0, norm2 = 0.0;
  
  mat F_cur(d, d), M_cur(d, d), M_prop(d, d), M_star(d, d);
  I_theta(theta_cur, st_cur, F_cur);
  M_cur = F_cur;
  
  mat L_M(d, d), invM(d, d);
  if(!arma::chol(L_M, M_cur, "lower")) 
  {
    non_conv_count++;  
    return theta_cur;
  }
  const vec z_norm = arma::randn<vec>(d);
  vec p_cur  = L_M * z_norm;
  
  vec theta_prop = arma::zeros<vec>(d);
  vec p_prop = theta_prop;
  
  vec theta_prev = theta_prop;
  vec p_prev = p_prop;
  
  vec g(d);
  
  int j = 0;
  bool converged = false;
  
  // Loop de Ponto Fixo (j-loop) para adaptação da matriz de massa M
  while(!converged && j < max_fixed_point_iters)
  {
    total_j_iters++;
    
    if(!arma::inv_sympd(invM, M_cur)) 
    {
      non_conv_count++;  
      return theta_cur;
    }
    
    M_star = M_cur;
    theta_prop = theta_cur;
    p_prop = p_cur;
    g = grad_theta(theta_prop, y, h_cur, l_cur, T);
    
    for(int k = 0; k < nlpfrog; ++k)
    {
      p_prop += (0.5 * epsilon) * g;
      theta_prop += epsilon * (invM * p_prop);
      g = grad_theta(theta_prop, y, h_cur, l_cur, T);
      p_prop += (0.5 * epsilon) * g;
    }
    
    if(!theta_prop.is_finite() || !p_prop.is_finite()) 
    {
      non_conv_count++;  
      return theta_cur;
    }
    
    I_theta(theta_prop, st_cur, M_prop);
    M_cur = 0.5 * (F_cur + M_prop);
    
    // Critério de convergência do Ponto Fixo
    norm1 = arma::norm(theta_prop-theta_prev);
    norm2 = arma::norm(p_prop-p_prev);
    
    if( j > 0 && norm1 <= tol && norm2 <= tol) 
    {
      converged = true;
    }
    
    theta_prev = theta_prop;
    p_prev = p_prop;
    
    j++;
  }
  
  if(!converged) {
    non_conv_count++;
    return theta_cur;
  }
  
  // Energia com M_star: invM é a inversa da M da última trajetória (não foi alterada após o laço)
  const double lp_cur  = lpost_theta(theta_cur, y, h_cur, l_cur, T);
  const double lp_prop = lpost_theta(theta_prop, y, h_cur, l_cur, T);
  const double H_init  = lp_cur  - 0.5 * arma::dot(p_cur,  invM * p_cur);
  const double H_final = lp_prop - 0.5 * arma::dot(p_prop, invM * p_prop);
  const double log_alpha = H_final - H_init;
  
  if (!std::isfinite(log_alpha)) return theta_cur;
  
  if( std::log(R::unif_rand()) < log_alpha ) 
  {
    acc += 1;
    return theta_prop;
  }else{
    return theta_cur;
  }
  
} 

// #############################################################################
// ############################ h ##############################################
// #############################################################################

double AUHMC::lpost_h(const vec& h, const vec& theta, const vec& y, const vec& l, 
                      const int T)
{
  // 1. Unpack dos parâmetros de theta
  const double b0    = theta[0];
  const double b1    = std::tanh(theta[1]); 
  const double b2    = theta[2];
  const double mu    = theta[3];
  const double phi   = std::tanh(theta[4]);
  const double delta = theta[5];
  const double tau   = std::exp(theta[6]);
  
  const double phi2  = phi * phi;
  const double tau2  = tau * tau;
  
  // 2. Alinhamento de y, h e l
  // y possui T+1 elementos: y[0]=y0, y[1..T]=y1..yT
  // h e l possuem T elementos: h[0..T-1]=h1..hT
  const vec y_obs  = y.subvec(1, T);       // y_1, ..., y_T
  const vec y_lag1 = y.subvec(0, T - 1);   // y_0, ..., y_{T-1}
  const vec l_sub  = l.subvec(0, T - 1);   // \lambda_1, ..., \lambda_T
  
  // 3. Resíduos e erros padronizados z_t (t = 1..T)
  vec exp_h = exp(h);                      // exp(h_t) para t = 1..T
  vec e_t   = y_obs - b0 - b1 * y_lag1 - b2 * exp_h;
  
  vec s_t   = sqrt(l_sub) % exp(-0.5 * h);
  vec zt    = e_t % s_t;                   // z_1, ..., z_T (índices 0..T-1)
  
  // 4. Componentes da log-posterior L(h_T)
  const double term_sum_h = -0.5 * sum(h);
  const double term_z     = -0.5 * sum(square(zt));
  
  // h_1 é o índice 0 de h
  const double diff_h1    = h[0] - mu;
  const double term_h1    = - ((1.0 - phi2) / (2.0 * tau2)) * (diff_h1 * diff_h1);
  
  // A_t = h_{t+1} - mu - phi*(h_t - mu) - delta*z_t (para t = 1..T-1, índices 0..T-2)
  vec At = h.subvec(1, T - 1) - mu - phi * (h.subvec(0, T - 2) - mu) - delta * zt.subvec(0, T - 2);
  const double term_A     = - (1.0 / (2.0 * tau2)) * sum(square(At));
  
  return term_sum_h + term_z + term_h1 + term_A;
}

vec AUHMC::grad_h(const vec& h, const vec& theta, const vec& y, const vec& l, 
                  const int T)
{
  // 1. Unpack dos parâmetros de theta
  const double b0    = theta[0];
  const double b1    = std::tanh(theta[1]); 
  const double b2    = theta[2];
  const double mu    = theta[3];
  const double phi   = std::tanh(theta[4]);
  const double delta = theta[5];
  const double tau   = std::exp(theta[6]);
  
  const double phi2  = phi * phi;
  const double tau2  = tau * tau;
  
  // 2. Ajuste dos vetores observados (y_0, y_1, ..., y_T)
  const vec y_obs  = y.subvec(1, T);       // y_1, ..., y_T
  const vec y_lag1 = y.subvec(0, T - 1);   // y_0, ..., y_{T-1}
  const vec l_sub  = l.subvec(0, T - 1);   // \lambda_1, ..., \lambda_T
  
  // 3. Média condicional \mu_t = \beta_0 + \beta_1 y_{t-1} + \beta_2 e^{h_t}
  vec exp_h = exp(h);
  vec mu_t  = b0 + b1 * y_lag1 + b2 * exp_h;
  
  // 4. Termos z_t (t = 1..T, índices 0..T-1)
  vec s_t   = sqrt(l_sub) % exp(-0.5 * h);
  vec e_t   = y_obs - mu_t;
  vec zt    = e_t % s_t;
  vec z2_t  = square(zt);
  
  // 5. Termos A_t e B_t para t = 1..T-1 (índices 0..T-2)
  vec At = h.subvec(1, T - 1) - mu - phi * (h.subvec(0, T - 2) - mu) - delta * zt.subvec(0, T - 2);
  vec Bt = delta * b2 * exp(0.5 * h.subvec(0, T - 2)) % sqrt(l_sub.subvec(0, T - 2)) 
    + 0.5 * delta * zt.subvec(0, T - 2) - phi;
  
  vec grad(T);
  
  // Base comum para todo t: -0.5 + b2 * lambda_t * (y_t - mu_t) + 0.5 * z_t^2
  vec base_term = -0.5 + b2 * l_sub % (y_obs - mu_t) + 0.5 * z2_t;
  
  // Para t = 1 (índice 0)
  const double diff_h1 = h[0] - mu;
  grad[0] = base_term[0] - ((1.0 - phi2) / tau2) * diff_h1 - (1.0 / tau2) * At[0] * Bt[0];
  
  // Para t = 2..T-1 (índices 1..T-2)
  for(int t = 1; t < T - 1; ++t) {
    grad[t] = base_term[t] - (1.0 / tau2) * (At[t] * Bt[t] + At[t - 1]);
  }
  
  // Para t = T (índice T-1)
  grad[T - 1] = base_term[T - 1] - (1.0 / tau2) * At[T - 2];
  
  
  return grad;
}

bool AUHMC::I_h(const vec& h, const vec& theta, const vec& l, const int T,
                vec& ldiag_I, vec& msubdiag_I, vec& ldiag, vec& msubdiag)
{
  // 1. Unpack dos parâmetros de theta
  const double b2    = theta[2];
  const double phi   = std::tanh(theta[4]);
  const double delta = theta[5];
  const double tau   = std::exp(theta[6]);
  
  const double phi2   = phi * phi;
  const double tau2   = tau * tau;
  const double delta2 = delta * delta;
  
  // 2. Pré-alocação dos vetores da matriz I(h_T)
  ldiag_I.zeros(T);
  msubdiag_I.zeros(T - 1);
  
  const vec l_sub = l.subvec(0, T - 1);
  vec exp_h = exp(h);
  vec ct    = b2 * sqrt(l_sub) % exp(0.5 * h); // c_t = b2 * lambda_t^{1/2} * e^{h_t/2}
  
  // 3. Montagem do vetor w
  vec w(T);
  
  // Para t = 1 (índice 0)
  const double term_c1 = delta * ct[0] - phi;
  w[0] = (1.0 / tau2) * (0.25 * delta2 + term_c1 * term_c1 + (1.0 - phi2));
  
  // Para t = 2..T-1 (índices 1..T-2)
  for(int t = 1; t < T - 1; ++t) 
  {
    const double term_ct = delta * ct[t] - phi;
    w[t] = (1.0 / tau2) * (0.25 * delta2 + term_ct * term_ct + 1.0);
  }
  
  // Para t = T (índice T-1)
  w[T - 1] = 1.0 / tau2;
  
  // 4. Diagonal principal de I(h_T): beta2^2 * lambda_t * e^{h_t} + 0.5 + w_t
  ldiag_I = b2*b2 * l_sub % exp_h + 0.5 + w;
  
  // 5. Subdiagonal de I(h_T) para t = 1..T-1 (índices 0..T-2)
  for(int t = 0; t < T - 1; ++t) 
  {
    msubdiag_I[t] = (1.0 / tau2) * (delta * ct[t] - phi);
  }
  
  // 6. Fatoração de Cholesky bidiagonal de I(h_T) = L L^T
  ldiag.zeros(T);
  msubdiag.zeros(T - 1);
  
  //if (ldiag_I[0] <= 0.0) return false;
  if (!(ldiag_I[0] > 0.0) || !std::isfinite(ldiag_I[0])) return false;
  ldiag[0] = std::sqrt(ldiag_I[0]);
  
  for (int j = 0; j < T - 1; ++j) {
    msubdiag[j] = msubdiag_I[j] / ldiag[j];
    const double val = ldiag_I[j + 1] - msubdiag[j] * msubdiag[j];
    if (val <= 0.0 || !std::isfinite(val)) return false; // Falha na Cholesky (não def. positiva)
    ldiag[j + 1] = std::sqrt(val);
  }
  
  return true;
}


vec AUHMC::updateh(const vec& theta_cur, const vec& y, const vec& h_cur, 
                   const vec& l_cur, const int T, const int nlpfrog, 
                   const double epsilon, int& acc, 
                   int& total_j_iters, int& non_conv_count)
{
  // 1. Checagens de tamanho de vetor
  if (y.n_elem != static_cast<uword>(T + 1))
    Rcpp::stop("y deve ter T+1 elementos: (y0, y1, ..., yT).");
  if (h_cur.n_elem != static_cast<uword>(T) || l_cur.n_elem < static_cast<uword>(T))
    Rcpp::stop("h deve ter exatamente T elementos e lambda pelo menos T.");
  
  const int    max_fixed_point_iters = 100;
  const double tol = 1e-5;
  
  // 2. F(h_r) tridiagonal (e sua Cholesky, usada para sorteiar p_cur)
  vec ldiag_F(T), msub_F(T - 1), ldiag_FL, msub_FL;
  if (!I_h(h_cur, theta_cur, l_cur, T, ldiag_F, msub_F, ldiag_FL, msub_FL)) {
    non_conv_count++;
    return h_cur;
  }
  
  // M_0 = F(h_r)
  vec ldiag_Mi = ldiag_F;                 // diagonal de M_cur
  vec msub_Mi  = msub_F;                  // subdiagonal de M_cur
  vec ldiag_M(T), msub_M(T - 1);          // Cholesky da M usada na trajetória (M_star)
  
  // Momentum inicial: gamma_0 ~ N(0, F(h_r)), sorteado uma única vez
  const vec z_norm = arma::randn<vec>(T);
  vec p_cur(T);
  apply_chol_tri(ldiag_FL, msub_FL, z_norm, p_cur, T);
  
  vec h_prop = arma::zeros<vec>(T);
  vec p_prop = h_prop;
  vec h_prev = h_prop;
  vec p_prev = h_prop;
  vec g_k(T);
  
  vec ldiag_P(T), msub_P(T - 1), ldiag_dummy, msub_dummy;   // F(h_prop)
  
  int  j = 0;
  bool converged = false;
  
  // 3. Loop do ponto fixo (j-loop) para adaptação da matriz de massa M
  while (!converged && j < max_fixed_point_iters)
  {
    total_j_iters++;
    
    // Cholesky tridiagonal da M desta trajetória, em O(T)
    chol_tridiagonal(ldiag_Mi, msub_Mi, ldiag_M, msub_M, T);
    if (!ldiag_M.is_finite() || !msub_M.is_finite()) {
      non_conv_count++;
      return h_cur;
    }
    
    // Cada trajetória recomeça em h_r, com o mesmo momentum inicial
    h_prop = h_cur;
    p_prop = p_cur;
    
    // 4. Leapfrog com gradiente reaproveitado (nlpfrog + 1 chamadas)
    g_k = grad_h(h_prop, theta_cur, y, l_cur, T);
    for (int k = 0; k < nlpfrog; ++k)
    {
      p_prop += (0.5 * epsilon) * g_k;
      h_prop += epsilon * solve_M(ldiag_M, msub_M, p_prop, T);
      g_k     = grad_h(h_prop, theta_cur, y, l_cur, T);
      p_prop += (0.5 * epsilon) * g_k;
    }
    
    if (!h_prop.is_finite() || !p_prop.is_finite()) {
      non_conv_count++;
      return h_cur;
    }
    
    // M_{j+1} = 0.5 * [F(h_r) + F(h^{L,j})]
    if (!I_h(h_prop, theta_cur, l_cur, T, ldiag_P, msub_P, ldiag_dummy, msub_dummy)) {
      non_conv_count++;
      return h_cur;
    }
    ldiag_Mi = 0.5 * (ldiag_F + ldiag_P);
    msub_Mi  = 0.5 * (msub_F  + msub_P);
    
    // Convergência do ponto fixo: compara com a iteração anterior (tolerância relativa)
    if (j > 0 && arma::norm(h_prop - h_prev) <= tol && arma::norm(p_prop - p_prev) <= tol);
    {
      converged = true;
    }
    
    h_prev = h_prop;
    p_prev = p_prop;
    
    j++;
  }
  
  if (!converged) {
    non_conv_count++;
    return h_cur;
  }
  
  // 5. Aceitação Metropolis-Hastings: energia com M_star (a M da última trajetória)
  // ldiag_M / msub_M não foram recalculados depois do laço, logo são a Cholesky de M_star
  const double lp_cur  = lpost_h(h_cur,  theta_cur, y, l_cur, T);
  const double lp_prop = lpost_h(h_prop, theta_cur, y, l_cur, T);
  
  const vec u_c = solve_L(ldiag_M, msub_M, p_cur,  T);
  const vec u_p = solve_L(ldiag_M, msub_M, p_prop, T);
  
  const double H_init  = lp_cur  - 0.5 * arma::dot(u_c, u_c);
  const double H_final = lp_prop - 0.5 * arma::dot(u_p, u_p);
  const double log_alpha = H_final - H_init;
  
  if (!std::isfinite(log_alpha)) return h_cur;
  
  if (std::log(R::unif_rand()) < log_alpha) {
    acc += 1;
    return h_prop;
  }
  return h_cur;
}
  
// #############################################################################
// ############################ l ##############################################
// #############################################################################

double AUHMC::lpost_l(const vec& l_star, const vec& h, const vec& theta, 
                      const vec& y, const int T)
{
  // 1. Unpack dos parâmetros de theta (com reparametrização logística para nu em (2, 40))
  const double b0          = theta[0];
  const double b1          = std::tanh(theta[1]);
  const double b2          = theta[2];
  const double mu          = theta[3];
  const double phi         = std::tanh(theta[4]);
  const double delta       = theta[5];
  const double tau         = std::exp(theta[6]);
  const double logistic_nu = 1.0 / (1.0 + std::exp(-theta[7]));
  const double nu          = 40.0 * logistic_nu + 2.0 * (1.0 - logistic_nu);
  
  const double tau2 = tau * tau;
  
  // 2. Sub-vetores dos dados observados
  const vec l = exp(l_star);
  const vec y_obs  = y.subvec(1, T);       // y_1, ..., y_T
  const vec y_lag1 = y.subvec(0, T - 1);   // y_0, ..., y_{T-1}
  const vec l_sub  = l.subvec(0, T - 1);   // \lambda_1, ..., \lambda_T
  
  // 3. Média condicional \mu_t = \beta_0 + \beta_1 y_{t-1} + \beta_2 e^{h_t}
  vec exp_h = exp(h);
  vec mu_t  = b0 + b1 * y_lag1 + b2 * exp_h;
  
  // 4. Resíduos padronizados z_t = (y_t - \mu_t) * \sqrt{\lambda_t} * exp(-0.5 * h_t)
  vec e_t = y_obs - mu_t;
  vec s_t = sqrt(l_sub) % exp(-0.5 * h);
  vec zt  = e_t % s_t;
  
  // 5. Termo A_t para t = 1..T-1 (índices 0..T-2)
  vec At = h.subvec(1, T - 1) - mu - phi * (h.subvec(0, T - 2) - mu) - delta * zt.subvec(0, T - 2);
  
  // 6. Acumulação da Log-Posterior \mathcal{L}(\lambda)
  const double term_logl    = 0.5 * (nu + 1.0) * arma::sum(arma::log(l_sub));
  const double term_z2      = -0.5 * arma::dot(zt, zt);
  const double term_A2      = - (0.5 / tau2) * arma::dot(At, At);
  const double term_prior_l = -0.5 * nu * arma::sum(l_sub);
  
  return term_logl + term_z2 + term_A2 + term_prior_l;
}

vec AUHMC::grad_l(const vec& l_star, const vec& h, const vec& theta, const vec& y, 
                  const int T)
{
  const double b0          = theta[0];
  const double b1          = std::tanh(theta[1]);
  const double b2          = theta[2];
  const double mu          = theta[3];
  const double phi         = std::tanh(theta[4]);
  const double delta       = theta[5];
  const double tau         = std::exp(theta[6]);
  const double logistic_nu = 1.0 / (1.0 + std::exp(-theta[7]));
  const double nu          = 40.0 * logistic_nu + 2.0 * (1.0 - logistic_nu);
  const double tau2        = tau * tau;
  
  
  const vec l = exp(l_star);
  const vec y_obs  = y.subvec(1, T);
  const vec y_lag1 = y.subvec(0, T - 1);
  const vec l_sub  = l.subvec(0, T - 1);
  
  const vec mu_t = b0 + b1 * y_lag1 + b2 * exp(h);
  const vec s_t  = sqrt(l_sub) % exp(-0.5 * h);
  const vec zt   = (y_obs - mu_t) % s_t;
  const vec At   = h.subvec(1, T - 1) - mu - phi * (h.subvec(0, T - 2) - mu)
    - delta * zt.subvec(0, T - 2);
  
  // dL/dl*_t = (nu+1)/2 - z_t^2/2 - nu*lambda_t/2 + delta z_t A_t/(2 tau^2)  (t <= T-1)
  vec grad = 0.5 * (nu + 1.0) - 0.5 * square(zt) - 0.5 * nu * l_sub;
  grad.subvec(0, T-2) += (0.5 * delta / tau2) * zt.subvec(0, T - 2) % At;
  
  return grad;
}

vec AUHMC::I_l(const vec& l_star, const vec& theta, const int T)
{
  const double delta       = theta[5];
  const double tau         = std::exp(theta[6]);
  const double logistic_nu = 1.0 / (1.0 + std::exp(-theta[7]));
  const double nu          = 40.0 * logistic_nu + 2.0 * (1.0 - logistic_nu);
  const double tau2        = tau * tau;
  
  const vec l_sub = exp(l_star.subvec(0, T - 1));
  
  // Diagonal principal de I(l*): 0.5 + 0.5 * nu * lambda_t + k_t
  vec diag_I = 0.5 + 0.5 * nu * l_sub;
  
  // k_t = delta^2 / (4 * tau^2) para t = 1..T-1 (índices 0..T-2)
  diag_I.subvec(0, T - 2) += (0.25 * delta * delta) / tau2;
  
  return diag_I;
}

vec AUHMC::updatel(const vec& l_star_cur, const vec& h_cur, const vec& theta_cur,
                   const vec& y, const int T, const int nlpfrog,
                   const double epsilon, int& acc,
                   int& total_j_iters, int& non_conv_count)
{
  // 1. Checagens de tamanho de vetor
  if (y.n_elem != static_cast<uword>(T + 1))
    Rcpp::stop("y deve ter T+1 elementos: (y0, y1, ..., yT).");
  if (l_star_cur.n_elem < static_cast<uword>(T) || h_cur.n_elem < static_cast<uword>(T))
    Rcpp::stop("l_star e h devem ter pelo menos T elementos.");
  
  // 2. Matriz de Informação de Fisher diagonal inicial I(l_star_cur)
  vec diag_I_curr = I_l(l_star_cur, theta_cur, T);
  if (!diag_I_curr.is_finite() || arma::any(diag_I_curr <= 0.0)) {
    non_conv_count++;
    return l_star_cur;
  }
  
  const int max_fixed_point_iters = 1e2;
  const double tol = 1e-5;
  double norm1 = 0.0, norm2 = 0.0;
  
  const vec F_cur = I_l(l_star_cur, theta_cur, T);
  if (!F_cur.is_finite() || arma::any(F_cur <= 0.0)) {
    non_conv_count++;
    return l_star_cur;
  }
  
  vec M_cur = F_cur;
  vec M_prop(T), M_star(T), inv_M(T);
  
  const vec z_norm = arma::randn<vec>(T);
  const vec p_cur  = sqrt(M_cur) % z_norm;
  
  vec l_star_prop = arma::zeros<vec>(T);
  vec p_prop      = l_star_prop;
  
  vec l_star_prev = l_star_prop;
  vec p_prev      = p_prop;
  
  int  j = 0;
  bool converged = false;
  
  while (!converged && j < max_fixed_point_iters)
  {
    total_j_iters++;
    
    if (!M_cur.is_finite() || arma::any(M_cur <= 0.0)) {
      non_conv_count++;
      return l_star_cur;
    }
    
    inv_M  = 1.0 / M_cur;
    M_star = M_cur;                        // M usada nesta trajetória
    
    // Cada trajetória recomeça em l*_r, com o mesmo momentum inicial
    l_star_prop = l_star_cur;
    p_prop      = p_cur;
    
    for (int k = 0; k < nlpfrog; ++k)
    {
      p_prop      += (0.5 * epsilon) * grad_l(l_star_prop, h_cur, theta_cur, y, T);
      l_star_prop += epsilon * (inv_M % p_prop);
      p_prop      += (0.5 * epsilon) * grad_l(l_star_prop, h_cur, theta_cur, y, T);
    }
    
    if (!l_star_prop.is_finite() || !p_prop.is_finite()) {
      non_conv_count++;
      return l_star_cur;
    }
    
    // M_{j+1} = 0.5 * [F(l*_r) + F(l*^{L,j})]
    M_prop = I_l(l_star_prop, theta_cur, T);
    M_cur  = 0.5 * (F_cur + M_prop);
    
    // Convergência do ponto fixo: compara com a iteração anterior
    norm1 = arma::norm(l_star_prop - l_star_prev);
    norm2 = arma::norm(p_prop - p_prev);
    
    if (j > 0 && norm1 <= tol && norm2 <= tol) {
      converged = true;
    }
    
    l_star_prev = l_star_prop;
    p_prev      = p_prop;
    
    j++;
  }
  
  if (!converged) {
    non_conv_count++;
    return l_star_cur;
  }
  
  // 4. Aceitação Metropolis-Hastings: energia com M_star nos dois lados
  const double lp_cur  = lpost_l(l_star_cur,  h_cur, theta_cur, y, T);
  const double lp_prop = lpost_l(l_star_prop, h_cur, theta_cur, y, T);
  const double H_init  = lp_cur  - 0.5 * arma::dot(p_cur,  inv_M % p_cur);
  const double H_final = lp_prop - 0.5 * arma::dot(p_prop, inv_M % p_prop);
  const double log_alpha = H_final - H_init;
  
  if (!std::isfinite(log_alpha)) return l_star_cur;
  
  if (std::log(R::unif_rand()) < log_alpha) {
    acc += 1;
    return l_star_prop;
  }else{
    return l_star_cur;
  }
  
}


