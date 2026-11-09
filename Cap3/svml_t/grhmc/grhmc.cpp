#include "grhmc.h"
#include <cmath>

// Parameters priori
double mu_b0 = 0.0, s_b0 = 3.2;
double a_b1 = 5.0, b_b1 = 1.5;
double mu_b2 = 0.0, s_b2 = 3.2;
double mu_0 = 0.0, s_0 = 3.2;
double a_phi = 20.0, b_phi = 1.5;
double mu_delta = 0.0, s_delta = 0.1;
double a_tau = 2.5, b_tau = 0.02;      
double mu_nu = -1.0, s_nu = 1.0;

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

SvmStats GRHMC::make_stats(const vec& y, const vec& h, const vec& l, const int T)
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

void GRHMC::chol_tridiagonal(const vec& b, const vec& cc, vec& l_diag, 
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

void GRHMC::apply_chol_tri(const vec& ldiag, const vec& msubdiag, const vec& x,
                             vec& vP, const int T)
{
  
  vP.set_size(T);
  vP[0] = ldiag[0] * x[0];
  for(int i = 1; i < T; ++i) 
  {
    vP[i] = msubdiag[i-1] * x[i-1] + ldiag[i] * x[i];
  }
  
}


double GRHMC::sample_r(const double m, const double c, const int d)
{
  const double mc    = m * c;
  const double mc2   = m * c * c;
  const double shape = static_cast<double>(d);
  const double scale = 1.0 / c;
  
  // Gamma: shape = d, scale = 1/c (rate = c)
  //const distr_param param(static_cast<double>(d), scale);
  
  while(true) 
  {
    double r = randg( distr_param(shape, scale) );
    const double x = r / mc;
    const double log_ratio = -mc2 / (std::sqrt(x * x + 1.0) + x);
    if (std::log(randu<double>()) <= log_ratio) return r;
  }
}

vec GRHMC::sample_pI(const double m, const double c, const int d)
{
  const double r = sample_r(m, c, d);
  const vec z = randn<vec>(d);
  return (r / norm(z, 2)) * z;
}


// #############################################################################
// ############################ theta ##########################################
// #############################################################################
double GRHMC::lpost_theta(const vec& theta, const vec& y, const vec& h, 
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

vec GRHMC::grad_theta(const vec& theta, const vec& y, const vec& h, const vec& l,
                     const int T, const mat& Lchol, const mat& Iparb1, 
                     const mat& Iparmu, const mat& Iparphi, const mat& Ipardelta, 
                     const mat& Ipartau, const mat& Iparnu)
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
  
  // Termo -0.5 * tr(G^{-1} dG/dtheta_i)
  const mat Xb1    = solve(trimatu(Lchol.t()), solve(trimatl(Lchol), Iparb1));
  const mat Xmu    = solve(trimatu(Lchol.t()), solve(trimatl(Lchol), Iparmu));
  const mat Xphi   = solve(trimatu(Lchol.t()), solve(trimatl(Lchol), Iparphi));
  const mat Xdelta = solve(trimatu(Lchol.t()), solve(trimatl(Lchol), Ipardelta));
  const mat Xtau   = solve(trimatu(Lchol.t()), solve(trimatl(Lchol), Ipartau));
  const mat Xnu    = solve(trimatu(Lchol.t()), solve(trimatl(Lchol), Iparnu));
  
  //grad[0] -= 0   (dI/d beta_0 = 0)
  grad[1] -= 0.5 * trace(Xb1);
  //grad[2] -= 0   (dI/d beta_2 = 0)
  grad[3] -= 0.5 * trace(Xmu);
  grad[4] -= 0.5 * trace(Xphi);
  grad[5] -= 0.5 * trace(Xdelta);
  grad[6] -= 0.5 * trace(Xtau);
  grad[7] -= 0.5 * trace(Xnu);
  
  return grad;
}

vec GRHMC::nuHtheta(const mat& Lchol, const mat& Iparb1, const mat& Iparmu, 
                   const mat& Iparphi, const mat& Ipardelta, const mat& Ipartau, 
                   const mat& Iparnu, const vec& p, const double m, 
                   const double c) 
{
  vec u   = solve(trimatl(Lchol), p);       // L u = p
  vec iGp = solve(trimatu(Lchol.t()), u);   // L^T x = u
  
  vec nuH = zeros<vec>(8);
  
  //nuH[0] = 0
  nuH[1] = dot(iGp, Iparb1 * iGp);
  //nuH[2] = 0
  nuH[3] = dot(iGp, Iparmu * iGp);
  nuH[4] = dot(iGp, Iparphi * iGp);
  nuH[5] = dot(iGp, Ipardelta * iGp);
  nuH[6] = dot(iGp, Ipartau * iGp);
  nuH[7] = dot(iGp, Iparnu * iGp);
  
  double M = m * std::sqrt(dot(u, u) / (m * m * c * c) + 1.0);
  
  return (1/M) * nuH ;
}

vec GRHMC::gradHmomtheta(const mat& Lchol, const vec& p, const double m, 
                         const double c)
{
  vec u = solve(trimatl(Lchol), p);
  double M = m * std::sqrt(dot(u, u) / (m * m * c * c) + 1.0);
  
  return (1/M) * solve(trimatu(Lchol.t()), u);
}

void GRHMC::I_theta_ref(const vec& theta, const vec& y, const vec& h, const vec& l, 
                       const int T, mat& I, mat& Iparb1, mat& Iparmu, mat& Iparphi, 
                       mat& Ipardelta, mat& Ipartau, mat& Iparnu)
{
  I.zeros(8, 8);
  Iparb1.zeros(8, 8);
  Iparmu.zeros(8, 8);
  Iparphi.zeros(8, 8);
  Ipardelta.zeros(8, 8);
  Ipartau.zeros(8, 8);
  Iparnu.zeros(8, 8);
  
  // 1. Unpack dos parametros
  const double theta1      = theta[1]; // beta_1*
  const double b1          = std::tanh(theta1);
  const double mu          = theta[3];
  const double theta4      = theta[4]; // phi*
  const double phi         = std::tanh(theta4);
  const double delta       = theta[5];
  const double tau         = std::exp(theta[6]);
  const double theta7      = theta[7]; // nu*
  const double logistic_nu = 1.0/(1+std::exp(-theta7));
  const double nu          = 40.0 * logistic_nu + 2.0 * (1.0 - logistic_nu);
  
  // 2. Termos auxiliares
  const double phi2        = phi * phi;
  const double tau2        = tau * tau;
  const double delta2      = delta * delta;
  const double b1_2        = b1 * b1;
  const double d_b1        = 1.0 - phi2;
  const double term_b1_sq  = (1.0 - b1_2) * (1.0 - b1_2);
  const double ratio_tau_d = (tau2 + delta2) / tau2;
  
  const vec y_lag1         = y.subvec(0, T - 1);          // y_0, ..., y_{T-1}
  vec exp_h                = exp(h.subvec(0, T - 1));
  vec s_t                  = sqrt(l.subvec(0, T - 1)) % exp(-0.5 * h.subvec(0, T - 1));
  vec s2_t                 = square(s_t);                 // s_t^2 = lambda_t e^{-h_t}
  vec lambda_t             = l.subvec(0, T - 1);          // lambda_t cru
  
  vec h_diff               = h.subvec(0, T - 2) - mu;
  const double diff_h1     = h[0] - mu;
  const double diff_h1_sq  = diff_h1 * diff_h1;
  
  // Somatorios (t = 1, ..., T-1)
  double sum_s2            = sum(s2_t.subvec(0, T - 2));
  double sum_s2_ylag       = sum(s2_t.subvec(0, T - 2) % y_lag1.subvec(0, T - 2));
  double sum_s2_ylag2      = sum(s2_t.subvec(0, T - 2) % square(y_lag1.subvec(0, T - 2)));
  double sum_lambda        = sum(lambda_t.subvec(0, T - 2));
  double sum_lambda_ylag   = sum(lambda_t.subvec(0, T - 2) % y_lag1.subvec(0, T - 2));
  double sum_lambda_exph   = sum(lambda_t.subvec(0, T - 2) % exp_h.subvec(0, T - 2));
  double sum_s             = sum(s_t.subvec(0, T - 2));
  double sum_s_ylag        = sum(s_t.subvec(0, T - 2) % y_lag1.subvec(0, T - 2));
  double sum_s_ylag_hdiff  = sum(s_t.subvec(0, T - 2) % y_lag1.subvec(0, T - 2) % h_diff);
  double sum_s_hdiff       = sum(s_t.subvec(0, T - 2) % h_diff);
  double sum_s_exph        = sum(s_t.subvec(0, T - 2) % exp_h.subvec(0, T - 2));
  double sum_s_exph_hdiff  = sum(s_t.subvec(0, T - 2) % exp_h.subvec(0, T - 2) % h_diff);
  double sum_hdiff         = sum(h_diff);
  double sum_hdiff2        = sum(square(h_diff));
  
  // =============================================================
  // 1. Matriz de informacao I(theta*)
  // =============================================================
  
  // Bloco (b0, .)
  I(0, 0) = s2_t[T - 1] + ratio_tau_d * sum_s2 + (1.0 / var_b0);
  I(0, 1) = (1.0 - b1_2) * (s2_t[T - 1] * y_lag1[T - 1] + ratio_tau_d * sum_s2_ylag);
  I(0, 2) = lambda_t[T - 1] + ratio_tau_d * sum_lambda;
  I(0, 3) = - (1.0 - phi) * delta / tau2 * sum_s;
  I(0, 4) = - d_b1 * delta / tau2 * sum_s_hdiff;
  
  // Bloco (b1*, .)
  I(1, 0) = I(0, 1);
  I(1, 1) = term_b1_sq * s2_t[T - 1] * (y_lag1[T - 1] * y_lag1[T - 1])
    + term_b1_sq * ratio_tau_d * sum_s2_ylag2
  + (a_b1 + b_b1) * (1.0 - b1_2);
  I(1, 2) = (1.0 - b1_2) * (y_lag1[T - 1] * lambda_t[T - 1] + ratio_tau_d * sum_lambda_ylag);
  I(1, 3) = - (1.0 - b1_2) * delta * (1.0 - phi) / tau2 * sum_s_ylag;
  I(1, 4) = - (1.0 - b1_2) * delta * d_b1 / tau2 * sum_s_ylag_hdiff;
  
  // Bloco (b2, .)
  I(2, 0) = I(0, 2);
  I(2, 1) = I(1, 2);
  I(2, 2) = lambda_t[T - 1] * exp_h[T - 1] + ratio_tau_d * sum_lambda_exph + (1.0 / var_b2);
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
  double sum_log_l_minus_l  = sum(log(l.subvec(0, T - 1)) - l.subvec(0, T - 1));
  double digamma_half_nu    = R::digamma(half_nu);
  double trigamma_half_nu   = R::trigamma(half_nu);
  
  double d_nu               = (T / 2.0) * (1.0 + std::log(half_nu)) - (T / 2.0) * digamma_half_nu + 0.5 * sum_log_l_minus_l;
  double d_prime_nu         = (T / (2.0 * nu)) - (T / 4.0) * trigamma_half_nu;
  
  double c_nu               = ((40.0 - nu) * (nu - 2.0)) / 38.0;
  double c_prime_nu         = (21.0 - nu) / 19.0;
  double term_nu_deriv      = c_prime_nu * d_nu + c_nu * d_prime_nu;
  
  I(7, 7) = - (c_nu * term_nu_deriv - (1.0 / var_nu));
  
  // =============================================================
  // 2. Derivada dI / d beta_1* (Iparb1)
  // =============================================================
  
  const double term_chain_b1   = -2.0 * b1 * (1.0 - b1_2);
  const double term_chain_b1_2 = -4.0 * b1 * term_b1_sq;
  
  double dI_b0_b1star = term_chain_b1 * (s2_t[T - 1] * y_lag1[T - 1] + ratio_tau_d * sum_s2_ylag);
  Iparb1(0, 1) = dI_b0_b1star;
  Iparb1(1, 0) = dI_b0_b1star;
  
  double dI_b1star_b1star = term_chain_b1_2 * (s2_t[T - 1] * (y_lag1[T - 1] * y_lag1[T - 1]) + ratio_tau_d * sum_s2_ylag2)
    + term_chain_b1 * (a_b1 + b_b1);
  Iparb1(1, 1) = dI_b1star_b1star;
  
  double dI_b1star_b2 = term_chain_b1 * (y_lag1[T - 1] * lambda_t[T - 1] + ratio_tau_d * sum_lambda_ylag);
  Iparb1(1, 2) = dI_b1star_b2;
  Iparb1(2, 1) = dI_b1star_b2;
  
  double dI_b1star_mu = - term_chain_b1 * (delta * (1.0 - phi) / tau2) * sum_s_ylag;
  Iparb1(1, 3) = dI_b1star_mu;
  Iparb1(3, 1) = dI_b1star_mu;
  
  double dI_b1star_phi = - term_chain_b1 * (delta * d_b1 / tau2) * sum_s_ylag_hdiff;
  Iparb1(1, 4) = dI_b1star_phi;
  Iparb1(4, 1) = dI_b1star_phi;
  
  // =============================================================
  // 3. Derivada dI / d mu (Iparmu)
  // =============================================================
  
  double dI_b0_phi_mu = (d_b1 * delta / tau2) * sum_s;
  Iparmu(0, 4) = dI_b0_phi_mu;
  Iparmu(4, 0) = dI_b0_phi_mu;
  
  double dI_b1star_phi_mu = ((1.0 - b1_2) * delta * d_b1 / tau2) * sum_s_ylag;
  Iparmu(1, 4) = dI_b1star_phi_mu;
  Iparmu(4, 1) = dI_b1star_phi_mu;
  
  double dI_b2_phi_mu = (d_b1 * delta / tau2) * sum_s_exph;
  Iparmu(2, 4) = dI_b2_phi_mu;
  Iparmu(4, 2) = dI_b2_phi_mu;
  
  double dI_mu_phi_mu = - (2.0 * phi * d_b1 / tau2) - (d_b1 * (1.0 - phi) * (T - 1.0) / tau2);
  Iparmu(3, 4) = dI_mu_phi_mu;
  Iparmu(4, 3) = dI_mu_phi_mu;
  
  double dI_mu_tau_mu = - (2.0 * d_b1 / tau2);
  Iparmu(3, 6) = dI_mu_tau_mu;
  Iparmu(6, 3) = dI_mu_tau_mu;
  
  double dI_phistar_phistar_mu = (2.0 * d_b1 * d_b1 * diff_h1 / tau2) 
    - (4.0 * phi2 * d_b1 * diff_h1 / tau2) 
    - (2.0 * d_b1 * d_b1 / tau2) * sum_hdiff;
    Iparmu(4, 4) = dI_phistar_phistar_mu;
    
    double dI_phistar_taustar_mu = - (4.0 * phi * d_b1 * diff_h1 / tau2);
    Iparmu(4, 6) = dI_phistar_taustar_mu;
    Iparmu(6, 4) = dI_phistar_taustar_mu;
    
    double dI_taustar_taustar_mu = - (4.0 * d_b1 * diff_h1 / tau2);
    Iparmu(6, 6) = dI_taustar_taustar_mu;
    
    // =============================================================
    // 4. Derivada dI / d phi* (Iparphi)
    // =============================================================
    
    double dI_b0_mu_phi = (d_b1 * delta / tau2) * sum_s;
    Iparphi(0, 3) = dI_b0_mu_phi;
    Iparphi(3, 0) = dI_b0_mu_phi;
    
    double dI_b0_phi_phi = (2.0 * phi * d_b1 * delta / tau2) * sum_s_hdiff;
    Iparphi(0, 4) = dI_b0_phi_phi;
    Iparphi(4, 0) = dI_b0_phi_phi;
    
    double dI_b1star_mu_phi = ((1.0 - b1_2) * d_b1 * delta / tau2) * sum_s_ylag;
    Iparphi(1, 3) = dI_b1star_mu_phi;
    Iparphi(3, 1) = dI_b1star_mu_phi;
    
    double dI_b1star_phi_phi = (2.0 * phi * (1.0 - b1_2) * d_b1 * delta / tau2) * sum_s_ylag_hdiff;
    Iparphi(1, 4) = dI_b1star_phi_phi;
    Iparphi(4, 1) = dI_b1star_phi_phi;
    
    double dI_b2_mu_phi = (d_b1 * delta / tau2) * sum_s_exph;
    Iparphi(2, 3) = dI_b2_mu_phi;
    Iparphi(3, 2) = dI_b2_mu_phi;
    
    double dI_b2_phi_phi = (2.0 * phi * d_b1 * delta / tau2) * sum_s_exph_hdiff;
    Iparphi(2, 4) = dI_b2_phi_phi;
    Iparphi(4, 2) = dI_b2_phi_phi;
    
    double dI_mu_mu_phi = - (2.0 * phi * d_b1 / tau2) - (2.0 * d_b1 * (1.0 - phi) * (T - 1.0) / tau2);
    Iparphi(3, 3) = dI_mu_mu_phi;
    
    double dI_mu_phi_phi = (d_b1 * diff_h1 * (2.0 - 6.0 * phi2) / tau2) 
      - (d_b1 * (1.0 + 2.0 * phi - 3.0 * phi2) / tau2) * sum_hdiff;
    Iparphi(3, 4) = dI_mu_phi_phi;
    Iparphi(4, 3) = dI_mu_phi_phi;
    
    double dI_mu_tau_phi = - (4.0 * phi * d_b1 * diff_h1 / tau2);
    Iparphi(3, 6) = dI_mu_tau_phi;
    Iparphi(6, 3) = dI_mu_tau_phi;
    
    double dI_phistar_phistar_phi = (4.0 * phi * d_b1 * (2.0 - 3.0 * phi2) * diff_h1_sq / tau2)
      - (4.0 * phi * d_b1 * d_b1 / tau2) * sum_hdiff2
    - 2.0 * phi * d_b1 * (a_phi + b_phi + 1.0);
    Iparphi(4, 4) = dI_phistar_phistar_phi;
    
    double dI_phistar_taustar_phi = (2.0 * d_b1 * (1.0 - 3.0 * phi2) * diff_h1_sq / tau2);
    Iparphi(4, 6) = dI_phistar_taustar_phi;
    Iparphi(6, 4) = dI_phistar_taustar_phi;
    
    double dI_taustar_taustar_phi = - (4.0 * phi * d_b1 * diff_h1_sq / tau2);
    Iparphi(6, 6) = dI_taustar_taustar_phi;
    
    // =============================================================
    // 5. Derivada dI / d delta (Ipardelta)
    // =============================================================
    
    const double two_delta_over_tau2 = (2.0 * delta) / tau2;
    
    double dI_b0_b0_delta = two_delta_over_tau2 * sum_s2;
    Ipardelta(0, 0) = dI_b0_b0_delta;
    
    double dI_b0_b1star_delta = two_delta_over_tau2 * (1.0 - b1_2) * sum_s2_ylag;
    Ipardelta(0, 1) = dI_b0_b1star_delta;
    Ipardelta(1, 0) = dI_b0_b1star_delta;
    
    double dI_b0_b2_delta = two_delta_over_tau2 * sum_lambda;
    Ipardelta(0, 2) = dI_b0_b2_delta;
    Ipardelta(2, 0) = dI_b0_b2_delta;
    
    double dI_b0_mu_delta = - ((1.0 - phi) / tau2) * sum_s;
    Ipardelta(0, 3) = dI_b0_mu_delta;
    Ipardelta(3, 0) = dI_b0_mu_delta;
    
    double dI_b0_phi_delta = - (d_b1 / tau2) * sum_s_hdiff;
    Ipardelta(0, 4) = dI_b0_phi_delta;
    Ipardelta(4, 0) = dI_b0_phi_delta;
    
    double dI_b1star_b1star_delta = two_delta_over_tau2 * term_b1_sq * sum_s2_ylag2;
    Ipardelta(1, 1) = dI_b1star_b1star_delta;
    
    double dI_b1star_b2_delta = two_delta_over_tau2 * (1.0 - b1_2) * sum_lambda_ylag;
    Ipardelta(1, 2) = dI_b1star_b2_delta;
    Ipardelta(2, 1) = dI_b1star_b2_delta;
    
    double dI_b1star_mu_delta = - ((1.0 - b1_2) * (1.0 - phi) / tau2) * sum_s_ylag;
    Ipardelta(1, 3) = dI_b1star_mu_delta;
    Ipardelta(3, 1) = dI_b1star_mu_delta;
    
    double dI_b1star_phi_delta = - ((1.0 - b1_2) * d_b1 / tau2) * sum_s_ylag_hdiff;
    Ipardelta(1, 4) = dI_b1star_phi_delta;
    Ipardelta(4, 1) = dI_b1star_phi_delta;
    
    double dI_b2_b2_delta = two_delta_over_tau2 * sum_lambda_exph;
    Ipardelta(2, 2) = dI_b2_b2_delta;
    
    double dI_b2_mu_delta = - ((1.0 - phi) / tau2) * sum_s_exph;
    Ipardelta(2, 3) = dI_b2_mu_delta;
    Ipardelta(3, 2) = dI_b2_mu_delta;
    
    double dI_b2_phi_delta = - (d_b1 / tau2) * sum_s_exph_hdiff;
    Ipardelta(2, 4) = dI_b2_phi_delta;
    Ipardelta(4, 2) = dI_b2_phi_delta;
    
    // =============================================================
    // 6. Derivada dI / d tau* (Ipartau)
    // =============================================================
    
    const double neg_two_delta2_over_tau2 = - (2.0 * delta2) / tau2;
    const double two_delta_over_tau2_pos  = (2.0 * delta) / tau2;
    
    double dI_b0_b0_tau = neg_two_delta2_over_tau2 * sum_s2;
    Ipartau(0, 0) = dI_b0_b0_tau;
    
    double dI_b0_b1star_tau = neg_two_delta2_over_tau2 * (1.0 - b1_2) * sum_s2_ylag;
    Ipartau(0, 1) = dI_b0_b1star_tau;
    Ipartau(1, 0) = dI_b0_b1star_tau;
    
    double dI_b0_b2_tau = neg_two_delta2_over_tau2 * sum_lambda;
    Ipartau(0, 2) = dI_b0_b2_tau;
    Ipartau(2, 0) = dI_b0_b2_tau;
    
    double dI_b0_mu_tau = two_delta_over_tau2_pos * (1.0 - phi) * sum_s;
    Ipartau(0, 3) = dI_b0_mu_tau;
    Ipartau(3, 0) = dI_b0_mu_tau;
    
    double dI_b0_phi_tau = two_delta_over_tau2_pos * d_b1 * sum_s_hdiff;
    Ipartau(0, 4) = dI_b0_phi_tau;
    Ipartau(4, 0) = dI_b0_phi_tau;
    
    double dI_b1star_b1star_tau = neg_two_delta2_over_tau2 * term_b1_sq * sum_s2_ylag2;
    Ipartau(1, 1) = dI_b1star_b1star_tau;
    
    double dI_b1star_b2_tau = neg_two_delta2_over_tau2 * (1.0 - b1_2) * sum_lambda_ylag;
    Ipartau(1, 2) = dI_b1star_b2_tau;
    Ipartau(2, 1) = dI_b1star_b2_tau;
    
    double dI_b1star_mu_tau = two_delta_over_tau2_pos * (1.0 - b1_2) * (1.0 - phi) * sum_s_ylag;
    Ipartau(1, 3) = dI_b1star_mu_tau;
    Ipartau(3, 1) = dI_b1star_mu_tau;
    
    double dI_b1star_phi_tau = two_delta_over_tau2_pos * (1.0 - b1_2) * d_b1 * sum_s_ylag_hdiff;
    Ipartau(1, 4) = dI_b1star_phi_tau;
    Ipartau(4, 1) = dI_b1star_phi_tau;
    
    double dI_b2_b2_tau = neg_two_delta2_over_tau2 * sum_lambda_exph;
    Ipartau(2, 2) = dI_b2_b2_tau;
    
    double dI_b2_mu_tau = two_delta_over_tau2_pos * (1.0 - phi) * sum_s_exph;
    Ipartau(2, 3) = dI_b2_mu_tau;
    Ipartau(3, 2) = dI_b2_mu_tau;
    
    double dI_b2_phi_tau = two_delta_over_tau2_pos * d_b1 * sum_s_exph_hdiff;
    Ipartau(2, 4) = dI_b2_phi_tau;
    Ipartau(4, 2) = dI_b2_phi_tau;
    
    double dI_mu_mu_tau = - (2.0 * d_b1 / tau2) - (2.0 * (1.0 - phi) * (1.0 - phi) * (T - 1.0) / tau2);
    Ipartau(3, 3) = dI_mu_mu_tau;
    
    double dI_mu_phi_tau = - (4.0 * phi * d_b1 * diff_h1 / tau2) - (2.0 * d_b1 * (1.0 - phi) / tau2) * sum_hdiff;
    Ipartau(3, 4) = dI_mu_phi_tau;
    Ipartau(4, 3) = dI_mu_phi_tau;
    
    double dI_mu_tau_tau = - (4.0 * d_b1 * diff_h1 / tau2);
    Ipartau(3, 6) = dI_mu_tau_tau;
    Ipartau(6, 3) = dI_mu_tau_tau;
    
    double dI_phistar_phistar_tau = (2.0 * d_b1 * d_b1 * diff_h1_sq / tau2) 
      - (4.0 * phi2 * d_b1 * diff_h1_sq / tau2) 
      - (2.0 * d_b1 * d_b1 / tau2) * sum_hdiff2;
      Ipartau(4, 4) = dI_phistar_phistar_tau;
      
      double dI_phistar_tau_tau = - (4.0 * phi * d_b1 * diff_h1_sq / tau2);
      Ipartau(4, 6) = dI_phistar_tau_tau;
      Ipartau(6, 4) = dI_phistar_tau_tau;
      
      double dI_delta_delta_tau = - (2.0 * (T - 1.0) / tau2);
      Ipartau(5, 5) = dI_delta_delta_tau;
      
      double dI_taustar_taustar_tau = - (4.0 * d_b1 * diff_h1_sq / tau2) - (8.0 * b_tau / tau2);
      Ipartau(6, 6) = dI_taustar_taustar_tau;
      
      // =============================================================
      // 7. Derivada dI / d nu* (Iparnu)
      // =============================================================
      
      double tetragamma_half_nu = R::tetragamma(half_nu); // psi^{(2)}(nu/2)
      
      // dI_nunu/dnu* = c^2 d/19 - c c'(c' d + c d') - 2 c^2 c' d' - c^3 d''
      // com d'' = -T/(2 nu^2) - (T/8) psi^{(2)}(nu/2)
      double term1_nu = (c_nu * c_nu * d_nu) / 19.0;
      double term2_nu = - c_nu * c_prime_nu * (c_prime_nu * d_nu + c_nu * d_prime_nu);
      double term3_nu = - 2.0 * c_nu * c_nu * c_prime_nu * d_prime_nu;
      double term4_nu = std::pow(c_nu, 3.0) * ((T / (2.0 * nu * nu)) + (T / 8.0) * tetragamma_half_nu);
      
      Iparnu(7, 7) = term1_nu + term2_nu + term3_nu + term4_nu;
}

void GRHMC::I_theta(const vec& theta, const SvmStats& st, 
                   mat& I, mat& Iparb1, mat& Iparmu, mat& Iparphi, 
                   mat& Ipardelta, mat& Ipartau, mat& Iparnu)
{
  I.zeros(8, 8);
  Iparb1.zeros(8, 8);
  Iparmu.zeros(8, 8);
  Iparphi.zeros(8, 8);
  Ipardelta.zeros(8, 8);
  Ipartau.zeros(8, 8);
  Iparnu.zeros(8, 8);
  
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
  I.diag() += 1e-4;
  
  // =============================================================
  // 2. dI / d beta_1*
  // =============================================================
  const double term_chain_b1   = -2.0 * b1 * (1.0 - b1_2);
  const double term_chain_b1_2 = -4.0 * b1 * term_b1_sq;
  
  Iparb1(0, 1) = Iparb1(1, 0) = term_chain_b1 * (s2T * yT + ratio_tau_d * sum_s2_ylag);
  Iparb1(1, 1) = term_chain_b1_2 * (s2T * yT * yT + ratio_tau_d * sum_s2_ylag2)
    + term_chain_b1 * (a_b1 + b_b1);
  Iparb1(1, 2) = Iparb1(2, 1) = term_chain_b1 * (yT * lamT + ratio_tau_d * sum_lambda_ylag);
  Iparb1(1, 3) = Iparb1(3, 1) = - term_chain_b1 * (delta * (1.0 - phi) / tau2) * sum_s_ylag;
  Iparb1(1, 4) = Iparb1(4, 1) = - term_chain_b1 * (delta * d_b1 / tau2) * sum_s_ylag_hdiff;
  
  // =============================================================
  // 3. dI / d mu
  // =============================================================
  Iparmu(0, 4) = Iparmu(4, 0) = (d_b1 * delta / tau2) * sum_s;
  Iparmu(1, 4) = Iparmu(4, 1) = ((1.0 - b1_2) * delta * d_b1 / tau2) * sum_s_ylag;
  Iparmu(2, 4) = Iparmu(4, 2) = (d_b1 * delta / tau2) * sum_s_exph;
  Iparmu(3, 4) = Iparmu(4, 3) = - (2.0 * phi * d_b1 / tau2) - (d_b1 * (1.0 - phi) * (T - 1.0) / tau2);
  Iparmu(3, 6) = Iparmu(6, 3) = - (2.0 * d_b1 / tau2);
  Iparmu(4, 4) = (2.0 * d_b1 * d_b1 * diff_h1 / tau2) 
    - (4.0 * phi2 * d_b1 * diff_h1 / tau2) 
    - (2.0 * d_b1 * d_b1 / tau2) * sum_hdiff;
    Iparmu(4, 6) = Iparmu(6, 4) = - (4.0 * phi * d_b1 * diff_h1 / tau2);
    Iparmu(6, 6) = - (4.0 * d_b1 * diff_h1 / tau2);
    
    // =============================================================
    // 4. dI / d phi*
    // =============================================================
    Iparphi(0, 3) = Iparphi(3, 0) = (d_b1 * delta / tau2) * sum_s;
    Iparphi(0, 4) = Iparphi(4, 0) = (2.0 * phi * d_b1 * delta / tau2) * sum_s_hdiff;
    Iparphi(1, 3) = Iparphi(3, 1) = ((1.0 - b1_2) * d_b1 * delta / tau2) * sum_s_ylag;
    Iparphi(1, 4) = Iparphi(4, 1) = (2.0 * phi * (1.0 - b1_2) * d_b1 * delta / tau2) * sum_s_ylag_hdiff;
    Iparphi(2, 3) = Iparphi(3, 2) = (d_b1 * delta / tau2) * sum_s_exph;
    Iparphi(2, 4) = Iparphi(4, 2) = (2.0 * phi * d_b1 * delta / tau2) * sum_s_exph_hdiff;
    Iparphi(3, 3) = - (2.0 * phi * d_b1 / tau2) - (2.0 * d_b1 * (1.0 - phi) * (T - 1.0) / tau2);
    Iparphi(3, 4) = Iparphi(4, 3) = (d_b1 * diff_h1 * (2.0 - 6.0 * phi2) / tau2) 
      - (d_b1 * (1.0 + 2.0 * phi - 3.0 * phi2) / tau2) * sum_hdiff;
    Iparphi(3, 6) = Iparphi(6, 3) = - (4.0 * phi * d_b1 * diff_h1 / tau2);
    Iparphi(4, 4) = (4.0 * phi * d_b1 * (2.0 - 3.0 * phi2) * diff_h1_sq / tau2)
      - (4.0 * phi * d_b1 * d_b1 / tau2) * sum_hdiff2
    - 2.0 * phi * d_b1 * (a_phi + b_phi + 1.0);
    Iparphi(4, 6) = Iparphi(6, 4) = (2.0 * d_b1 * (1.0 - 3.0 * phi2) * diff_h1_sq / tau2);
    Iparphi(6, 6) = - (4.0 * phi * d_b1 * diff_h1_sq / tau2);
    
    // =============================================================
    // 5. dI / d delta
    // =============================================================
    const double two_delta_over_tau2 = (2.0 * delta) / tau2;
    
    Ipardelta(0, 0) = two_delta_over_tau2 * sum_s2;
    Ipardelta(0, 1) = Ipardelta(1, 0) = two_delta_over_tau2 * (1.0 - b1_2) * sum_s2_ylag;
    Ipardelta(0, 2) = Ipardelta(2, 0) = two_delta_over_tau2 * sum_lambda;
    Ipardelta(0, 3) = Ipardelta(3, 0) = - ((1.0 - phi) / tau2) * sum_s;
    Ipardelta(0, 4) = Ipardelta(4, 0) = - (d_b1 / tau2) * sum_s_hdiff;
    Ipardelta(1, 1) = two_delta_over_tau2 * term_b1_sq * sum_s2_ylag2;
    Ipardelta(1, 2) = Ipardelta(2, 1) = two_delta_over_tau2 * (1.0 - b1_2) * sum_lambda_ylag;
    Ipardelta(1, 3) = Ipardelta(3, 1) = - ((1.0 - b1_2) * (1.0 - phi) / tau2) * sum_s_ylag;
    Ipardelta(1, 4) = Ipardelta(4, 1) = - ((1.0 - b1_2) * d_b1 / tau2) * sum_s_ylag_hdiff;
    Ipardelta(2, 2) = two_delta_over_tau2 * sum_lambda_exph;
    Ipardelta(2, 3) = Ipardelta(3, 2) = - ((1.0 - phi) / tau2) * sum_s_exph;
    Ipardelta(2, 4) = Ipardelta(4, 2) = - (d_b1 / tau2) * sum_s_exph_hdiff;
    
    // =============================================================
    // 6. dI / d tau*
    // =============================================================
    const double neg_two_delta2_over_tau2 = - (2.0 * delta2) / tau2;
    
    Ipartau(0, 0) = neg_two_delta2_over_tau2 * sum_s2;
    Ipartau(0, 1) = Ipartau(1, 0) = neg_two_delta2_over_tau2 * (1.0 - b1_2) * sum_s2_ylag;
    Ipartau(0, 2) = Ipartau(2, 0) = neg_two_delta2_over_tau2 * sum_lambda;
    Ipartau(0, 3) = Ipartau(3, 0) = two_delta_over_tau2 * (1.0 - phi) * sum_s;
    Ipartau(0, 4) = Ipartau(4, 0) = two_delta_over_tau2 * d_b1 * sum_s_hdiff;
    Ipartau(1, 1) = neg_two_delta2_over_tau2 * term_b1_sq * sum_s2_ylag2;
    Ipartau(1, 2) = Ipartau(2, 1) = neg_two_delta2_over_tau2 * (1.0 - b1_2) * sum_lambda_ylag;
    Ipartau(1, 3) = Ipartau(3, 1) = two_delta_over_tau2 * (1.0 - b1_2) * (1.0 - phi) * sum_s_ylag;
    Ipartau(1, 4) = Ipartau(4, 1) = two_delta_over_tau2 * (1.0 - b1_2) * d_b1 * sum_s_ylag_hdiff;
    Ipartau(2, 2) = neg_two_delta2_over_tau2 * sum_lambda_exph;
    Ipartau(2, 3) = Ipartau(3, 2) = two_delta_over_tau2 * (1.0 - phi) * sum_s_exph;
    Ipartau(2, 4) = Ipartau(4, 2) = two_delta_over_tau2 * d_b1 * sum_s_exph_hdiff;
    Ipartau(3, 3) = - (2.0 * d_b1 / tau2) - (2.0 * (1.0 - phi) * (1.0 - phi) * (T - 1.0) / tau2);
    Ipartau(3, 4) = Ipartau(4, 3) = - (4.0 * phi * d_b1 * diff_h1 / tau2) 
      - (2.0 * d_b1 * (1.0 - phi) / tau2) * sum_hdiff;
    Ipartau(3, 6) = Ipartau(6, 3) = - (4.0 * d_b1 * diff_h1 / tau2);
    Ipartau(4, 4) = (2.0 * d_b1 * d_b1 * diff_h1_sq / tau2) 
      - (4.0 * phi2 * d_b1 * diff_h1_sq / tau2) 
      - (2.0 * d_b1 * d_b1 / tau2) * sum_hdiff2;
      Ipartau(4, 6) = Ipartau(6, 4) = - (4.0 * phi * d_b1 * diff_h1_sq / tau2);
      Ipartau(5, 5) = - (2.0 * (T - 1.0) / tau2);
      Ipartau(6, 6) = - (4.0 * d_b1 * diff_h1_sq / tau2) - (8.0 * b_tau / tau2);
      
      // =============================================================
      // 7. dI / d nu*
      // =============================================================
      const double tetragamma_half_nu = R::tetragamma(half_nu);
      
      const double term1_nu = (c_nu * c_nu * d_nu) / 19.0;
      const double term2_nu = - c_nu * c_prime_nu * (c_prime_nu * d_nu + c_nu * d_prime_nu);
      const double term3_nu = - 2.0 * c_nu * c_nu * c_prime_nu * d_prime_nu;
      const double term4_nu = std::pow(c_nu, 3.0) * ((T / (2.0 * nu * nu)) + (T / 8.0) * tetragamma_half_nu);
      
      Iparnu(7, 7) = term1_nu + term2_nu + term3_nu + term4_nu;
}

vec GRHMC::updatetheta(const vec& theta_cur, const vec& y, const vec& h_cur, 
                       const vec& l_cur, const int T, const int nfix,
                       const int nlpfrog, const double epsilon, const double m, 
                       const double c, int& acc)
{
  // y = (y_0, y_1, ..., y_T): T+1 elementos; h e l: T elementos
  if (y.n_elem != (uword)(T + 1))
    Rcpp::stop("y deve ter T+1 elementos: (y0, y1, ..., yT).");
  if (h_cur.n_elem < (uword)T || l_cur.n_elem < (uword)T)
    Rcpp::stop("h e lambda devem ter pelo menos T elementos.");
  
  // Estatisticas O(T), calculadas uma unica vez
  const SvmStats st = make_stats(y, h_cur, l_cur, T);
  
  mat I(8, 8, fill::zeros);
  mat Iparb1(8, 8, fill::zeros);
  mat Iparmu(8, 8, fill::zeros);
  mat Iparphi(8, 8, fill::zeros);
  mat Ipardelta(8, 8, fill::zeros);
  mat Ipartau(8, 8, fill::zeros);
  mat Iparnu(8, 8, fill::zeros);
  
  I_theta(theta_cur, st, I, Iparb1, Iparmu, Iparphi, Ipardelta, Ipartau, Iparnu);
  mat Lchol;
  if(!I.is_finite() || !chol(Lchol, I, "lower")) return theta_cur;
  vec pcur = Lchol * sample_pI(m, c, 8);
  
  vec u;
  if(!solve(u, trimatl(Lchol), pcur)) return theta_cur;
  double M = m * std::sqrt(dot(u, u) / (m * m * c * c) + 1.0);
  double logdetG = 2.0 * sum(log(Lchol.diag()));
  
  double H_init = -lpost_theta(theta_cur, y, h_cur, l_cur, T);
  H_init += 0.5 * logdetG;
  H_init += c * c * M;
  
  vec thetaa = theta_cur, thetab = theta_cur, pa = pcur, pb = pcur;
  vec gdpth1(8), gradHMCth1(8), gdmom1(8), gdmom2(8);
  
  // Gradiente em thetaa (reaproveitado entre passos de leapfrog)
  gdpth1 = grad_theta(thetaa, y, h_cur, l_cur, T, Lchol, Iparb1, Iparmu, 
                      Iparphi, Ipardelta, Ipartau, Iparnu);
  
  for(int jgl = 1; jgl <= nlpfrog; ++jgl) 
  {
    // dH/dtheta = -grad_theta - 0.5 * nuH   (Lchol, Ipar* e gdpth1 ja em thetaa)
    gradHMCth1 = - gdpth1 - 0.5 * nuHtheta(Lchol, Iparb1, Iparmu, Iparphi, 
                                           Ipardelta, Ipartau, Iparnu, pa, m,c);
    pb = pa;
    
    // Meio passo no momento (ponto fixo)
    for(int jfix = 1; jfix <= nfix; ++jfix)
    {
      pb = pa - 0.5 * epsilon * gradHMCth1;
      
      gradHMCth1 = - gdpth1 - 0.5 * nuHtheta(Lchol, Iparb1, Iparmu, Iparphi, 
                                             Ipardelta, Ipartau, Iparnu, pb, m, c);
    }
    
    // Passo completo na posicao (ponto fixo)
    gdmom1 = gradHmomtheta(Lchol, pb, m, c);
    gdmom2 = gdmom1;
    thetab = thetaa;
    
    for(int jfix = 1; jfix <= nfix; ++jfix) 
    {
      thetab = thetaa + 0.5 * epsilon * gdmom1 + 0.5 * epsilon * gdmom2;
      if(!thetab.is_finite()) return theta_cur;
      
      I_theta(thetab, st, I, Iparb1, Iparmu, Iparphi, Ipardelta, Ipartau, Iparnu);
      if(!I.is_finite() || !chol(Lchol, I, "lower")) return theta_cur;
      
      gdmom2 = gradHmomtheta(Lchol, pb, m, c);
    }
    
    // Meio passo final no momento, em thetab
    gdpth1 = grad_theta(thetab, y, h_cur, l_cur, T, Lchol, Iparb1, Iparmu, 
                        Iparphi, Ipardelta, Ipartau, Iparnu);
    
    gradHMCth1 = - gdpth1 - 0.5 * nuHtheta(Lchol, Iparb1, Iparmu, Iparphi, 
                                           Ipardelta, Ipartau, Iparnu, pb, m, c);
    
    pb = pb - 0.5 * epsilon * gradHMCth1;
    
    if(!pb.is_finite()) return theta_cur;
    
    pa = pb;
    thetaa = thetab;
  }
  
  if(!solve(u, trimatl(Lchol), pb)) return theta_cur;
  M = m * std::sqrt(dot(u, u) / (m * m * c * c) + 1.0);
  logdetG = 2.0 * sum(log(Lchol.diag()));
  double H_final = -lpost_theta(thetab, y, h_cur, l_cur, T);
  H_final += 0.5 * logdetG;
  H_final += c * c * M;
  
  const double log_alpha = -H_final + H_init;
  
  if(!std::isfinite(log_alpha)) return theta_cur;
  
  const double aMH = std::min(0.0, log_alpha);
  const double uMH = std::log(randu());
  
  if(uMH < aMH){
    acc += 1;
    return thetab;
  } else {
    return theta_cur;
  }
}

// #############################################################################
// ############################ h ##############################################
// #############################################################################
double GRHMC::lpost_h(const vec& h, const vec& theta, const vec& y, const vec& l, 
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

vec GRHMC::grad_h(const vec& h, const vec& theta, const vec& y, const vec& l, 
                  const int T, const vec& ldiag, const vec& msubdiag)
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
  
  // 2. Cálculo do vetor c_t = b2 * sqrt(lambda_t) * exp(h_t / 2)
  vec ct = b2 * sqrt(l_sub) % exp(0.5 * h);
  
  // 3. Vetores para armazenar a diagonal e superdiagonal de I^{-1}
  vec omega_diag(T);
  vec omega_super(T - 1);
  
  // --- Recorrência para trás (Backward Recurrence) ---
  // t = T (índice T-1)
  omega_diag[T - 1] = 1.0 / (ldiag[T - 1] * ldiag[T - 1]);
  
  for(int t = T - 2; t >= 0; --t) 
  {
    const double d_t   = ldiag[t];
    const double d_next= ldiag[t + 1];
    const double l_t   = msubdiag[t];
    
    //omega_super[t] = - (l_t / d_next) * omega_diag[t + 1];
    omega_super[t] = - (l_t / d_t) * omega_diag[t + 1];
    omega_diag[t]  = (1.0 / (d_t * d_t)) + ((l_t / d_t) * (l_t / d_t)) * omega_diag[t + 1];
  }
  
  // 4. Montagem do vetor v
  // vec v(T);
  
  // t = 1..T-1 (índices 0..T-2)
  for(int t = 0; t < T - 1; ++t) 
  {
    const double c_t  = ct[t];
    const double term1 = c_t * c_t + (delta * c_t * (delta * c_t - phi)) / tau2;
    const double term2 = (delta * c_t) / tau2;
    
    //v[t] = term1 * omega_diag[t] + term2 * omega_super[t];
    grad[t] += - 0.5*( term1 * omega_diag[t] + term2 * omega_super[t] );
  }
  
  // t = T (índice T-1)
  //v[T-1] = ct[T - 1] * ct[T - 1] * omega_diag[T - 1];
  grad[T-1] += - 0.5*( ct[T-1] * ct[T-1] * omega_diag[T-1] );
  
  return grad;
}

vec GRHMC::nuHh(const vec& h, const vec& theta, const vec& l, const int T,
                const vec& ldiag, const vec& msubdiag, const vec& p, const double m, 
                const double c)
{
  // 1. Unpack dos parâmetros
  const double b2    = theta[2];
  const double phi   = std::tanh(theta[4]);
  const double delta = theta[5];
  const double tau   = std::exp(theta[6]);
  const double tau2  = tau * tau;
  
  // 2. Resolução do sistema I x = p (x = I^{-1} p) em O(T) via Cholesky bidiagonal
  vec u(T);
  // Forward substitution: L u = p
  u[0] = p[0] / ldiag[0];
  for(int i = 1; i < T; ++i) 
  {
    u[i] = (p[i] - msubdiag[i - 1] * u[i - 1]) / ldiag[i];
  }
  double M = m * std::sqrt(dot(u, u) / (m * m * c * c) + 1.0);
  
  vec x(T);
  // Backward substitution: L^T x = u
  x[T - 1] = u[T - 1] / ldiag[T - 1];
  for(int i = T - 2; i >= 0; --i) 
  {
    x[i] = (u[i] - msubdiag[i] * x[i + 1]) / ldiag[i];
  }
  
  // 3. Vetor c_t = b2 * sqrt(lambda_t) * exp(0.5 * h_t)
  const vec l_sub = l.subvec(0, T - 1);
  vec ct = b2 * sqrt(l_sub) % exp(0.5 * h);
  
  // 4. Montagem de nuH
  vec nuH(T);
  
  // t = 1..T-1 (índices 0..T-2)
  for(int t = 0; t < T - 1; ++t) 
  {
    const double c_t   = ct[t];
    const double term1 = c_t * c_t + (delta * c_t * (delta * c_t - phi)) / tau2;
    const double term2 = (delta * c_t) / tau2;
    
    const double x_t    = x[t];
    const double x_next = x[t + 1];
    
    nuH[t] = term1 * (x_t * x_t) + term2 * (x_t * x_next);
  }
  
  // t = T (índice T-1)
  nuH[T - 1] = (ct[T - 1] * ct[T - 1]) * (x[T - 1] * x[T - 1]);
  
  
  return (1/M) * nuH;
}

vec GRHMC::gradHmomh(const vec& ldiag, const vec& msubdiag, const vec& p, 
                     const int T, const double m, 
                     const double c)
{
  vec u(T);
  // 1. Forward substitution: L u = p
  u[0] = p[0] / ldiag[0];
  for(int i = 1; i < T; ++i) 
  {
    u[i] = (p[i] - msubdiag[i - 1] * u[i - 1]) / ldiag[i];
  }
  double M = m * std::sqrt(dot(u, u) / (m * m * c * c) + 1.0);
  
  vec x(T);
  // 2. Backward substitution: L^T x = u
  x[T - 1] = u[T - 1] / ldiag[T - 1];
  for(int i = T - 2; i >= 0; --i) 
  {
    x[i] = (u[i] - msubdiag[i] * x[i + 1]) / ldiag[i];
  }
  
  return (1/M) * x;
}

bool GRHMC::I_h(const vec& h, const vec& theta, const vec& l, const int T,
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

vec GRHMC::updateh(const vec& h_cur, const vec& theta_cur, const vec& y, 
                   const vec& l_cur, const int T, const int nfix, 
                   const int nlpfrog, const double epsilon, const double m, 
                   const double c, int& acc)
{
  // Checagens de tamanho de vetor
  if(y.n_elem != (uword)(T + 1))
    Rcpp::stop("y deve ter T+1 elementos: (y0, y1, ..., yT).");
  if(h_cur.n_elem < (uword)T || l_cur.n_elem < (uword)T)
    Rcpp::stop("h e lambda devem ter pelo menos T elementos.");
  
  // 1. Matriz de Informacao de Fisher inicial I(h_cur)
  vec ldiag_I, msubdiag_I, ldiag, msubdiag;
  if(!I_h(h_cur, theta_cur, l_cur, T, ldiag_I, msubdiag_I, ldiag, msubdiag)){
    return h_cur; // Falha na Cholesky (nao definida positiva)
  }
  
  // 2. Amostragem do momento inicial p ~ N(0, I(h_cur)) -> p = L * eta
  //vec eta = arma::randn<vec>(T);
  vec eta = sample_pI(m, c, T);
  vec pcur(T);
  apply_chol_tri(ldiag, msubdiag, eta, pcur, T);
  
  // u = L^{-1} * pcur = eta
  vec u = eta;
  double M = m * std::sqrt(dot(u, u) / (m * m * c * c) + 1.0);
  // 3. Log-determinante de I(h_cur) e Hamiltoniano inicial H_init
  double logdetG = 2.0 * sum(log(ldiag));
  double H_init = -lpost_h(h_cur, theta_cur, y, l_cur, T);
  H_init += 0.5 * logdetG;
  H_init += c * c * M;
  
  // 4. Copia dos estados para os iteradores do integrador de ponto fixo
  vec ha = h_cur, hb = h_cur, pa = pcur, pb = pcur;
  vec gdph1(T), gradHMCh1(T), gdmom1(T), gdmom2(T);
  
  // Gradiente inicial d(-lpost_h)/dh + 0.5 * tr[I^{-1} dI/dh]
  gdph1 = grad_h(ha, theta_cur, y, l_cur, T, ldiag, msubdiag);
  
  // 5. Integrador de Ponto Fixo (Generalized Leapfrog)
  for(int jgl = 1; jgl <= nlpfrog; ++jgl)
  {
    // dH/dh = -grad_h + 0.5 * nuHh
    gradHMCh1 = - gdph1 - 0.5 * nuHh(ha, theta_cur, l_cur, T, ldiag, msubdiag, 
                                     pa, m, c);
    pb = pa;
    
    // Meio passo no momento (ponto fixo)
    for(int jfix = 1; jfix <= nfix; ++jfix)
    {
      pb = pa - 0.5 * epsilon * gradHMCh1;
      gradHMCh1 = - gdph1 - 0.5 * nuHh(ha, theta_cur, l_cur, T, ldiag, msubdiag, 
                                       pb, m, c);
    }
    
    // Passo completo na posicao h (ponto fixo)
    gdmom1 = gradHmomh(ldiag, msubdiag, pb, T, m, c);
    gdmom2 = gdmom1;
    hb = ha;
    
    for(int jfix = 1; jfix <= nfix; ++jfix)
    {
      hb = ha + 0.5 * epsilon * gdmom1 + 0.5 * epsilon * gdmom2;
      
      if (!hb.is_finite()) return h_cur;
      
      // Recomputa I(hb) e a Cholesky bidiagonal
      if (!I_h(hb, theta_cur, l_cur, T, ldiag_I, msubdiag_I, ldiag, msubdiag)) {
        return h_cur;
      }
      
      gdmom2 = gradHmomh(ldiag, msubdiag, pb, T, m, c);
    }
    
    // Meio passo final no momento em hb
    gdph1 = grad_h(hb, theta_cur, y, l_cur, T, ldiag, msubdiag);
    gradHMCh1 = - gdph1 - 0.5 * nuHh(hb, theta_cur, l_cur, T, ldiag, msubdiag, 
                                     pb, m, c);
    pb = pb - 0.5 * epsilon * gradHMCh1;
    
    if (!pb.is_finite()) return h_cur;
    
    pa = pb;
    ha = hb;
  }
  
  // 6. Hamiltoniano final H_final
  // Resolvendo L * u_final = pb
  vec u_final(T);
  u_final[0] = pb[0] / ldiag[0];
  for (int i = 1; i < T; ++i) {
    u_final[i] = (pb[i] - msubdiag[i - 1] * u_final[i - 1]) / ldiag[i];
  }
  M = m * std::sqrt(dot(u_final, u_final) / (m * m * c * c) + 1.0);
  logdetG = 2.0 * sum(log(ldiag));
  double H_final = -lpost_h(hb, theta_cur, y, l_cur, T);
  H_final += 0.5 * logdetG;
  H_final += c * c * M;
  
  // 7. Teste de Aceitacao de Metropolis-Hastings
  const double log_alpha = -H_final + H_init;
  
  if (!std::isfinite(log_alpha)) return h_cur;
  
  const double aMH = std::min(0.0, log_alpha);
  const double uMH = std::log(randu());
  
  if (uMH < aMH) {
    acc += 1;
    return hb;
  } else {
    return h_cur;
  }
  
}


// #############################################################################
// ############################ l ##############################################
// #############################################################################

double GRHMC::lpost_l(const vec& l_star, const vec& h, const vec& theta, 
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

vec GRHMC::grad_l(const vec& l_star, const vec& h, const vec& theta, const vec& y, 
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
  grad.subvec(0, T - 2) += (0.5 * delta / tau2) * zt.subvec(0, T - 2) % At;
  
  // -0.5 * tr(I^{-1} dI/dl*_t) = -0.5 * (nu*lambda_t/2) / I_tt  (I diagonal)
  //const vec I_diag = I_l(l, theta, T);
  const vec I_diag = I_l(l_star, theta, T);
  return grad - 0.5 * (0.5 * nu * l_sub) / I_diag;
}

vec GRHMC::nuHl(const vec& l_star, const vec& theta, const vec& p, const int T, 
                const double m, const double c)
{
  // 1. Unpack dos parâmetros de theta
  const double delta       = theta[5];
  const double tau         = std::exp(theta[6]);
  const double logistic_nu = 1.0 / (1.0 + std::exp(-theta[7]));
  const double nu          = 40.0 * logistic_nu + 2.0 * (1.0 - logistic_nu);
  
  const double tau2  = tau * tau;
  const double k_val = (0.25 * delta * delta) / tau2;
  
  // 2. Sub-vetores de lambda e do momento p
  const vec l = exp(l_star);
  const vec l_sub = l.subvec(0, T - 1);
  const vec p_sub = p.subvec(0, T - 1);
  
  // 3. Montagem do vetor k_T (k_t = delta^2 / (4*tau^2) para t=1..T-1; k_T = 0)
  vec k(T, fill::zeros);
  k.subvec(0, T - 2).fill(k_val);
  
  // 4. Numerador: (nu / 2) * lambda_t * p_t^2
  vec num = (0.5 * nu * l_sub) % square(p_sub);
  
  // 5. Denominador: I_{tt}^2 = (0.5 + (nu / 2) * lambda_t + k_t)^2
  vec I_diag = 0.5 + 0.5 * nu * l_sub + k;
  vec den    = square(I_diag);
  
  // 6. nuH_t = num_t / den_t
  vec nuH = num / den;
  
  double M = m * std::sqrt(dot(p, (1/I_diag) % p)/(m*m * c*c) +1);
  
  return (1/M) * nuH;
}

vec GRHMC::gradHmoml(const vec& l_star, const vec& theta, const vec& p, const int T, 
                     const double m, const double c)
{
  // 1. Unpack dos parâmetros de theta
  const double delta       = theta[5];
  const double tau         = std::exp(theta[6]);
  const double logistic_nu = 1.0 / (1.0 + std::exp(-theta[7]));
  const double nu          = 40.0 * logistic_nu + 2.0 * (1.0 - logistic_nu);
  
  const double tau2  = tau * tau;
  const double k_val = (0.25 * delta * delta) / tau2;
  
  // 2. Sub-vetores de lambda e do momento p
  const vec l = exp(l_star);
  const vec l_sub = l.subvec(0, T - 1);
  const vec p_sub = p.subvec(0, T - 1);
  
  // 3. Montagem da diagonal de I(lambda*_T): 0.5 + (nu/2) * lambda_t + k_t
  vec I_diag(T);
  I_diag = 0.5 + 0.5 * nu * l_sub;
  I_diag.subvec(0, T - 2) += k_val; // Adiciona k_val para t = 1..T-1 (k_T = 0)
  
  // 4. Velocidade: I(lambda*)^(-1) * p (elemento a elemento)
  vec v_l = p_sub / I_diag;
  
  double M = m * std::sqrt(dot(p, (1/I_diag) % p)/(m*m * c*c) +1);
  
  return (1/M) * v_l;
}

vec GRHMC::I_l(const vec& l_star, const vec& theta, const int T)
{
  // 1. Unpack dos parâmetros de theta (com reparametrização logística para nu)
  const double delta       = theta[5];
  const double tau         = std::exp(theta[6]);
  const double logistic_nu = 1.0 / (1.0 + std::exp(-theta[7]));
  const double nu          = 40.0 * logistic_nu + 2.0 * (1.0 - logistic_nu);
  
  const double tau2  = tau * tau;
  const double k_val = (0.25 * delta * delta) / tau2;
  
  // 2. Sub-vetor de lambda (t = 1..T)
  const vec l = exp(l_star);
  const vec l_sub = l.subvec(0, T - 1);
  
  // 3. Diagonal base: 0.5 + (nu / 2) * lambda_t
  vec I_diag = 0.5 + 0.5 * nu * l_sub;
  
  // 4. Adiciona k_t = delta^2 / (4 * tau^2) para t = 1..T-1 (índices 0..T-2); k_T = 0
  I_diag.subvec(0, T - 2) += k_val;
  
  return I_diag;
}

vec GRHMC::updatel(const vec& lstar_cur, const vec& h_cur, const vec& theta_cur, 
                   const vec& y, const int T, const int nfix, 
                   const int nlpfrog, const double epsilon, 
                   const double m, const double c, int& acc)
{
  if (y.n_elem != (uword)(T + 1))
    Rcpp::stop("y deve ter T+1 elementos: (y0, y1, ..., yT).");
  if (h_cur.n_elem != (uword)T || lstar_cur.n_elem != (uword)T || T < 2)
    Rcpp::stop("h e lambda* devem ter T elementos, e T >= 2.");
  
  // Diagonal da metrica em lstar_cur (>= 1/2 se finita)
  vec I_diag = I_l(lstar_cur, theta_cur, T);
  if (!I_diag.is_finite()) return lstar_cur;
  
  vec pcur = sqrt(I_diag) % sample_pI(m, c, T);
  double M = m * std::sqrt(dot(pcur, (1/I_diag) % pcur)/(m*m * c*c) +1);
  double logdetG = sum(log(I_diag));
  
  double H_init  = -lpost_l(lstar_cur, h_cur, theta_cur, y, T);
  H_init += 0.5 * logdetG;
  H_init += c * c * M;
  
  vec la = lstar_cur, lb = lstar_cur, pa = pcur, pb = pcur;
  vec gdpl1(T), gradHMCl1(T), gdmom1(T), gdmom2(T);
  
  gdpl1 = grad_l(la, h_cur, theta_cur, y, T);
  
  for (int jgl = 1; jgl <= nlpfrog; ++jgl)
  {
    // dH/dl* = -grad_l - 0.5 * nuHl
    gradHMCl1 = - gdpl1 - 0.5 * nuHl(la, theta_cur, pa, T, m, c);
    pb = pa;
    
    // Meio passo no momento (ponto fixo)
    for (int jfix = 1; jfix <= nfix; ++jfix)
    {
      pb = pa - 0.5 * epsilon * gradHMCl1;
      gradHMCl1 = - gdpl1 - 0.5 * nuHl(la, theta_cur, pb, T, m, c);
    }
    
    // Passo completo na posicao l* (ponto fixo)
    gdmom1 = gradHmoml(la, theta_cur, pb, T, m, c);
    gdmom2 = gdmom1;
    lb = la;
    
    for (int jfix = 1; jfix <= nfix; ++jfix)
    {
      lb = la + 0.5 * epsilon * gdmom1 + 0.5 * epsilon * gdmom2;
      if (!lb.is_finite()) return lstar_cur;
      
      I_diag = I_l(lb, theta_cur, T);          // exp(lb) pode estourar
      if (!I_diag.is_finite()) return lstar_cur;
      
      gdmom2 = gradHmoml(lb, theta_cur, pb, T, m, c);
    }
    
    // Meio passo final no momento, em lb
    gdpl1 = grad_l(lb, h_cur, theta_cur, y, T);
    gradHMCl1 = - gdpl1 - 0.5 * nuHl(lb, theta_cur, pb, T, m, c);
    pb = pb - 0.5 * epsilon * gradHMCl1;
    
    if (!pb.is_finite() || !gdpl1.is_finite()) return lstar_cur;
    
    pa = pb;
    la = lb;
  }
  
  // Hamiltoniano final
  I_diag = I_l(lb, theta_cur, T);
  //u = pb / sqrt(I_diag);
  M = m * std::sqrt(dot(pb, (1/I_diag) % pb)/(m*m * c*c) +1);
  logdetG = sum(log(I_diag));
  
  double H_final = -lpost_l(lb, h_cur, theta_cur, y, T);
  H_final +=  + 0.5 * logdetG;
  H_final += c * c * M;
  
  const double log_alpha = -H_final + H_init;
  if (!std::isfinite(log_alpha)) return lstar_cur;
  
  if (std::log(randu()) < std::min(0.0, log_alpha)) {
    acc += 1;
    return lb;               // em escala lambda*
  } else {
    return lstar_cur;
  }
}


