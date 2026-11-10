#include "rmlmc.h"
#include <cmath>

// Parameters priori
double mu_b0 = 0.0, s_b0 = 3.2, sig2_b0 = s_b0 * s_b0;
double a_b1 = 5.0, b_b1 = 1.5;
double mu_b2 = 0.0, s_b2 = 3.2, sig2_b2 = s_b2 * s_b2;
double mu_0 = 0.0, s_0 = 3.2, sig2_mu = s_0 * s_0;
double a_phi = 20.0, b_phi = 1.5;
double mu_delta = 0.0, s_delta = 0.1, sig2_delta = s_delta * s_delta;
double a_tau = 2.5, b_tau = 0.02;      
double mu_nu = -1.0, s_nu = 1.0, sig2_nu = s_nu * s_nu;

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

SvmStats RMLMC::make_stats(const vec& y, const vec& h, const vec& l, const int T)
{
  SvmStats st;
  st.T  = T;
  st.h1 = h[0];
  
  // vetores completos (y tem T+1 elementos: y_0..y_T)
  st.yt   = y.subvec(1, T);          // y_1..y_T
  st.ylag = y.subvec(0, T - 1);      // y_0..y_{T-1}
  st.h    = h.subvec(0, T - 1);      // h_1..h_T
  st.eh   = exp(st.h);
  const vec ll = l.subvec(0, T - 1);
  st.s    = sqrt(ll) % exp(-0.5 * st.h);
  const vec s2 = square(st.s);
  
  // termos t = T
  st.s2_T     = s2[T - 1];
  st.ylag_T   = st.ylag[T - 1];
  st.lam_T    = ll[T - 1];
  st.lamexp_T = ll[T - 1] * st.eh[T - 1];
  
  // somas em t = 1..T-1 (indices 0..T-2)
  const vec sa  = st.s.head(T - 1);
  const vec s2a = s2.head(T - 1);
  const vec la  = ll.head(T - 1);
  const vec ya  = st.ylag.head(T - 1);
  const vec eha = st.eh.head(T - 1);
  const vec ha  = st.h.head(T - 1);
  
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



// #############################################################################
// ############################ theta ##########################################
// #############################################################################
double RMLMC::lpost_theta(const vec& theta, const vec& y, const vec& h, 
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

vec RMLMC::grad_theta(const vec& theta, const vec& y, const vec& h, const vec& l,
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


void RMLMC::H_theta(const vec& theta, const SvmStats& st, mat& H, const int T)
{
  H.zeros(8, 8);
  const double Td = static_cast<double>(T);
  const uword  n  = T - 1;                       // t = 1..T-1
  
  // 1. Unpack dos parametros
  const double b0          = theta[0];
  const double b1          = std::tanh(theta[1]);
  const double b2          = theta[2];
  const double mu          = theta[3];
  const double phi         = std::tanh(theta[4]);
  const double delta       = theta[5];
  const double tau         = std::exp(theta[6]);
  const double logistic_nu = 1.0 / (1.0 + std::exp(-theta[7]));
  const double nu          = 40.0 * logistic_nu + 2.0 * (1.0 - logistic_nu);
  
  // 2. Termos auxiliares
  const double phi2       = phi * phi;
  const double tau2       = tau * tau;
  const double delta2     = delta * delta;
  const double b1_2       = b1 * b1;
  const double omb1       = 1.0 - b1_2;            // 1 - beta1^2
  const double omb1_sq    = omb1 * omb1;
  const double omphi      = 1.0 - phi;             // 1 - phi
  const double omphi2     = 1.0 - phi2;            // 1 - phi^2
  const double ratio_tau_d = (tau2 + delta2) / tau2;
  
  const double diff_h1    = st.h1 - mu;
  const double diff_h1_sq = diff_h1 * diff_h1;
  
  // termos t = T
  const double s2T     = st.s2_T;
  const double yT      = st.ylag_T;
  const double lamT    = st.lam_T;
  const double lamexpT = st.lamexp_T;
  
  // somas que nao dependem de theta (t = 1..T-1)
  const double sum_s2          = st.S_s2;
  const double sum_s2_ylag     = st.S_s2y;
  const double sum_s2_ylag2    = st.S_s2y2;
  const double sum_lambda      = st.S_lam;
  const double sum_lambda_ylag = st.S_lamy;
  const double sum_lambda_exph = st.S_lamexp;
  const double sum_s           = st.S_s;
  const double sum_s_ylag      = st.S_sy;
  const double sum_s_exph      = st.S_sexp;
  
  // somas que dependem de mu (expansao polinomial)
  const double sum_s_hdiff      = st.S_sh    - mu * st.S_s;
  const double sum_s_ylag_hdiff = st.S_syh   - mu * st.S_sy;
  const double sum_s_exph_hdiff = st.S_sexph - mu * st.S_sexp;
  const double sum_hdiff        = st.S_h     - mu * (Td - 1.0);
  const double sum_hdiff2       = st.S_h2 - 2.0 * mu * st.S_h + mu * mu * (Td - 1.0);
  
  // 3. z_t (t = 1..T) e A_t (t = 1..T-1)
  const vec z  = st.s % (st.yt - b0 - b1 * st.ylag - b2 * st.eh);
  const vec zc = z.head(n);
  const vec hc = st.h.head(n);
  const vec hd = hc - mu;
  const vec A  = st.h.subvec(1, n) - mu - phi * hd - delta * zc;
  const vec Ad = A - delta * zc;                    // A_t - delta z_t
  
  const vec sc  = st.s.head(n);
  const vec ylc = st.ylag.head(n);
  const vec ehc = st.eh.head(n);
  
  const double S_zsy_T = accu(z % st.s % st.ylag);   // sum_{t=1}^T z s y_{t-1}
  const double S_sA    = dot(sc, A);
  const double S_sAd   = dot(sc, Ad);
  const double S_syA   = accu(sc % ylc % A);
  const double S_syAd  = accu(sc % ylc % Ad);
  const double S_seA   = accu(sc % ehc % A);
  const double S_seAd  = accu(sc % ehc % Ad);
  const double S_A     = accu(A);
  const double S_Ahd   = dot(A, hd);
  const double S_z     = accu(zc);
  const double S_zhd   = dot(zc, hd);
  const double S_z2    = dot(zc, zc);
  const double S_Az    = dot(A, zc);
  const double S_A2    = dot(A, A);
  
  // =============================================================
  // 4. Matriz Hessiana H(theta*)  (indices: 0 b0, 1 b1*, 2 b2, 3 mu,
  //    4 phi*, 5 delta, 6 tau*, 7 nu*)
  // =============================================================
  // beta0
  H(0,0) = -s2T - ratio_tau_d * sum_s2 - 1.0 / sig2_b0;
  H(0,1) = omb1 * (-s2T * yT - ratio_tau_d * sum_s2_ylag);
  H(0,2) = -lamT - ratio_tau_d * sum_lambda;
  H(0,3) = omphi * delta / tau2 * sum_s;
  H(0,4) = omphi2 * delta / tau2 * sum_s_hdiff;
  H(0,5) = -S_sAd / tau2;
  H(0,6) = 2.0 * delta / tau2 * S_sA;
  
  // beta1*
  H(1,1) = -2.0 * b1 * omb1 * S_zsy_T
  - omb1_sq * (sum_s2_ylag2 + s2T * yT * yT)
    + 2.0 * b1 * omb1 * delta / tau2 * S_syA
    - omb1_sq * delta2 / tau2 * sum_s2_ylag2
    - (a_b1 + b_b1) * omb1;
    H(1,2) = -omb1 * (yT * lamT + ratio_tau_d * sum_lambda_ylag);
    H(1,3) = omb1 * delta * omphi / tau2 * sum_s_ylag;
    H(1,4) = omb1 * delta * omphi2 / tau2 * sum_s_ylag_hdiff;
    H(1,5) = -omb1 / tau2 * S_syAd;
    H(1,6) = 2.0 * delta * omb1 / tau2 * S_syA;
    
    // beta2
    H(2,2) = -lamexpT - ratio_tau_d * sum_lambda_exph - 1.0 / sig2_b2;
    H(2,3) = omphi * delta / tau2 * sum_s_exph;
    H(2,4) = omphi2 * delta / tau2 * sum_s_exph_hdiff;
    H(2,5) = -S_seAd / tau2;
    H(2,6) = 2.0 * delta / tau2 * S_seA;
    
    // mu
    H(3,3) = -omphi2 / tau2 - omphi * omphi * (Td - 1.0) / tau2 - 1.0 / sig2_mu;
    H(3,4) = -2.0 * phi * omphi2 * diff_h1 / tau2
    - omphi2 / tau2 * S_A
    - omphi2 * omphi / tau2 * sum_hdiff;
    H(3,5) = -omphi / tau2 * S_z;
    H(3,6) = -2.0 * omphi2 * diff_h1 / tau2 - 2.0 * omphi / tau2 * S_A;
    
    // phi*
    H(4,4) = diff_h1_sq * omphi2 * omphi2 / tau2
    - 2.0 * phi2 * omphi2 * diff_h1_sq / tau2
    - 2.0 * phi * omphi2 / tau2 * S_Ahd
    - omphi2 * omphi2 / tau2 * sum_hdiff2
    - omphi2 * (a_phi + b_phi + 1.0);
    H(4,5) = -omphi2 / tau2 * S_zhd;
    H(4,6) = -2.0 * phi * omphi2 * diff_h1_sq / tau2
    - 2.0 * omphi2 / tau2 * S_Ahd;
    
    // delta, tau*
    H(5,5) = -S_z2 / tau2 - 1.0 / sig2_delta;
    H(5,6) = -2.0 / tau2 * S_Az;
    H(6,6) = -2.0 / tau2 * (omphi2 * diff_h1_sq + S_A2 + 2.0 * b_tau);
    
    // nu* (bloco isolado)
    const double J  = (40.0 - nu) * (nu - 2.0) / 38.0;     // d nu / d nu*
    const double c2 = (21.0 - nu) / 19.0;                  // (1 - 2L), d2nu/dnu*2 = J*c2
    const double d_nu = 0.5 * Td * (1.0 + std::log(0.5 * nu))
      - 0.5 * Td * R::digamma(0.5 * nu)
      + 0.5 * st.sum_logl_minus_l;
      const double dd_nu = Td / (2.0 * nu) - 0.25 * Td * R::trigamma(0.5 * nu);
      H(7,7) = J * (c2 * d_nu + J * dd_nu) - 1.0 / sig2_nu;
      
      // 5. Simetriza (preencheu so o triangulo superior)
      H = symmatu(H);
}


void RMLMC::sqrtG_theta(const vec& theta, const vec& y, const vec& h, const vec& l,
                        const int T, const double alpha, mat& sqrtG)
{
  
  vec gradtheta = grad_theta(theta, y, h, l, T);
  double norm2_grad = dot(gradtheta, gradtheta); 
  double Lalpha = 1.0 + alpha * alpha * norm2_grad;
  
  sqrtG = eye<mat>(8, 8); 
  
  const double tol = 1e-8;
  
  if (norm2_grad < tol) {
    sqrtG -= (0.5 * alpha * alpha) * (gradtheta * gradtheta.t());
  } else {
    double scale = (1.0 / norm2_grad) * ((1.0 / std::sqrt(Lalpha)) - 1.0);
    sqrtG += scale * (gradtheta * gradtheta.t());
  }
  
}


double RMLMC::Etheta(const vec& theta, const vec& y, const vec& h, 
                     const vec& l, const int T, const vec& v, const double alpha)
{
  double lposttheta = lpost_theta(theta, y, h, l, T);
  vec gradtheta = grad_theta(theta, y, h, l, T);
  
  double norm2_grad = dot(gradtheta, gradtheta);
  double dotgradv = dot(gradtheta, v);
  
  double energy = -lposttheta - 0.5 * std::log(1.0 + alpha * alpha * norm2_grad);
  energy += 0.5 * dot(v, v) + 0.5 * alpha * alpha * dotgradv * dotgradv;
  
  return energy;
}

vec RMLMC::updatetheta(const vec& theta_cur, const vec& y, const vec& h_cur, 
                       const vec& l_cur, const int T, const double alpha,
                       const int nlpfrog, const double epsilon, 
                       int& acc, int& div)
{
  // 1. Validação de entradas
  if (y.n_elem != (uword)(T + 1))
    Rcpp::stop("y deve ter T+1 elementos: (y0, y1, ..., yT).");
  if (h_cur.n_elem < (uword)T || l_cur.n_elem < (uword)T)
    Rcpp::stop("h e lambda devem ter pelo menos T elementos.");
  
  const SvmStats st = make_stats(y, h_cur, l_cur, T);
  const double alph2 = alpha * alpha;
  
  // 2. Amostragem da velocidade inicial v ~ N(0, G^-1)
  mat sqrtG(8, 8, fill::zeros);
  sqrtG_theta(theta_cur, y, h_cur, l_cur, T, alpha, sqrtG);
  vec vcur = sqrtG * randn<vec>(8);
  
  // 3. Energia Inicial
  double Einit = Etheta(theta_cur, y, h_cur, l_cur, T, vcur, alpha);
  if (!std::isfinite(Einit)) {
    div += 1;
    return theta_cur;
  }
  
  double deltadet = 0.0;
  vec thetaprop = theta_cur;
  vec vprop = vcur;
  
  vec gradtheta(8);
  mat H(8, 8, fill::zeros);
  
  // 4. Integrador Explícito Leapfrog
  for (int jgl = 1; jgl <= nlpfrog; ++jgl) 
  {
    // =============================================================
    // MEIO-PASSO 1: Atualização da Velocidade na posição atual
    // =============================================================
    gradtheta = grad_theta(thetaprop, y, h_cur, l_cur, T);
    H_theta(thetaprop, st, H, T);
    
    double norm2_grad = dot(gradtheta, gradtheta);
    double Lalpha = 1.0 + alph2 * norm2_grad;
    vec Hv = H * vprop;
    double gHv = dot(gradtheta, Hv);
    
    // Proteção e Ajuste de Determinante (Entrada)
    double det_in_1 = 1.0 + (0.5 * alph2 * epsilon / Lalpha) * gHv;
    if (!std::isfinite(det_in_1) || std::abs(det_in_1) < 1e-12) {
      div += 1;
      return theta_cur;
    }
    deltadet -= std::log(std::abs(det_in_1));
    
    // Vetores auxiliares e atualização v_half
    vec w1 = gradtheta + 0.5 * epsilon * Hv;
    vec Hg = H * gradtheta;
    double gv = dot(gradtheta, vprop);
    vec F2_1 = (alph2 * gv + 0.5 * epsilon) * gradtheta - (0.5 * alph2 * epsilon / Lalpha) * Hg + vprop;
    
    double d1 = dot(w1, gradtheta) + (1.0 / alph2);
    if (!std::isfinite(d1) || std::abs(d1) < 1e-12) {
      div += 1;
      return theta_cur;
    }
    vprop = F2_1 - (dot(w1, F2_1) / d1) * gradtheta;
    
    if (!vprop.is_finite()) {
      div += 1;
      return theta_cur;
    }
    
    // Proteção e Ajuste de Determinante (Saída)
    gHv = dot(gradtheta, H * vprop);
    double det_out_1 = 1.0 - (0.5 * alph2 * epsilon / Lalpha) * gHv;
    if (!std::isfinite(det_out_1) || std::abs(det_out_1) < 1e-12) {
      div += 1;
      return theta_cur;
    }
    deltadet += std::log(std::abs(det_out_1));
    
    // =============================================================
    // PASSO INTEIRO: Atualização da Posição
    // =============================================================
    thetaprop += epsilon * vprop;
    
    if (!thetaprop.is_finite()) {
      div += 1;
      return theta_cur;
    }
    
    // =============================================================
    // MEIO-PASSO 2: Atualização da Velocidade na nova posição
    // =============================================================
    gradtheta = grad_theta(thetaprop, y, h_cur, l_cur, T);
    H_theta(thetaprop, st, H, T);
    
    norm2_grad = dot(gradtheta, gradtheta);
    Lalpha = 1.0 + alph2 * norm2_grad;
    Hv = H * vprop;
    gHv = dot(gradtheta, Hv);
    
    // Proteção e Ajuste de Determinante (Entrada)
    double det_in_2 = 1.0 + (0.5 * alph2 * epsilon / Lalpha) * gHv;
    if (!std::isfinite(det_in_2) || std::abs(det_in_2) < 1e-12) {
      div += 1;
      return theta_cur;
    }
    deltadet -= std::log(std::abs(det_in_2));
    
    vec w2 = gradtheta + 0.5 * epsilon * Hv;
    Hg = H * gradtheta;
    gv = dot(gradtheta, vprop);
    vec F2_2 = (alph2 * gv + 0.5 * epsilon) * gradtheta - (0.5 * alph2 * epsilon / Lalpha) * Hg + vprop;
    
    double d2 = dot(w2, gradtheta) + (1.0 / alph2);
    if (!std::isfinite(d2) || std::abs(d2) < 1e-12) {
      div += 1;
      return theta_cur;
    }
    vprop = F2_2 - (dot(w2, F2_2) / d2) * gradtheta;
    
    if (!vprop.is_finite()) {
      div += 1;
      return theta_cur;
    }
    
    // Proteção e Ajuste de Determinante (Saída)
    gHv = dot(gradtheta, H * vprop);
    double det_out_2 = 1.0 - (0.5 * alph2 * epsilon / Lalpha) * gHv;
    if (!std::isfinite(det_out_2) || std::abs(det_out_2) < 1e-12) {
      div += 1;
      return theta_cur;
    }
    deltadet += std::log(std::abs(det_out_2));
  }
  
  // 5. Avaliação da Energia Final e Aceitação MH
  double Efinal = Etheta(thetaprop, y, h_cur, l_cur, T, vprop, alpha);
  
  double log_alpha = Einit - Efinal + deltadet;
  if (!std::isfinite(log_alpha)) {
    div += 1;
    return theta_cur;
  }
  
  double aMH = std::min(0.0, log_alpha);
  double uMH = std::log(R::unif_rand());
  
  if (uMH < aMH) {
    acc += 1;
    return thetaprop;
  } else {
    return theta_cur;
  }
}

// #############################################################################
// ############################ h ##############################################
// #############################################################################


double RMLMC::lpost_h(const vec& h, const vec& theta, const vec& y, const vec& l, 
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

vec RMLMC::grad_h(const vec& h, const vec& theta, const vec& y, const vec& l, 
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

void RMLMC::H_h(const vec& h_star, const vec& l_cur, const vec& theta, 
                const vec& y, const int T, vec& diagH, vec& subH)
{
  diagH.zeros(T);
  subH.zeros(T - 1);
  
  // 1. Unpack correto de theta (alinhado com lpost_h e grad_h)
  const double b0    = theta[0];
  const double b1    = std::tanh(theta[1]); 
  const double b2    = theta[2];
  const double mu    = theta[3];
  const double phi   = std::tanh(theta[4]);
  const double delta = theta[5];
  const double tau   = std::exp(theta[6]);
  
  const double phi2 = phi * phi;
  const double tau2 = tau * tau;
  const double inv_tau2 = 1.0 / tau2;
  
  // 2. Alinhamento de vetores
  const vec y_obs  = y.subvec(1, T);       // y_1, ..., y_T
  const vec y_lag1 = y.subvec(0, T - 1);   // y_0, ..., y_{T-1}
  const vec l_sub  = l_cur.subvec(0, T - 1);
  
  vec exp_h = exp(h_star);
  vec mu_t  = b0 + b1 * y_lag1 + b2 * exp_h;
  
  vec s_t   = sqrt(l_sub) % exp(-0.5 * h_star);
  vec e_t   = y_obs - mu_t;
  vec zt    = e_t % s_t;
  vec z2_t  = square(zt);
  
  // 3. Vetores auxiliares A_t e B_t para t = 1..T-1 (indices 0..T-2)
  vec At = h_star.subvec(1, T - 1) - mu - phi * (h_star.subvec(0, T - 2) - mu) - delta * zt.subvec(0, T - 2);
  vec Bt = delta * b2 * exp(0.5 * h_star.subvec(0, T - 2)) % sqrt(l_sub.subvec(0, T - 2)) 
    + 0.5 * delta * zt.subvec(0, T - 2) - phi;
  
  // 4. Preenchimento da Subdiagonal (indices 0..T-2)
  for (int t = 0; t < T - 1; ++t) {
    subH(t) = -inv_tau2 * Bt(t);
  }
  
  // 5. Preenchimento da Diagonal Principal
  for (int t = 0; t < T; ++t) {
    double v_t = 0.0;
    
    if (t == 0) {
      v_t = inv_tau2 * (Bt(0) * Bt(0) - 0.25 * delta * zt(0) * At(0) + (1.0 - phi2));
    } else if (t < T - 1) {
      v_t = inv_tau2 * (Bt(t) * Bt(t) - 0.25 * delta * zt(t) * At(t) + 1.0);
    } else { // t == T - 1
      v_t = inv_tau2;
    }
    
    diagH(t) = - b2 * b2 * l_sub(t) * exp_h(t) 
      - b2 * l_sub(t) * (y_obs(t) - mu_t(t)) 
      - 0.5 * z2_t(t) 
      - v_t;
  }
}

vec RMLMC::mult_H_tridiag(const vec& diagH, const vec& subH, const vec& x, const int T)
{
  vec Hx(T);
  
  // Primeiro elemento (t = 0)
  Hx(0) = diagH(0) * x(0) + subH(0) * x(1);
  
  // Elementos intermediários (t = 1 até T-2)
  for (int t = 1; t < T - 1; ++t) {
    Hx(t) = subH(t - 1) * x(t - 1) + diagH(t) * x(t) + subH(t) * x(t + 1);
  }
  
  // Último elemento (t = T - 1)
  Hx(T - 1) = subH(T - 2) * x(T - 2) + diagH(T - 1) * x(T - 1);
  
  return Hx;
}

vec RMLMC::sample_v_h(const vec& h_cur, const vec& l_cur, const vec& theta_cur, 
                      const vec& y, const int T, const double alpha)
{
  // 1. Vetor de gradiente em relação a h_T
  vec gradh = grad_h(h_cur, theta_cur, y,  l_cur, T);
  
  // 2. Norma ao quadrado do gradiente e fator de escala de Monge (L_alpha)
  double norm2_grad = dot(gradh, gradh);
  double Lalpha = 1.0 + alpha * alpha * norm2_grad;
  
  // 3. Amostragem de ruído gaussiano padrão z ~ N(0, I_T)
  vec z = randn<vec>(T);
  double gz = dot(gradh, z);
  
  // 4. Limiar para evitar divisão por zero próximo do modo (||grad|| -> 0)
  const double tol = 1e-8;
  double scale;
  
  if (norm2_grad < tol) {
    // Aproximação limite quando ||grad||^2 -> 0: -alpha^2 / 2
    scale = -0.5 * alpha * alpha;
  } else {
    // Fórmula exata: (1 / ||grad||^2) * (1 / sqrt(L_alpha) - 1)
    scale = (1.0 / norm2_grad) * ((1.0 / std::sqrt(Lalpha)) - 1.0);
  }
  
  // v = z + scale * (grad_h^T * z) * grad_h  --> Complexidade O(T)
  return z + (scale * gz) * gradh;
}

double RMLMC::Eh(const vec& h_cur, const vec& l_cur, const vec& theta, const vec& y, 
                 const int T, const vec& v, const double alpha)
{
  // 1. Log-posterior em relação a h_T
  double lposth = lpost_h(h_cur, theta, y, l_cur, T);
  
  // 2. Gradiente do log-alvo em relação a h_T
  vec gradh = grad_h(h_cur, theta, y,  l_cur, T);
  
  // 3. Produtos escalares em O(T)
  double norm2_grad = dot(gradh, gradh);
  double dotgradv = dot(gradh, v);
  
  // 4. Energia = -log_posterior - 0.5*log(det G_M) + Kinetic_Energy
  double energy = -lposth - 0.5 * std::log(1.0 + alpha * alpha * norm2_grad);
  energy += 0.5 * dot(v, v) + 0.5 * alpha * alpha * dotgradv * dotgradv;
  
  return energy;
}


vec RMLMC::updateh(const vec& h_cur, const vec& l_cur, const vec& theta_cur,
                   const vec& y, const int T, const double alpha,
                   const int nlpfrog, const double epsilon,
                   int& acc, int& div)
{
  if (y.n_elem != (uword)(T + 1))
    Rcpp::stop("y deve ter T+1 elementos.");
  if (h_cur.n_elem < (uword)T || l_cur.n_elem < (uword)T)
    Rcpp::stop("h e lambda devem ter pelo menos T elementos.");
  
  const double alph2 = alpha * alpha;
  
  // 1. Amostragem da velocidade inicial v ~ N(0, G^-1) em O(T)
  vec vcur = sample_v_h(h_cur, l_cur, theta_cur, y, T, alpha);
  
  // 2. Energia Inicial
  double Einit = Eh(h_cur, l_cur, theta_cur, y, T, vcur, alpha);
  if (!std::isfinite(Einit)) {
    div += 1;
    return h_cur;
  }
  
  double deltadet = 0.0;
  vec hprop = h_cur;
  vec vprop = vcur;
  
  vec gradh(T), diagH(T), subH(T - 1);
  
  for (int jgl = 1; jgl <= nlpfrog; ++jgl)
  {
    // =============================================================
    // MEIO-PASSO 1: Atualização da Velocidade em hprop
    // =============================================================
    gradh = grad_h(hprop, theta_cur, y,  l_cur, T);
    H_h(hprop, l_cur, theta_cur, y, T, diagH, subH);
    
    double norm2_grad = dot(gradh, gradh);
    double Lalpha = 1.0 + alph2 * norm2_grad;
    
    // H * v via operador tridiagonal em O(T)
    vec Hv = mult_H_tridiag(diagH, subH, vprop, T);
    double gHv = dot(gradh, Hv);
    
    // Ajuste de Determinante (Entrada)
    double det_in_1 = 1.0 + (0.5 * alph2 * epsilon / Lalpha) * gHv;
    if (!std::isfinite(det_in_1) || std::abs(det_in_1) < 1e-12) {
      div += 1;
      return h_cur;
    }
    deltadet -= std::log(std::abs(det_in_1));
    
    vec w1 = gradh + 0.5 * epsilon * Hv;
    vec Hg = mult_H_tridiag(diagH, subH, gradh, T); // O(T)
    double gh_v = dot(gradh, vprop);
    
    vec F2_1 = (alph2 * gh_v + 0.5 * epsilon) * gradh 
    - (0.5 * alph2 * epsilon / Lalpha) * Hg + vprop;
    
    double d1 = dot(w1, gradh) + (1.0 / alph2);
    if (!std::isfinite(d1) || std::abs(d1) < 1e-12) {
      div += 1;
      return h_cur;
    }
    vprop = F2_1 - (dot(w1, F2_1) / d1) * gradh;
    
    if (!vprop.is_finite()) {
      div += 1;
      return h_cur;
    }
    
    // Ajuste de Determinante (Saída)
    Hv = mult_H_tridiag(diagH, subH, vprop, T);
    gHv = dot(gradh, Hv);
    double det_out_1 = 1.0 - (0.5 * alph2 * epsilon / Lalpha) * gHv;
    if (!std::isfinite(det_out_1) || std::abs(det_out_1) < 1e-12) {
      div += 1;
      return h_cur;
    }
    deltadet += std::log(std::abs(det_out_1));
    
    // =============================================================
    // PASSO INTEIRO: Atualização da Posição h
    // =============================================================
    hprop += epsilon * vprop;
    
    if (!hprop.is_finite()) {
      div += 1;
      return h_cur;
    }
    
    // =============================================================
    // MEIO-PASSO 2: Atualização da Velocidade na nova posição
    // =============================================================
    gradh = grad_h(hprop, theta_cur, y,  l_cur, T);
    H_h(hprop, l_cur, theta_cur, y, T, diagH, subH);
    
    norm2_grad = dot(gradh, gradh);
    Lalpha = 1.0 + alph2 * norm2_grad;
    Hv = mult_H_tridiag(diagH, subH, vprop, T);
    gHv = dot(gradh, Hv);
    
    double det_in_2 = 1.0 + (0.5 * alph2 * epsilon / Lalpha) * gHv;
    if (!std::isfinite(det_in_2) || std::abs(det_in_2) < 1e-12) {
      div += 1;
      return h_cur;
    }
    deltadet -= std::log(std::abs(det_in_2));
    
    vec w2 = gradh + 0.5 * epsilon * Hv;
    Hg = mult_H_tridiag(diagH, subH, gradh, T);
    gh_v = dot(gradh, vprop);
    
    vec F2_2 = (alph2 * gh_v + 0.5 * epsilon) * gradh 
    - (0.5 * alph2 * epsilon / Lalpha) * Hg + vprop;
    
    double d2 = dot(w2, gradh) + (1.0 / alph2);
    if (!std::isfinite(d2) || std::abs(d2) < 1e-12) {
      div += 1;
      return h_cur;
    }
    vprop = F2_2 - (dot(w2, F2_2) / d2) * gradh;
    
    if (!vprop.is_finite()) {
      div += 1;
      return h_cur;
    }
    
    Hv = mult_H_tridiag(diagH, subH, vprop, T);
    gHv = dot(gradh, Hv);
    double det_out_2 = 1.0 - (0.5 * alph2 * epsilon / Lalpha) * gHv;
    if (!std::isfinite(det_out_2) || std::abs(det_out_2) < 1e-12) {
      div += 1;
      return h_cur;
    }
    deltadet += std::log(std::abs(det_out_2));
  }
  
  // 3. Avaliação da Energia Final e Aceitação MH
  double Efinal = Eh(hprop, l_cur, theta_cur, y, T, vprop, alpha);
  
  double log_alpha = Einit - Efinal + deltadet;
  if (!std::isfinite(log_alpha)) {
    div += 1;
    return h_cur;
  }
  
  double aMH = std::min(0.0, log_alpha);
  double uMH = std::log(randu());
  
  if (uMH < aMH) {
    acc += 1;
    return hprop;
  } else {
    return h_cur;
  }
}



// #############################################################################
// ############################ l ##############################################
// #############################################################################

double RMLMC::lpost_l(const vec& l_star, const vec& h, const vec& theta, 
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

vec RMLMC::grad_l(const vec& l_star, const vec& h, const vec& theta, const vec& y, 
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


vec RMLMC::Hdiag_l(const vec& l_star, const vec& h, const vec& theta,
                   const vec& y, const int T)
{
  const double b0 = theta[0], b1 = std::tanh(theta[1]), b2 = theta[2];
  const double mu = theta[3], phi = std::tanh(theta[4]), delta = theta[5];
  const double tau2 = std::exp(2.0 * theta[6]);
  const double Lg = 1.0 / (1.0 + std::exp(-theta[7]));
  const double nu = 40.0 * Lg + 2.0 * (1.0 - Lg);
  
  const vec l  = exp(l_star.head(T));
  const vec hh = h.head(T);
  const vec z  = (y.subvec(1, T) - b0 - b1 * y.subvec(0, T - 1) - b2 * exp(hh))
    % sqrt(l) % exp(-0.5 * hh);
  const vec zc = z.head(T - 1);
  const vec A  = hh.subvec(1, T - 1) - mu - phi * (hh.subvec(0, T - 2) - mu) - delta * zc;
  
  vec Hd = -0.5 * square(z) - 0.5 * nu * l;
  Hd.head(T - 1) += (delta / (4.0 * tau2)) * zc % (A - delta * zc);
  return Hd;
}

vec RMLMC::sample_v(const vec& l_star, const vec& h, const vec& theta, 
                    const vec& y, const int T, const double alpha)
{
  
  vec gradl = grad_l(l_star, h, theta, y, T);
  double norm2_grad = dot(gradl, gradl);
  double Lalpha = 1.0 + alpha * alpha * norm2_grad;
  
  vec z = randn<vec>(T);
  double gz = dot(gradl, z);
  
  const double tol = 1e-8;
  double scale;
  
  if (norm2_grad < tol) {
    scale = -0.5 * alpha * alpha;
  } else {
    scale = (1.0 / norm2_grad) * ((1.0 / std::sqrt(Lalpha)) - 1.0);
  }
  
  // v = z + scale * (grad_l^T * z) * grad_l  --> Operação O(T)!
  return z + (scale * gz) * gradl;
  
}


double RMLMC::El(const vec& l_star, const vec& h, const vec& theta, const vec& y, 
                 const int T, const vec& v, const double alpha)
{
  double lpostl = lpost_l(l_star,h, theta, y, T);
  vec gradl = grad_l(l_star, h, theta, y, T);
  
  double norm2_grad = dot(gradl, gradl);
  double dotgradv = dot(gradl, v);
  
  double energy = -lpostl - 0.5 * std::log(1.0 + alpha * alpha * norm2_grad);
  energy += 0.5 * dot(v, v) + 0.5 * alpha * alpha * dotgradv * dotgradv;
  
  return energy;
}


vec RMLMC::updatel(const vec& l_cur, const vec& h_cur, const vec& theta_cur,
                   const vec& y, const int T, const double alpha,
                   const int nlpfrog, const double epsilon,
                   int& acc, int& div)
{
  // 1. Validação de entradas
  if (y.n_elem != (uword)(T + 1))
    Rcpp::stop("y deve ter T+1 elementos: (y0, y1, ..., yT).");
  if (h_cur.n_elem < (uword)T || l_cur.n_elem < (uword)T)
    Rcpp::stop("h e lambda devem ter pelo menos T elementos.");
  
  const double alph2 = alpha * alpha;
  
  // 2. Velocidade inicial v ~ N(0, G^-1), O(T)
  vec vcur = sample_v(l_cur, h_cur, theta_cur, y, T, alpha);
  
  // 3. Energia inicial
  const double Einit = El(l_cur, h_cur, theta_cur, y, T, vcur, alpha);
  if (!std::isfinite(Einit)) {
    div += 1;
    return l_cur;
  }
  
  double deltadet = 0.0;
  vec lprop = l_cur;
  vec vprop = vcur;
  
  // Meio-passo de velocidade na posição 'pos' (H diagonal avaliada em 'pos').
  // Atualiza v e acumula o ajuste do determinante em deltadet.
  // Retorna false se houver divergência.
  auto half_step = [&](const vec& pos, vec& v) -> bool
  {
    const vec g  = grad_l(pos, h_cur, theta_cur, y, T);
    const vec Hd = Hdiag_l(pos, h_cur, theta_cur, y, T);
    
    const double Lalpha = 1.0 + alph2 * dot(g, g);
    const double c  = 0.5 * alph2 * epsilon / Lalpha;
    const vec    Hv = Hd % v;
    
    // determinante (entrada)
    const double det_in = 1.0 + c * dot(g, Hv);
    if (!std::isfinite(det_in) || std::abs(det_in) < 1e-12) return false;
    deltadet -= std::log(std::abs(det_in));
    
    // v_novo = F - g (w'F) / (w'g + 1/alpha^2)
    const vec w = g + 0.5 * epsilon * Hv;
    const vec F = (alph2 * dot(g, v) + 0.5 * epsilon) * g - c * (Hd % g) + v;
    const double d = dot(w, g) + 1.0 / alph2;
    if (!std::isfinite(d) || std::abs(d) < 1e-12) return false;
    v = F - (dot(w, F) / d) * g;
    if (!v.is_finite()) return false;
    
    // determinante (saída)
    const double det_out = 1.0 - c * dot(g, Hd % v);
    if (!std::isfinite(det_out) || std::abs(det_out) < 1e-12) return false;
    deltadet += std::log(std::abs(det_out));
    return true;
  };
  
  // 4. Integrador explícito (meio-passo de v, passo de posição, meio-passo de v)
  for (int jgl = 0; jgl < nlpfrog; ++jgl)
  {
    if (!half_step(lprop, vprop)) { div += 1; return l_cur; }
    
    lprop += epsilon * vprop;
    if (!lprop.is_finite()) { div += 1; return l_cur; }
    
    if (!half_step(lprop, vprop)) { div += 1; return l_cur; }
  }
  
  // 5. Energia final e aceitação MH
  const double Efinal = El(lprop, h_cur, theta_cur, y, T, vprop, alpha);
  const double log_alpha = Einit - Efinal + deltadet;
  if (!std::isfinite(log_alpha)) {
    div += 1;
    return l_cur;
  }
  
  const double aMH = std::min(0.0, log_alpha);
  const double uMH = std::log(randu());
  if (uMH < aMH) {
    acc += 1;
    return lprop;
  }
  return l_cur;
}






