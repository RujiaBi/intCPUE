#ifndef TMB_LIB_INIT
#define TMB_LIB_INIT R_init_intCPUE_template
#endif
#ifndef INTCPUE_ENABLE_OMEGA
#define INTCPUE_ENABLE_OMEGA 1
#endif
#ifndef INTCPUE_ENABLE_EPSILON
#define INTCPUE_ENABLE_EPSILON 1
#endif
#ifndef INTCPUE_ENABLE_FLAG_S
#define INTCPUE_ENABLE_FLAG_S 1
#endif
#define EIGEN_DONT_PARALLELIZE
#include <TMB.hpp>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace Eigen;
using namespace tmbutils;

template<class Type>
struct LOSM_t : vector<SparseMatrix<Type> > {
  LOSM_t(SEXP x){
    this->resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asSparseMatrix<Type>(sm);
    }
  }
};

// Bias-corrected log-normal density.
template<class Type>
Type dlnorm_bc(const Type& x, const Type& meanlog, const Type& sdlog, int give_log = 0) {
  Type sd2 = sdlog * sdlog;
  Type adjusted_meanlog = meanlog - sd2 / 2;
  Type logres = dnorm(log(x), adjusted_meanlog, sdlog, true) - log(x);
  if (give_log) return logres;
  return exp(logres);
}

// PC prior for Matern range and marginal standard deviation.
template <class Type>
Type pc_prior_matern(
  Type range,
  Type sigma,
  Type matern_range,
  Type matern_sigma,
  Type range_prob,
  Type sigma_prob,
  int give_log = 0,
  int share_range = 0
){
  Type d = 2.;
  Type dhalf = d / 2.;
  Type lam1 = -log(range_prob) * pow(matern_range, dhalf);
  Type lam2 = -log(sigma_prob) / matern_sigma;
  Type range_ll = log(dhalf) + log(lam1) +
    (-1.0 - dhalf) * log(range) - lam1 * pow(range, -dhalf);
  Type sigma_ll = log(lam2) - lam2 * sigma;
  Type penalty = sigma_ll;
  if(!share_range) penalty += range_ll;
  if(give_log) return penalty;
  return exp(penalty);
}

// Main TMB objective function.
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // Dimensions and model switches
  DATA_INTEGER(n_i);
  DATA_INTEGER(n_t);
  DATA_INTEGER(n_v);
  DATA_INTEGER(n_f);
  DATA_INTEGER(n_g);
  DATA_INTEGER(n_region);
  DATA_INTEGER(use_pop_spatiotemporal_rw);
  DATA_INTEGER(use_pop_spatiotemporal_ar1);
  DATA_INTEGER(use_q_diffs_time);
  DATA_INTEGER(use_flag_sd);
  DATA_INTEGER(has_component1);
  DATA_INTEGER(all_obs_use_component1);
  DATA_INTEGER(use_split_component1_design);
  DATA_INTEGER(use_poisson_link);
  DATA_INTEGER(use_flag_f_random);

  // Observation-level data and component-1 indexing
  DATA_VECTOR(b_i);
  DATA_IVECTOR(e_i);
  DATA_IVECTOR(t_i);
  DATA_IVECTOR(v_i);
  DATA_IVECTOR(f_i);
  DATA_IVECTOR(obs_uses_component1);
  DATA_IVECTOR(v_component1_keep);
  DATA_IVECTOR(v_component1_col);
  DATA_IVECTOR(flag_component1_keep);
  DATA_IVECTOR(flag_component1_col);
  DATA_IMATRIX(flag_t_index_1);
  DATA_IMATRIX(flag_t_index_2);
  DATA_VECTOR(area_g);
  DATA_IVECTOR(region_g);

  // Catchability smooth and factor design objects
  DATA_INTEGER(has_smooths_catch);
  DATA_MATRIX(Xs_catch);
  DATA_STRUCT(Zs_catch, LOSM_t);
  DATA_IVECTOR(b_smooth_start_catch);
  DATA_INTEGER(has_fixed_factors_catch);
  DATA_INTEGER(has_random_factors_catch);
  DATA_MATRIX(Xf_catch);
  DATA_STRUCT(Zf_catch, LOSM_t);
  DATA_IVECTOR(b_factor_start_catch);
  DATA_INTEGER(has_smooths_catch_1);
  DATA_MATRIX(Xs_catch_1);
  DATA_STRUCT(Zs_catch_1, LOSM_t);
  DATA_IVECTOR(b_smooth_start_catch_1);
  DATA_INTEGER(has_fixed_factors_catch_1);
  DATA_INTEGER(has_random_factors_catch_1);
  DATA_MATRIX(Xf_catch_1);
  DATA_STRUCT(Zf_catch_1, LOSM_t);
  DATA_IVECTOR(b_factor_start_catch_1);

  // Population smooth and factor design objects
  DATA_INTEGER(has_smooths_pop);
  DATA_MATRIX(Xs_pop_i);
  DATA_STRUCT(Zs_pop_i, LOSM_t);
  DATA_MATRIX(Xs_pop_g);
  DATA_STRUCT(Zs_pop_g, LOSM_t);
  DATA_IVECTOR(b_smooth_start_pop);
  DATA_INTEGER(has_fixed_factors_pop);
  DATA_INTEGER(has_random_factors_pop);
  DATA_INTEGER(has_pop_intercept);
  DATA_MATRIX(Xf_pop_i);
  DATA_STRUCT(Zf_pop_i, LOSM_t);
  DATA_MATRIX(Xf_pop_g);
  DATA_STRUCT(Zf_pop_g, LOSM_t);
  DATA_IVECTOR(b_factor_start_pop);
  DATA_INTEGER(has_smooths_pop_1);
  DATA_MATRIX(Xs_pop_i_1);
  DATA_STRUCT(Zs_pop_i_1, LOSM_t);
  DATA_MATRIX(Xs_pop_g_1);
  DATA_STRUCT(Zs_pop_g_1, LOSM_t);
  DATA_IVECTOR(b_smooth_start_pop_1);
  DATA_INTEGER(has_fixed_factors_pop_1);
  DATA_INTEGER(has_random_factors_pop_1);
  DATA_MATRIX(Xf_pop_i_1);
  DATA_STRUCT(Zf_pop_i_1, LOSM_t);
  DATA_MATRIX(Xf_pop_g_1);
  DATA_STRUCT(Zf_pop_g_1, LOSM_t);
  DATA_IVECTOR(b_factor_start_pop_1);

  // Parameters
  // Observation model
  PARAMETER(ln_sd);
  PARAMETER_VECTOR(ln_sd_flag);

  // SPDE hyperparameters
#if INTCPUE_ENABLE_OMEGA || INTCPUE_ENABLE_EPSILON || INTCPUE_ENABLE_FLAG_S
  PARAMETER_VECTOR(ln_H_input);
  PARAMETER_VECTOR(ln_range_1);
#endif
#if INTCPUE_ENABLE_OMEGA
  PARAMETER_VECTOR(ln_sigma_0_1);
#endif
#if INTCPUE_ENABLE_EPSILON
  PARAMETER_VECTOR(ln_sigma_t_1);
  PARAMETER_VECTOR(transf_rho_1);
#endif
#if INTCPUE_ENABLE_OMEGA || INTCPUE_ENABLE_EPSILON || INTCPUE_ENABLE_FLAG_S
  PARAMETER_VECTOR(ln_range_2);
#endif
#if INTCPUE_ENABLE_OMEGA
  PARAMETER_VECTOR(ln_sigma_0_2);
#endif
#if INTCPUE_ENABLE_EPSILON
  PARAMETER_VECTOR(ln_sigma_t_2);
  PARAMETER_VECTOR(transf_rho_2);
#endif

  // Vessel random effects
  PARAMETER_VECTOR(ves_v_1);
  PARAMETER_VECTOR(ves_v_2);
  PARAMETER_VECTOR(ves_ln_std_dev_1);
  PARAMETER(ves_ln_std_dev_2);

  // Temporal effects
  PARAMETER_VECTOR(yq_t_1);
  PARAMETER_VECTOR(yq_t_2);

  // Spatial and spatiotemporal fields
#if INTCPUE_ENABLE_OMEGA
  PARAMETER_VECTOR(omega_s_1);
  PARAMETER_VECTOR(omega_s_2);
#endif
#if INTCPUE_ENABLE_EPSILON
  PARAMETER_MATRIX(epsilon_st_1);
  PARAMETER_MATRIX(epsilon_st_2);
#endif

  // Flag systematic differences
  PARAMETER_VECTOR(flag_f_1);
  PARAMETER_VECTOR(flag_f_2);
  PARAMETER_VECTOR(flag_ln_std_dev_1);
  PARAMETER(flag_ln_std_dev_2);

  // Flag temporal differences
  PARAMETER_VECTOR(flag_t_1);
  PARAMETER_VECTOR(flag_t_2);
  PARAMETER_VECTOR(flag_t_ln_std_dev_1);
  PARAMETER_VECTOR(flag_t_ln_std_dev_2);

  // Flag spatial differences
#if INTCPUE_ENABLE_FLAG_S
  PARAMETER_MATRIX(flag_s_1);
  PARAMETER_MATRIX(flag_s_2);
  PARAMETER_VECTOR(ln_sigma_flag_1);
  PARAMETER_VECTOR(ln_sigma_flag_2);
#endif

  // Catchability smoothers and factors
  PARAMETER_MATRIX(bs_catch);
  PARAMETER_MATRIX(b_smooth_catch);
  PARAMETER_MATRIX(ln_smooth_sigma_catch);
  PARAMETER_VECTOR(bs_catch_1);
  PARAMETER_VECTOR(b_smooth_catch_1);
  PARAMETER_VECTOR(ln_smooth_sigma_catch_1);
  PARAMETER_MATRIX(bf_catch);
  PARAMETER_MATRIX(b_factor_catch);
  PARAMETER_MATRIX(ln_factor_sigma_catch);
  PARAMETER_VECTOR(bf_catch_1);
  PARAMETER_VECTOR(b_factor_catch_1);
  PARAMETER_VECTOR(ln_factor_sigma_catch_1);

  // Population smoothers and factors
  PARAMETER_MATRIX(bs_pop);
  PARAMETER_MATRIX(b_smooth_pop);
  PARAMETER_MATRIX(ln_smooth_sigma_pop);
  PARAMETER_VECTOR(bs_pop_1);
  PARAMETER_VECTOR(b_smooth_pop_1);
  PARAMETER_VECTOR(ln_smooth_sigma_pop_1);
  PARAMETER_VECTOR(pop_intercept);
  PARAMETER_MATRIX(bf_pop);
  PARAMETER_MATRIX(b_factor_pop);
  PARAMETER_MATRIX(ln_factor_sigma_pop);
  PARAMETER_VECTOR(bf_pop_1);
  PARAMETER_VECTOR(b_factor_pop_1);
  PARAMETER_VECTOR(ln_factor_sigma_pop_1);

  // Global variables
  Type nll = 0;
  Type nll_prior = 0;
  Type nll_penalty = 0;
  Type sd = exp(ln_sd);
  bool use_component1 = (has_component1 == 1);
  int shared_comp2_col = (use_component1 && use_split_component1_design == 0) ? 1 : 0;
  int pop_comp2_col = use_component1 ? 1 : 0;
  int pop_spatial_enabled = INTCPUE_ENABLE_OMEGA;
  int pop_spatiotemporal_enabled = INTCPUE_ENABLE_EPSILON;
  int q_diffs_spatial_enabled = INTCPUE_ENABLE_FLAG_S;

  vector<Type> sd_flag(n_f);
  for (int f=0; f<n_f; f++) sd_flag(f) = exp(ln_sd_flag(f));

  // Random-effect standard deviations
  Type flag_std_dev_1 = Type(0.0);
  Type flag_t_std_dev_1 = Type(0.0);
  Type ves_std_dev_1 = Type(0.0);
  if (use_component1) {
    if (flag_ln_std_dev_1.size() > 0) flag_std_dev_1 = exp(flag_ln_std_dev_1(0));
    if (flag_t_ln_std_dev_1.size() > 0) flag_t_std_dev_1 = exp(flag_t_ln_std_dev_1(0));
    if (ves_ln_std_dev_1.size() > 0) ves_std_dev_1 = exp(ves_ln_std_dev_1(0));
  }

  Type flag_t_std_dev_2 = Type(0.0);
  if (flag_t_ln_std_dev_2.size() > 0) flag_t_std_dev_2 = exp(flag_t_ln_std_dev_2(0));
  Type ves_std_dev_2 = exp(ves_ln_std_dev_2);
  Type flag_std_dev_2 = exp(flag_ln_std_dev_2);

#if INTCPUE_ENABLE_OMEGA || INTCPUE_ENABLE_EPSILON || INTCPUE_ENABLE_FLAG_S
  // SPDE precision matrices
  bool need_Q_1 = use_component1 &&
    (INTCPUE_ENABLE_OMEGA == 1 || INTCPUE_ENABLE_EPSILON == 1);
  bool need_Q_2 = (INTCPUE_ENABLE_OMEGA == 1 || INTCPUE_ENABLE_EPSILON == 1);
#if INTCPUE_ENABLE_FLAG_S
  need_Q_1 = need_Q_1 || (use_component1 && ln_sigma_flag_1.size() > 0);
  need_Q_2 = need_Q_2 || (ln_sigma_flag_2.size() > 0);
#endif

  int share_range_1 = 0;
  int share_range_2 = 0;

  SparseMatrix<Type> Q_1;
  SparseMatrix<Type> Q_2;
  if (need_Q_1 || need_Q_2) {
    DATA_STRUCT(spde, spde_aniso_t);
    // Anisotropy matrix, parameterized so det(H)=1.
    matrix<Type> H(2,2);
    H(0,0) = exp(ln_H_input(0));
    H(1,0) = ln_H_input(1);
    H(0,1) = ln_H_input(1);
    H(1,1) = (Type(1.0) + ln_H_input(1) * ln_H_input(1)) / exp(ln_H_input(0));
    if (need_Q_1) {
      Type range_1 = exp(ln_range_1(0));
      Type kappa_1 = sqrt(Type(8.0)) / range_1;
      Q_1 = Q_spde(spde, kappa_1, H);
    }
    if (need_Q_2) {
      Type range_2 = exp(ln_range_2(0));
      Type kappa_2 = sqrt(Type(8.0)) / range_2;
      Q_2 = Q_spde(spde, kappa_2, H);
    }
  }
#endif

#if INTCPUE_ENABLE_OMEGA
  // Encounter component: time-constant spatial field (omega)
  if (use_component1) {
    DATA_SCALAR(matern_range);
    DATA_SCALAR(range_prob);
    DATA_SCALAR(matern_sigma_0);
    DATA_SCALAR(sigma_prob);
    Type range_1 = exp(ln_range_1(0));
    Type sigma_0_1 = exp(ln_sigma_0_1(0));
    Type kappa_1 = sqrt(Type(8.0)) / range_1;
    nll_prior -= pc_prior_matern(range_1, sigma_0_1, matern_range, matern_sigma_0,
      range_prob, sigma_prob, 1, share_range_1);
    share_range_1 = 1;
    Type tau_0_1 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_1 * sigma_0_1);
    nll += SCALE(GMRF(Q_1), Type(1.0) / tau_0_1)(omega_s_1);
  }
#endif

#if INTCPUE_ENABLE_EPSILON
  // Encounter component: spatiotemporal field (epsilon)
  if (use_component1) {
    DATA_SCALAR(matern_range);
    DATA_SCALAR(range_prob);
    DATA_SCALAR(matern_sigma_t);
    DATA_SCALAR(sigma_prob);
    Type range_1 = exp(ln_range_1(0));
    Type sigma_t_1 = exp(ln_sigma_t_1(0));
    Type rho_1 = Type(2.0) / (Type(1.0) + exp(-Type(2.0) * transf_rho_1(0))) - Type(1.0);
    Type kappa_1 = sqrt(Type(8.0)) / range_1;
    nll_prior -= pc_prior_matern(range_1, sigma_t_1, matern_range, matern_sigma_t,
      range_prob, sigma_prob, 1, share_range_1);
    share_range_1 = 1;
    Type tau_t_1 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_1 * sigma_t_1);
    if (use_pop_spatiotemporal_rw == 1) {
      nll += SCALE(GMRF(Q_1), Type(1.0) / tau_t_1)(epsilon_st_1.col(0));
      for (int t=1; t<n_t; t++) {
        nll += SCALE(GMRF(Q_1), Type(1.0) / tau_t_1)(epsilon_st_1.col(t) - epsilon_st_1.col(t - 1));
      }
    } else if (use_pop_spatiotemporal_ar1 == 1) {
      nll += SCALE(GMRF(Q_1), Type(1.0) / (tau_t_1 * sqrt(Type(1.0) - rho_1 * rho_1)))(epsilon_st_1.col(0));
      for (int t=1; t<n_t; t++) {
        nll += SCALE(GMRF(Q_1), Type(1.0) / tau_t_1)(epsilon_st_1.col(t) - rho_1 * epsilon_st_1.col(t - 1));
      }
    } else {
      for (int t=0; t<n_t; t++) {
        nll += SCALE(GMRF(Q_1), Type(1.0) / tau_t_1)(epsilon_st_1.col(t));
      }
    }
  }
#endif

  // Vessel random-effect priors
  if (n_v > 1) {
    if (ves_v_1.size() == n_v) {
      for(int i=0; i<n_v; i++) {
        nll -= dnorm(ves_v_1(i), Type(0.0), ves_std_dev_1, true);
        nll -= dnorm(ves_v_2(i), Type(0.0), ves_std_dev_2, true);
      }
    } else {
      for(int i=0; i<ves_v_1.size(); i++) {
        nll -= dnorm(ves_v_1(i), Type(0.0), ves_std_dev_1, true);
      }
      for(int i=0; i<n_v; i++) {
        nll -= dnorm(ves_v_2(i), Type(0.0), ves_std_dev_2, true);
      }
    }
  }

  // Flag systematic-difference priors
  if (n_f > 1 && use_flag_f_random == 1) {
    if (flag_f_1.size() == n_f - 1) {
      for(int i=0; i<n_f-1; i++) {
        nll -= dnorm(flag_f_1(i), Type(0.0), flag_std_dev_1, true);
        nll -= dnorm(flag_f_2(i), Type(0.0), flag_std_dev_2, true);
      }
    } else {
      for(int i=0; i<flag_f_1.size(); i++) {
        nll -= dnorm(flag_f_1(i), Type(0.0), flag_std_dev_1, true);
      }
      for(int i=0; i<n_f-1; i++) {
        nll -= dnorm(flag_f_2(i), Type(0.0), flag_std_dev_2, true);
      }
    }
  }

  // Flag temporal-difference priors and reconstruction
  matrix<Type> flag_t_centered_1(use_component1 ? n_t : 0, use_component1 ? n_f - 1 : 0);
  matrix<Type> flag_t_centered_2(n_t, n_f - 1);
  flag_t_centered_1.setZero();
  flag_t_centered_2.setZero();
  if (use_q_diffs_time == 1) {
    for(int j=0; j<n_f-1; j++) {
      for(int t=0; t<n_t; t++) {
        int idx1 = flag_t_index_1(t, j);
        if (use_component1 && idx1 >= 0) {
          flag_t_centered_1(t, j) = flag_t_1(idx1);
          nll -= dnorm(flag_t_1(idx1), Type(0.0), flag_t_std_dev_1, true);
        }
        int idx2 = flag_t_index_2(t, j);
        if (idx2 >= 0) {
          flag_t_centered_2(t, j) = flag_t_2(idx2);
          nll -= dnorm(flag_t_2(idx2), Type(0.0), flag_t_std_dev_2, true);
        }
      }
    }
  }

#if INTCPUE_ENABLE_OMEGA
  // Positive component: time-constant spatial field (omega)
  {
    DATA_SCALAR(matern_range);
    DATA_SCALAR(range_prob);
    DATA_SCALAR(matern_sigma_0);
    DATA_SCALAR(sigma_prob);
    Type range_2 = exp(ln_range_2(0));
    Type sigma_0_2 = exp(ln_sigma_0_2(0));
    Type kappa_2 = sqrt(Type(8.0)) / range_2;
    nll_prior -= pc_prior_matern(range_2, sigma_0_2, matern_range, matern_sigma_0,
      range_prob, sigma_prob, 1, share_range_2);
    share_range_2 = 1;
    Type tau_0_2 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_2 * sigma_0_2);
    nll += SCALE(GMRF(Q_2), Type(1.0) / tau_0_2)(omega_s_2);
  }
#endif

#if INTCPUE_ENABLE_EPSILON
  // Positive component: spatiotemporal field (epsilon)
  {
    DATA_SCALAR(matern_range);
    DATA_SCALAR(range_prob);
    DATA_SCALAR(matern_sigma_t);
    DATA_SCALAR(sigma_prob);
    Type range_2 = exp(ln_range_2(0));
    Type sigma_t_2 = exp(ln_sigma_t_2(0));
    Type rho_2 = Type(2.0) / (Type(1.0) + exp(-Type(2.0) * transf_rho_2(0))) - Type(1.0);
    Type kappa_2 = sqrt(Type(8.0)) / range_2;
    nll_prior -= pc_prior_matern(range_2, sigma_t_2, matern_range, matern_sigma_t,
      range_prob, sigma_prob, 1, share_range_2);
    share_range_2 = 1;
    Type tau_t_2 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_2 * sigma_t_2);
    if (use_pop_spatiotemporal_rw == 1) {
      nll += SCALE(GMRF(Q_2), Type(1.0) / tau_t_2)(epsilon_st_2.col(0));
      for (int t=1; t<n_t; t++) {
        nll += SCALE(GMRF(Q_2), Type(1.0) / tau_t_2)(epsilon_st_2.col(t) - epsilon_st_2.col(t - 1));
      }
    } else if (use_pop_spatiotemporal_ar1 == 1) {
      nll += SCALE(GMRF(Q_2), Type(1.0) / (tau_t_2 * sqrt(Type(1.0) - rho_2 * rho_2)))(epsilon_st_2.col(0));
      for (int t=1; t<n_t; t++) {
        nll += SCALE(GMRF(Q_2), Type(1.0) / tau_t_2)(epsilon_st_2.col(t) - rho_2 * epsilon_st_2.col(t - 1));
      }
    } else {
      for (int t=0; t<n_t; t++) {
        nll += SCALE(GMRF(Q_2), Type(1.0) / tau_t_2)(epsilon_st_2.col(t));
      }
    }
  }
#endif

#if INTCPUE_ENABLE_FLAG_S
  // Flag-specific spatial differences
  {
    DATA_SCALAR(matern_range);
    DATA_SCALAR(range_prob);
    DATA_SCALAR(matern_sigma_flag);
    DATA_SCALAR(sigma_prob);
    Type tau_flag_1 = Type(0.0);
    Type tau_flag_2 = Type(0.0);
    if (use_component1 && flag_s_1.cols() > 0) {
      Type range_1 = exp(ln_range_1(0));
      Type sigma_flag_1 = exp(ln_sigma_flag_1(0));
      Type kappa_1 = sqrt(Type(8.0)) / range_1;
      nll_prior -= pc_prior_matern(range_1, sigma_flag_1, matern_range, matern_sigma_flag,
        range_prob, sigma_prob, 1, share_range_1);
      share_range_1 = 1;
      tau_flag_1 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_1 * sigma_flag_1);
    }
    if (flag_s_2.cols() > 0) {
      Type range_2 = exp(ln_range_2(0));
      Type sigma_flag_2 = exp(ln_sigma_flag_2(0));
      Type kappa_2 = sqrt(Type(8.0)) / range_2;
      nll_prior -= pc_prior_matern(range_2, sigma_flag_2, matern_range, matern_sigma_flag,
        range_prob, sigma_prob, 1, share_range_2);
      share_range_2 = 1;
      tau_flag_2 = Type(1.0) / (Type(2.0) * sqrt(M_PI) * kappa_2 * sigma_flag_2);
    }
    if (use_component1 && flag_s_1.cols() == flag_s_2.cols()) {
      for(int j=0; j<flag_s_2.cols(); j++) {
        nll += SCALE(GMRF(Q_1), Type(1.0) / tau_flag_1)(flag_s_1.col(j));
        nll += SCALE(GMRF(Q_2), Type(1.0) / tau_flag_2)(flag_s_2.col(j));
      }
    } else {
      for(int j=0; j<flag_s_1.cols(); j++) {
        nll += SCALE(GMRF(Q_1), Type(1.0) / tau_flag_1)(flag_s_1.col(j));
      }
      for(int j=0; j<flag_s_2.cols(); j++) {
        nll += SCALE(GMRF(Q_2), Type(1.0) / tau_flag_2)(flag_s_2.col(j));
      }
    }
  }
#endif

  // Linear predictor workspaces
  int n_i_1 = use_component1 ? n_i : 0;
  int n_g_1 = use_component1 ? n_g : 0;
  int n_gt = n_g * n_t;
  int n_gt_1 = use_component1 ? n_gt : 0;
  vector<Type> eta_smooth_catch_1(n_i_1);
  vector<Type> eta_smooth_catch_2(n_i);
  vector<Type> eta_smooth_pop_i_1(n_i_1);
  vector<Type> eta_smooth_pop_i_2(n_i);
  vector<Type> eta_smooth_pop_g_1(n_gt_1);
  vector<Type> eta_smooth_pop_g_2(n_gt);
  vector<Type> eta_factor_catch_1(n_i_1);
  vector<Type> eta_factor_catch_2(n_i);
  vector<Type> eta_factor_pop_i_1(n_i_1);
  vector<Type> eta_factor_pop_i_2(n_i);
  vector<Type> eta_factor_pop_g_1(n_gt_1);
  vector<Type> eta_factor_pop_g_2(n_gt);
  eta_smooth_catch_1.setZero();
  eta_smooth_catch_2.setZero();
  eta_smooth_pop_i_1.setZero();
  eta_smooth_pop_i_2.setZero();
  eta_smooth_pop_g_1.setZero();
  eta_smooth_pop_g_2.setZero();
  eta_factor_catch_1.setZero();
  eta_factor_catch_2.setZero();
  eta_factor_pop_i_1.setZero();
  eta_factor_pop_i_2.setZero();
  eta_factor_pop_g_1.setZero();
  eta_factor_pop_g_2.setZero();

  vector<Type> eta_hat_encounter_i(n_i);
  vector<Type> eta_hat_positive_i(n_i);
  vector<Type> encounter_prob_i(n_i);
  vector<Type> log_positive_mean_i(n_i);
  vector<Type> mu_hat_i(n_i);
  eta_hat_encounter_i.setZero();
  eta_hat_positive_i.setZero();
  encounter_prob_i.setZero();
  log_positive_mean_i.setZero();
  mu_hat_i.setZero();

  // Catchability smooth terms
  if (use_component1 && use_split_component1_design == 1 && has_smooths_catch_1) {
    int n_smooth = b_smooth_start_catch_1.size();
    for (int s=0; s<n_smooth; s++) {
      int k_s = Zs_catch_1(s).cols();
      int start = b_smooth_start_catch_1(s);
      Type smooth_sd0 = exp(ln_smooth_sigma_catch_1(s));
      vector<Type> beta0(k_s);
      for (int j=0; j<k_s; j++) {
        beta0(j) = b_smooth_catch_1(start + j);
        nll -= dnorm(beta0(j), Type(0.0), smooth_sd0, true);
      }
      eta_smooth_catch_1 += Zs_catch_1(s) * beta0;
    }
    if (Xs_catch_1.cols() > 0) {
      eta_smooth_catch_1 += Xs_catch_1 * bs_catch_1;
    }
  }

  if (has_smooths_catch) {
    int n_smooth = b_smooth_start_catch.size();
    for (int s=0; s<n_smooth; s++) {
      int k_s = Zs_catch(s).cols();
      int start = b_smooth_start_catch(s);
      Type smooth_sd0 = Type(1.0);
      if (use_component1 && use_split_component1_design == 0) smooth_sd0 = exp(ln_smooth_sigma_catch(s, 0));
      Type smooth_sd1 = exp(ln_smooth_sigma_catch(s, shared_comp2_col));
      vector<Type> beta0(k_s);
      vector<Type> beta1(k_s);
      for (int j=0; j<k_s; j++) {
        if (use_component1 && use_split_component1_design == 0) {
          beta0(j) = b_smooth_catch(start + j, 0);
          nll -= dnorm(beta0(j), Type(0.0), smooth_sd0, true);
        } else {
          beta0(j) = Type(0.0);
        }
        beta1(j) = b_smooth_catch(start + j, shared_comp2_col);
        nll -= dnorm(beta1(j), Type(0.0), smooth_sd1, true);
      }
      if (use_component1 && use_split_component1_design == 0) eta_smooth_catch_1 += Zs_catch(s) * beta0;
      eta_smooth_catch_2 += Zs_catch(s) * beta1;
    }
    if (Xs_catch.cols() > 0) {
      if (use_component1 && use_split_component1_design == 0) eta_smooth_catch_1 += Xs_catch * vector<Type>(bs_catch.col(0));
      eta_smooth_catch_2 += Xs_catch * vector<Type>(bs_catch.col(shared_comp2_col));
    }
  }

  // Catchability fixed factor terms
  if (use_component1 && use_split_component1_design == 1 && has_fixed_factors_catch_1 && Xf_catch_1.cols() > 0) {
    eta_factor_catch_1 += Xf_catch_1 * bf_catch_1;
  }
  if (has_fixed_factors_catch && Xf_catch.cols() > 0) {
    if (use_component1 && use_split_component1_design == 0) eta_factor_catch_1 += Xf_catch * vector<Type>(bf_catch.col(0));
    eta_factor_catch_2 += Xf_catch * vector<Type>(bf_catch.col(shared_comp2_col));
  }

  // Catchability random factor terms
  if (use_component1 && use_split_component1_design == 1 && has_random_factors_catch_1) {
    int n_factor = b_factor_start_catch_1.size();
    for (int s=0; s<n_factor; s++) {
      int k_s = Zf_catch_1(s).cols();
      int start = b_factor_start_catch_1(s);
      Type factor_sd0 = exp(ln_factor_sigma_catch_1(s));
      vector<Type> beta0(k_s);
      for (int j=0; j<k_s; j++) {
        beta0(j) = b_factor_catch_1(start + j);
        nll -= dnorm(beta0(j), Type(0.0), factor_sd0, true);
      }
      eta_factor_catch_1 += Zf_catch_1(s) * beta0;
    }
  }
  if (has_random_factors_catch) {
    int n_factor = b_factor_start_catch.size();
    for (int s=0; s<n_factor; s++) {
      int k_s = Zf_catch(s).cols();
      int start = b_factor_start_catch(s);
      vector<Type> beta0(k_s);
      if (use_component1 && use_split_component1_design == 0) {
        Type factor_sd0 = exp(ln_factor_sigma_catch(s, 0));
        for (int j=0; j<k_s; j++) {
          beta0(j) = b_factor_catch(start + j, 0);
          nll -= dnorm(beta0(j), Type(0.0), factor_sd0, true);
        }
        eta_factor_catch_1 += Zf_catch(s) * beta0;
      }
      Type factor_sd1 = exp(ln_factor_sigma_catch(s, shared_comp2_col));
      vector<Type> beta1(k_s);
      for (int j=0; j<k_s; j++) {
        beta1(j) = b_factor_catch(start + j, shared_comp2_col);
        nll -= dnorm(beta1(j), Type(0.0), factor_sd1, true);
      }
      eta_factor_catch_2 += Zf_catch(s) * beta1;
    }
  }

  // Population smooth terms on observations and projection grid
  if (use_component1 && use_split_component1_design == 1 && has_smooths_pop_1) {
    int n_smooth = b_smooth_start_pop_1.size();
    for (int s=0; s<n_smooth; s++) {
      int k_s = Zs_pop_i_1(s).cols();
      int start = b_smooth_start_pop_1(s);
      Type smooth_sd0 = exp(ln_smooth_sigma_pop_1(s));
      vector<Type> beta0(k_s);
      for (int j=0; j<k_s; j++) {
        beta0(j) = b_smooth_pop_1(start + j);
        nll -= dnorm(beta0(j), Type(0.0), smooth_sd0, true);
      }
      eta_smooth_pop_i_1 += Zs_pop_i_1(s) * beta0;
      eta_smooth_pop_g_1 += Zs_pop_g_1(s) * beta0;
    }
    if (Xs_pop_i_1.cols() > 0) {
      eta_smooth_pop_i_1 += Xs_pop_i_1 * bs_pop_1;
      eta_smooth_pop_g_1 += Xs_pop_g_1 * bs_pop_1;
    }
  }

  if (has_smooths_pop) {
    int n_smooth = b_smooth_start_pop.size();
    for (int s=0; s<n_smooth; s++) {
      int k_s = Zs_pop_i(s).cols();
      int start = b_smooth_start_pop(s);
      Type smooth_sd0 = Type(1.0);
      if (use_component1 && use_split_component1_design == 0) smooth_sd0 = exp(ln_smooth_sigma_pop(s, 0));
      Type smooth_sd1 = exp(ln_smooth_sigma_pop(s, shared_comp2_col));
      vector<Type> beta0(k_s);
      vector<Type> beta1(k_s);
      for (int j=0; j<k_s; j++) {
        if (use_component1 && use_split_component1_design == 0) {
          beta0(j) = b_smooth_pop(start + j, 0);
          nll -= dnorm(beta0(j), Type(0.0), smooth_sd0, true);
        } else {
          beta0(j) = Type(0.0);
        }
        beta1(j) = b_smooth_pop(start + j, shared_comp2_col);
        nll -= dnorm(beta1(j), Type(0.0), smooth_sd1, true);
      }
      if (use_component1 && use_split_component1_design == 0) {
        eta_smooth_pop_i_1 += Zs_pop_i(s) * beta0;
        eta_smooth_pop_g_1 += Zs_pop_g(s) * beta0;
      }
      eta_smooth_pop_i_2 += Zs_pop_i(s) * beta1;
      eta_smooth_pop_g_2 += Zs_pop_g(s) * beta1;
    }
    if (Xs_pop_i.cols() > 0) {
      if (use_component1 && use_split_component1_design == 0) {
        vector<Type> bs0 = vector<Type>(bs_pop.col(0));
        eta_smooth_pop_i_1 += Xs_pop_i * bs0;
        eta_smooth_pop_g_1 += Xs_pop_g * bs0;
      }
      vector<Type> bs1 = vector<Type>(bs_pop.col(shared_comp2_col));
      eta_smooth_pop_i_2 += Xs_pop_i * bs1;
      eta_smooth_pop_g_2 += Xs_pop_g * bs1;
    }
  }

  // Observation-space spatial projections
  vector<Type> s_effect_1(n_i_1);
  vector<Type> s_effect_2(n_i);
  vector<Type> st_effect_1(n_i_1);
  vector<Type> st_effect_2(n_i);
  vector<Type> flag_s_effect_1(n_i_1);
  vector<Type> flag_s_effect_2(n_i);
  s_effect_1.setZero();
  s_effect_2.setZero();
  st_effect_1.setZero();
  st_effect_2.setZero();
  flag_s_effect_1.setZero();
  flag_s_effect_2.setZero();

#if INTCPUE_ENABLE_OMEGA
  // Project omega from mesh knots to observations.
  DATA_SPARSE_MATRIX(A_is);
  if (use_component1) s_effect_1 = A_is * omega_s_1;
  s_effect_2 = A_is * omega_s_2;
#endif

#if INTCPUE_ENABLE_EPSILON || INTCPUE_ENABLE_FLAG_S
  // Sparse triplet projection for epsilon and/or flag_s.
  DATA_IMATRIX(Ais_ij);
  DATA_VECTOR(Ais_x);

#if INTCPUE_ENABLE_EPSILON && INTCPUE_ENABLE_FLAG_S
  for(int r=0; r<Ais_ij.rows(); r++) {
    int i = Ais_ij(r, 0);
    int s = Ais_ij(r, 1);
    int t_id = t_i(i);
    int f_id = f_i(i);
    Type a_is = Ais_x(r);
    bool use_component1_i = (use_component1 && obs_uses_component1(i) == 1);

    if (use_component1_i) st_effect_1(i) += a_is * epsilon_st_1(s, t_id);
    st_effect_2(i) += a_is * epsilon_st_2(s, t_id);

    if (f_id > 0) {
      if (use_component1_i && flag_component1_keep(f_id) == 1) {
        int f1_col = flag_component1_col(f_id);
        flag_s_effect_1(i) += a_is * flag_s_1(s, f1_col);
      }
      flag_s_effect_2(i) += a_is * flag_s_2(s, f_id - 1);
    }
  }
#elif INTCPUE_ENABLE_EPSILON
  for(int r=0; r<Ais_ij.rows(); r++) {
    int i = Ais_ij(r, 0);
    int s = Ais_ij(r, 1);
    int t_id = t_i(i);
    Type a_is = Ais_x(r);
    if (use_component1 && obs_uses_component1(i) == 1) st_effect_1(i) += a_is * epsilon_st_1(s, t_id);
    st_effect_2(i) += a_is * epsilon_st_2(s, t_id);
  }
#elif INTCPUE_ENABLE_FLAG_S
  for(int r=0; r<Ais_ij.rows(); r++) {
    int i = Ais_ij(r, 0);
    int s = Ais_ij(r, 1);
    int f_id = f_i(i);
    Type a_is = Ais_x(r);
    if (f_id <= 0) continue;
    if (use_component1 && obs_uses_component1(i) == 1 && flag_component1_keep(f_id) == 1) {
      int f1_col = flag_component1_col(f_id);
      flag_s_effect_1(i) += a_is * flag_s_1(s, f1_col);
    }
    flag_s_effect_2(i) += a_is * flag_s_2(s, f_id - 1);
  }
#endif
#endif

  // Population fixed factor terms
  if (use_component1 && use_split_component1_design == 1 && has_fixed_factors_pop_1 && Xf_pop_i_1.cols() > 0) {
    eta_factor_pop_i_1 += Xf_pop_i_1 * bf_pop_1;
    eta_factor_pop_g_1 += Xf_pop_g_1 * bf_pop_1;
  }
  if (has_fixed_factors_pop && Xf_pop_i.cols() > 0) {
    if (use_component1 && use_split_component1_design == 0) {
      vector<Type> bf0 = vector<Type>(bf_pop.col(0));
      eta_factor_pop_i_1 += Xf_pop_i * bf0;
      eta_factor_pop_g_1 += Xf_pop_g * bf0;
    }
    vector<Type> bf1 = vector<Type>(bf_pop.col(shared_comp2_col));
    eta_factor_pop_i_2 += Xf_pop_i * bf1;
    eta_factor_pop_g_2 += Xf_pop_g * bf1;
  }

  // Population random factor terms
  if (use_component1 && use_split_component1_design == 1 && has_random_factors_pop_1) {
    int n_factor = b_factor_start_pop_1.size();
    for (int s=0; s<n_factor; s++) {
      int k_s = Zf_pop_i_1(s).cols();
      int start = b_factor_start_pop_1(s);
      Type factor_sd0 = exp(ln_factor_sigma_pop_1(s));
      vector<Type> beta0(k_s);
      for (int j=0; j<k_s; j++) {
        beta0(j) = b_factor_pop_1(start + j);
        nll -= dnorm(beta0(j), Type(0.0), factor_sd0, true);
      }
      eta_factor_pop_i_1 += Zf_pop_i_1(s) * beta0;
      eta_factor_pop_g_1 += Zf_pop_g_1(s) * beta0;
    }
  }
  if (has_random_factors_pop) {
    int n_factor = b_factor_start_pop.size();
    for (int s=0; s<n_factor; s++) {
      int k_s = Zf_pop_i(s).cols();
      int start = b_factor_start_pop(s);
      vector<Type> beta0(k_s);
      if (use_component1 && use_split_component1_design == 0) {
        Type factor_sd0 = exp(ln_factor_sigma_pop(s, 0));
        for (int j=0; j<k_s; j++) {
          beta0(j) = b_factor_pop(start + j, 0);
          nll -= dnorm(beta0(j), Type(0.0), factor_sd0, true);
        }
        eta_factor_pop_i_1 += Zf_pop_i(s) * beta0;
        eta_factor_pop_g_1 += Zf_pop_g(s) * beta0;
      }
      Type factor_sd1 = exp(ln_factor_sigma_pop(s, shared_comp2_col));
      vector<Type> beta1(k_s);
      for (int j=0; j<k_s; j++) {
        beta1(j) = b_factor_pop(start + j, shared_comp2_col);
        nll -= dnorm(beta1(j), Type(0.0), factor_sd1, true);
      }
      eta_factor_pop_i_2 += Zf_pop_i(s) * beta1;
      eta_factor_pop_g_2 += Zf_pop_g(s) * beta1;
    }
  }

  // Optional population intercept used with f(tid).
  Type pop_intercept_1 = Type(0.0);
  Type pop_intercept_2 = Type(0.0);
  if (has_pop_intercept) {
    if (use_component1) pop_intercept_1 = pop_intercept(0);
    pop_intercept_2 = pop_intercept(pop_comp2_col);
  }

  // Encounter-component linear predictor for delta observations.
  auto eta_component1 = [&](int i, int tid, int vid, int fid) -> Type {
    Type ves_effect_1 = Type(0.0);
    if (n_v > 1 && v_component1_keep(vid) == 1) {
      int v1_col = v_component1_col(vid);
      ves_effect_1 = ves_v_1(v1_col);
    }
    Type yq_effect_1 = yq_t_1(tid);
    Type flag_effect_1 = Type(0.0);
    if (fid > 0 && flag_component1_keep(fid) == 1) {
      int f1_col = flag_component1_col(fid);
      flag_effect_1 = flag_f_1(f1_col);
    }

    Type flag_t_effect_1 = Type(0.0);
    if (use_q_diffs_time == 1 && fid > 0) {
      flag_t_effect_1 = flag_t_centered_1(tid, fid - 1);
    }

    return ves_effect_1 + yq_effect_1 + flag_effect_1 + flag_t_effect_1 +
      flag_s_effect_1(i) + pop_intercept_1 + s_effect_1(i) + st_effect_1(i) +
      eta_smooth_catch_1(i) + eta_smooth_pop_i_1(i) +
      eta_factor_catch_1(i) + eta_factor_pop_i_1(i);
  };

  // Positive-component linear predictor for all observations.
  auto eta_component2 = [&](int i, int tid, int vid, int fid) -> Type {
    Type ves_effect_2 = Type(0.0);
    if (n_v > 1) ves_effect_2 = ves_v_2(vid);
    Type yq_effect_2 = yq_t_2(tid);
    Type flag_effect_2 = Type(0.0);
    if (n_f > 1 && fid > 0) flag_effect_2 = flag_f_2(fid - 1);

    Type flag_t_effect_2 = Type(0.0);
    if (use_q_diffs_time == 1 && fid > 0) {
      int j = fid - 1;
      flag_t_effect_2 = flag_t_centered_2(tid, j);
    }

    return ves_effect_2 + yq_effect_2 + flag_effect_2 + flag_t_effect_2 +
      flag_s_effect_2(i) + pop_intercept_2 + s_effect_2(i) + st_effect_2(i) +
      eta_smooth_catch_2(i) + eta_smooth_pop_i_2(i) +
      eta_factor_catch_2(i) + eta_factor_pop_i_2(i);
  };

  // Positive biomass likelihood.
  auto add_positive_likelihood = [&](int i, int fid, Type logcat) {
    Type sd_i = sd;
    if (use_flag_sd == 1) sd_i = sd_flag(fid);
    nll -= dlnorm_bc(b_i(i), logcat, sd_i, true);
  };

  // Delta observation likelihood; Poisson-link and traditional links differ here.
  auto add_delta_observation = [&](int i, int tid, int vid, int fid) {
    Type eta1 = eta_component1(i, tid, vid, fid);
    Type eta2 = eta_component2(i, tid, vid, fid);

    Type log_one_minus_p = Type(0.0);
    Type logp = Type(0.0);
    Type p = Type(1.0);
    if (use_poisson_link == 1) {
      log_one_minus_p = -exp(eta1);
      logp = logspace_sub(Type(0.0), log_one_minus_p);
      p = exp(logp);
    } else {
      p = Type(1.0) / (Type(1.0) + exp(-eta1));
      logp = log(p);
      log_one_minus_p = log(Type(1.0) - p);
    }

    Type logcat = eta2;
    if (use_poisson_link == 1) logcat = eta1 + eta2 - logp;

    eta_hat_encounter_i(i) = eta1;
    eta_hat_positive_i(i) = eta2;
    encounter_prob_i(i) = p;
    log_positive_mean_i(i) = logcat;
    mu_hat_i(i) = (use_poisson_link == 1) ? exp(eta1 + eta2) : p * exp(eta2);

    nll -= e_i(i) * logp + (1 - e_i(i)) * log_one_minus_p;
    if(e_i(i) == 1) {
      add_positive_likelihood(i, fid, logcat);
    }
  };

  // Positive-only observation likelihood, with encounter probability fixed at 1.
  auto add_positive_observation = [&](int i, int tid, int vid, int fid) {
    Type eta2 = eta_component2(i, tid, vid, fid);
    eta_hat_encounter_i(i) = Type(0.0);
    eta_hat_positive_i(i) = eta2;
    encounter_prob_i(i) = Type(1.0);
    log_positive_mean_i(i) = eta2;
    mu_hat_i(i) = exp(eta2);
    add_positive_likelihood(i, fid, eta2);
  };

  // Route observations by data type: positive-only, delta-only, or mixed.
  if (!use_component1) {
    for(int i=0; i<n_i; i++) {
      add_positive_observation(i, t_i(i), v_i(i), f_i(i));
    }
  } else if (all_obs_use_component1 == 1) {
    for(int i=0; i<n_i; i++) {
      add_delta_observation(i, t_i(i), v_i(i), f_i(i));
    }
  } else {
    for(int i=0; i<n_i; i++) {
      int tid = t_i(i);
      int vid = v_i(i);
      int fid = f_i(i);
      if (obs_uses_component1(i) == 1) {
        add_delta_observation(i, tid, vid, fid);
      } else {
        add_positive_observation(i, tid, vid, fid);
      }
    }
  }

  // Projection-space spatial fields
  vector<Type> s_effect_proj_1(n_g_1);
  vector<Type> s_effect_proj_2(n_g);
  matrix<Type> st_effect_proj_1(n_g_1, use_component1 ? n_t : 0);
  matrix<Type> st_effect_proj_2(n_g, n_t);
  s_effect_proj_1.setZero();
  s_effect_proj_2.setZero();
  st_effect_proj_1.setZero();
  st_effect_proj_2.setZero();
#if INTCPUE_ENABLE_OMEGA || INTCPUE_ENABLE_EPSILON
  DATA_SPARSE_MATRIX(A_gs);
#if INTCPUE_ENABLE_OMEGA
  if (use_component1) s_effect_proj_1 = A_gs * omega_s_1;
  s_effect_proj_2 = A_gs * omega_s_2;
#endif
#if INTCPUE_ENABLE_EPSILON
  if (use_component1) st_effect_proj_1 = A_gs * epsilon_st_1;
  st_effect_proj_2 = A_gs * epsilon_st_2;
#endif
#endif

  // Projection and regional index calculations
  matrix<Type> cpue_density(n_g, n_t);
  vector<Type> mu_total(n_t);
  vector<Type> link_total(n_t);
  int n_index = n_region + 1;
  matrix<Type> mu_index(n_index, n_t);
  matrix<Type> link_index(n_index, n_t);
  mu_total.setZero();
  mu_index.setZero();

  for(int t=0; t<n_t; t++) {
    Type yq_effect_proj_1 = Type(0.0);
    if (use_component1) yq_effect_proj_1 = yq_t_1(t);
    Type yq_effect_proj_2 = yq_t_2(t);
    for(int g=0; g<n_g; g++) {
      int gt = g + n_g * t;
      Type eta1_proj = Type(0.0);
      if (use_component1) {
        eta1_proj = pop_intercept_1 + yq_effect_proj_1 + s_effect_proj_1(g) +
          st_effect_proj_1(g, t) + eta_smooth_pop_g_1(gt) + eta_factor_pop_g_1(gt);
      }
      Type eta2_proj = pop_intercept_2 + yq_effect_proj_2 + s_effect_proj_2(g) +
        st_effect_proj_2(g, t) + eta_smooth_pop_g_2(gt) + eta_factor_pop_g_2(gt);
      Type cpue = exp(eta2_proj);
      if (use_component1) {
        if (use_poisson_link == 1) {
          cpue = exp(eta1_proj + eta2_proj);
        } else {
          Type p_proj = Type(1.0) / (Type(1.0) + exp(-eta1_proj));
          cpue = p_proj * exp(eta2_proj);
        }
      }
      cpue_density(g, t) = cpue;
      Type adds = cpue * area_g(g);
      mu_total(t) += adds;
      int rg = region_g(g);
      if (rg >= 0 && rg < n_region) mu_index(rg, t) += adds;
    }
    mu_index(n_region, t) = mu_total(t);
    link_total(t) = log(mu_total(t));
    for (int r=0; r<n_index; r++) link_index(r, t) = log(mu_index(r, t));
  }

  // Bias-correction target for total or joint region + total index.
  PARAMETER_VECTOR(eps_index);
  if (eps_index.size() > 0) {
    Type S;
    if (eps_index.size() == n_t) {
      for (int t=0; t<n_t; t++) {
        S = mu_total(t);
        S = newton::Tag(S);
        nll_penalty += eps_index(t) * S;
      }
    } else {
      for (int t=0; t<n_t; t++) {
        for (int r=0; r<n_index; r++) {
          int idx = r + n_index * t;
          if (idx < eps_index.size()) {
            S = mu_index(r, t);
            S = newton::Tag(S);
            nll_penalty += eps_index(idx) * S;
          }
        }
      }
    }
  }

  // Final objective
  nll += nll_prior;
  nll += nll_penalty;

  // Reporting
  REPORT(nll_prior);
  REPORT(nll_penalty);
  REPORT(use_pop_spatiotemporal_rw);
  REPORT(use_pop_spatiotemporal_ar1);
  REPORT(pop_spatial_enabled);
  REPORT(pop_spatiotemporal_enabled);
  REPORT(use_q_diffs_time);
  REPORT(q_diffs_spatial_enabled);
  REPORT(use_flag_sd);
  REPORT(sd_flag);
  REPORT(eta_hat_encounter_i);
  REPORT(eta_hat_positive_i);
  REPORT(encounter_prob_i);
  REPORT(log_positive_mean_i);
  REPORT(mu_hat_i);
  REPORT(cpue_density);
  REPORT(mu_index);
  ADREPORT(link_total);
  ADREPORT(link_index);
  return nll;
}
