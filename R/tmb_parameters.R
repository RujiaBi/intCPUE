#' Internal: Build initial parameter list for the intCPUE TMB model
#'
#' This function defines the *parameter layout contract* between R and the
#' compiled TMB C++ template. The returned list must:
#' - match the exact parameter names used in the C++ code, and
#' - have the correct dimensions for vectors/matrices expected by the template.
#'
#' Notes for developers:
#' - If add/remove/rename a PARAMETER() or change a parameter dimension in C++,
#'   this function needs update accordingly.
#' - Any parameter that is "turned off" via `map` will be fixed at these initial
#'   values (usually 0 on the working scale), so choose defaults intentionally.
#' - Matrices are stored in column-major order when passed to TMB; keep the same
#'   dimension convention as the C++ template expects.
#'
#' @param n_t Number of time steps (e.g., years/quarters) after 0-based recoding.
#' @param n_v Number of vessels (or vessel IDs) after recoding.
#' @param n_f Number of flags after recoding (including baseline flag = level 0).
#' @param n_s Number of SPDE vertices (mesh nodes).
#' @param K_smooth Number of fixed-effect columns from smooth terms (Xs).
#' @param n_smooth Number of smooth terms (length of Zs list).
#' @param sum_k Total number of penalized coefficients across smooths
#'   (sum of ncol(Zs`k`) over k), i.e., length of b_smooth.
#'
#' @return Named list of initial parameter values for TMB::MakeADFun().
#' @keywords internal
.make_parameters_intCPUE <- function(n_t, n_v, n_f, n_s, K_smooth, n_smooth, sum_k) {
  
  # Helper for allocating numeric matrices on the working scale.
  # Using 0.0 is standard: it corresponds to neutral effects on log/logit scales.
  empty_mat <- function(nr, nc) {
    matrix(0.0, nrow = nr, ncol = nc)
  }
  
  # ------------------------------------------------------------
  # IMPORTANT: parameter names below must match C++ PARAMETER()'s.
  # ------------------------------------------------------------
  list(
    # ==========================================================
    # Observation model parameters
    # ==========================================================
    # Residual SD (often used for positive-component likelihood or observation noise).
    # Stored on log scale in C++ for positivity.
    ln_sd = 0.0,
    
    # ==========================================================
    # Anisotropy parameters
    # ==========================================================
    # ln_H_input typically parameterizes an anisotropy transform H (2 params).
    ln_H_input = c(0.0, 0.0),
    
    # ==========================================================
    # SPDE hyperparameters (two components, e.g. encounter vs positive)
    # ==========================================================
    # Range parameters (log scale). Often shared across spatial + spatiotemporal fields
    # within a component; the C++ defines how they enter Q.
    ln_range_1 = 0.0,
    # Marginal SDs on log scale:
    ln_sigma_0_1    = 0.0,   # spatial (omega) SD
    ln_sigma_t_1    = 0.0,   # spatiotemporal (epsilon) SD
    
    ln_range_2 = 0.0,
    ln_sigma_0_2    = 0.0,
    ln_sigma_t_2    = 0.0,
    
    # ==========================================================
    # Vessel random effects (two components)
    # ==========================================================
    # Vessel-level intercept deviations (one per vessel), typically iid N(0, sd^2).
    # If vessel_effect="off", these are fixed at 0 by map in fit().
    ves_v_1 = rep(0.0, n_v),
    ves_v_2 = rep(0.0, n_v),
    # Log SD for vessel effects (one per component).
    ves_ln_std_dev_1 = 0.0,
    ves_ln_std_dev_2 = 0.0,
    
    # ==========================================================
    # Time effects (global, not flag-specific)
    # ==========================================================
    # yq_t_* are time-varying effects common to all flags.
    # NOTE: name suggests "year-quarter"; ensure semantics match data.
    yq_t_1 = rep(0.0, n_t),
    yq_t_2 = rep(0.0, n_t),
    
    # ==========================================================
    # Spatial + spatiotemporal random fields (SPDE vertices)
    # ==========================================================
    # omega_s_* : spatial random field at mesh vertices (length n_s).
    omega_s_1 = rep(0.0, n_s),
    omega_s_2 = rep(0.0, n_s),
    
    # epsilon_st_* : spatiotemporal random field at mesh vertices across time.
    # Dimension convention here is (n_s x n_t). This must match C++.
    epsilon_st_1 = empty_mat(n_s, n_t),
    epsilon_st_2 = empty_mat(n_s, n_t),
    
    # ==========================================================
    # Flag systematic differences (relative to baseline flag = 0)
    # ==========================================================
    # flag_f_* : one coefficient per non-baseline flag (length n_f-1).
    # If n_f == 1, store as length-0 numeric to keep shapes consistent.
    flag_f_1 = if (n_f > 1L) rep(0.0, n_f - 1L) else numeric(0),
    flag_f_2 = if (n_f > 1L) rep(0.0, n_f - 1L) else numeric(0),
    
    # Log SD for flag systematic differences (if modeled as random effects).
    flag_ln_std_dev_1 = 0.0,
    flag_ln_std_dev_2 = 0.0,
    
    # ==========================================================
    # Flag temporal differences (relative to baseline flag = 0)
    # ==========================================================
    # flag_t_* : (n_t x (n_f-1)) matrix. Each column is a flag-specific time series.
    # In fit(), partially fix missing (t, flag) cells to 0 using has_tf.
    flag_t_1 = empty_mat(n_t, max(0L, n_f - 1L)),
    flag_t_2 = empty_mat(n_t, max(0L, n_f - 1L)),
    
    # Log SD for flag temporal differences (per component).
    # If no flags are estimable, fit() may map these to NA to fix variance at 0.
    flag_t_ln_std_dev_1 = 0.0,
    flag_t_ln_std_dev_2 = 0.0,
    
    # ==========================================================
    # Flag spatial differences (relative to baseline flag = 0)
    # ==========================================================
    # flag_s_* : (n_s x (n_f-1)) matrix. Each column is a flag-specific spatial field.
    # SD handled via ln_sigma_flag_* above (component-specific).
    flag_s_1 = empty_mat(n_s, max(0L, n_f - 1L)),
    flag_s_2 = empty_mat(n_s, max(0L, n_f - 1L)),
    
    # flag spatial-diff field SD
    ln_sigma_flag_1 = 0.0,  
    ln_sigma_flag_2 = 0.0,
    
    # ==========================================================
    # Smooth terms (mgcv::s())
    # ==========================================================
    # bs: fixed-effect coefficients for smooth "parametric" parts (Xs),
    # one column per model component (2 columns = encounter + positive).
    bs = empty_mat(K_smooth, 2L),
    
    # b_smooth: penalized coefficients (stacked across smooths; length = sum_k),
    # stored as a (sum_k x 2) matrix for two components.
    b_smooth = empty_mat(sum_k, 2L),
    
    # ln_smooth_sigma: log SD hyperparameters for each smooth (n_smooth x 2).
    ln_smooth_sigma = empty_mat(n_smooth, 2L),
    
    # ==========================================================
    # Optional "epsilon trick" / stabilization term
    # ==========================================================
    # eps_index is currently empty; if used, it should be length n_t or similar
    # and must match the way C++ consumes it.
    eps_index = numeric(0)
  )
}
