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
#' - Structural modules that are off should use length-0 vectors or 0x0
#'   matrices, rather than `map`, so they do not enter the AD graph.
#' - Matrices are stored in column-major order when passed to TMB; keep the same
#'   dimension convention as the C++ template expects.
#'
#' @param n_t Number of time steps (e.g., years/quarters) after 0-based recoding.
#' @param n_f Number of flags after recoding (including baseline flag = level 0).
#' @param n_s Number of SPDE vertices (mesh nodes).
#' @param K_smooth_catch Number of fixed-effect columns from catchability smooths.
#' @param n_smooth_catch Number of catchability smooth terms (length of Zs list).
#' @param sum_k_catch Total number of penalized coefficients across catchability smooths.
#' @param K_smooth_pop Number of fixed-effect columns from population smooths.
#' @param n_smooth_pop Number of population smooth terms (length of Zs list).
#' @param sum_k_pop Total number of penalized coefficients across population smooths.
#' @param K_factor_catch Number of fixed-effect columns from catchability factor terms.
#' @param n_factor_random_catch Number of catchability random-factor terms.
#' @param sum_k_factor_catch Total number of random-factor coefficients for catchability.
#' @param K_factor_pop Number of fixed-effect columns from population factor terms.
#' @param n_factor_random_pop Number of population random-factor terms.
#' @param sum_k_factor_pop Total number of random-factor coefficients for population.
#' @param pop_spatial_on Logical; include the omega spatial population field.
#' @param pop_spatiotemporal_on Logical; include the epsilon spatiotemporal
#'   population field.
#' @param K_smooth Legacy alias for `K_smooth_catch`.
#' @param n_smooth Legacy alias for `n_smooth_catch`.
#' @param sum_k Legacy alias for `sum_k_catch`.
#'
#' @return Named list of initial parameter values for TMB::MakeADFun().
#' @keywords internal
.make_parameters_intCPUE <- function(
    n_t,
    n_v,
    n_f,
    n_s,
    data_type = c("delta", "positive", "mixed"),
    n_v_component1 = n_v,
    n_f_component1 = max(0L, n_f - 1L),
    n_flag_t_free_1 = 0L,
    n_flag_t_free_2 = 0L,
    n_flag_t_free = NULL,
    pop_spatial_on = TRUE,
    pop_spatiotemporal_on = TRUE,
    use_q_diffs_time = TRUE,
    q_diffs_spatial_on = TRUE,
    K_smooth_catch = 0L,
    n_smooth_catch = 0L,
    sum_k_catch = 0L,
    K_smooth_pop = 0L,
    n_smooth_pop = 0L,
    sum_k_pop = 0L,
    K_factor_catch = 0L,
    n_factor_random_catch = 0L,
    sum_k_factor_catch = 0L,
    K_factor_pop = 0L,
    n_factor_random_pop = 0L,
    sum_k_factor_pop = 0L,
    split_component1_design = FALSE,
    K_smooth_catch_1 = 0L,
    n_smooth_catch_1 = 0L,
    sum_k_catch_1 = 0L,
    K_factor_catch_1 = 0L,
    n_factor_random_catch_1 = 0L,
    sum_k_factor_catch_1 = 0L,
    K_smooth_pop_1 = 0L,
    n_smooth_pop_1 = 0L,
    sum_k_pop_1 = 0L,
    K_factor_pop_1 = 0L,
    n_factor_random_pop_1 = 0L,
    sum_k_factor_pop_1 = 0L,
    K_smooth = NULL,
    n_smooth = NULL,
    sum_k = NULL
) {
  data_type <- match.arg(data_type)
  if (!is.null(K_smooth)) K_smooth_catch <- K_smooth
  if (!is.null(n_smooth)) n_smooth_catch <- n_smooth
  if (!is.null(sum_k)) sum_k_catch <- sum_k
  if (!is.null(n_flag_t_free)) {
    n_flag_t_free_1 <- n_flag_t_free
    n_flag_t_free_2 <- n_flag_t_free
  }
  has_component1 <- !identical(data_type, "positive")
  split_component1_design <- isTRUE(split_component1_design)
  n_comp <- if (has_component1 && !split_component1_design) 2L else 1L
  n_pop_comp <- if (has_component1) 2L else 1L
  pop_spatial_on <- isTRUE(pop_spatial_on)
  pop_spatiotemporal_on <- isTRUE(pop_spatiotemporal_on)
  n_f_component1 <- if (has_component1) as.integer(n_f_component1) else 0L
  q_diffs_spatial_on <- isTRUE(q_diffs_spatial_on)
  use_any_spde_1 <- has_component1 && (pop_spatial_on || pop_spatiotemporal_on || (q_diffs_spatial_on && n_f_component1 > 0L))
  use_any_spde_2 <- pop_spatial_on || pop_spatiotemporal_on || (q_diffs_spatial_on && n_f > 1L)
  
  # Helper for allocating numeric matrices on the working scale.
  # Using 0.0 is standard: it corresponds to neutral effects on log/logit scales.
  empty_mat <- function(nr, nc) {
    matrix(0.0, nrow = nr, ncol = nc)
  }
  zero_or_vec <- function(flag, n) {
    if (isTRUE(flag) && n > 0L) rep(0.0, n) else numeric(0)
  }
  zero_or_mat <- function(flag, nr, nc) {
    if (isTRUE(flag) && nr > 0L && nc > 0L) empty_mat(nr, nc) else matrix(0.0, nrow = 0L, ncol = 0L)
  }
  comp1_vec <- function(n) {
    if (has_component1 && n > 0L) rep(0.0, n) else numeric(0)
  }
  comp1_mat <- function(nr, nc) {
    if (has_component1 && nr > 0L && nc > 0L) empty_mat(nr, nc) else matrix(0.0, nrow = 0L, ncol = 0L)
  }
  flag_t_sd <- function(component_active) {
    if (isTRUE(component_active)) 0.0 else numeric(0)
  }
  
  # ------------------------------------------------------------
  # IMPORTANT: parameter names below must match C++ PARAMETER()'s.
  # ------------------------------------------------------------
  params <- list(
    # ==========================================================
    # Observation model parameters
    # ==========================================================
    # Residual SD (often used for positive-component likelihood or observation noise).
    # Stored on log scale in C++ for positivity.
    ln_sd = 0.0,
    ln_sd_flag = rep(0.0, n_f),
    
    # ==========================================================
    # Anisotropy parameters
    # ==========================================================
    # ln_H_input typically parameterizes an anisotropy transform H (2 params).
    ln_H_input = if (use_any_spde_1 || use_any_spde_2) c(0.0, 0.0) else numeric(0),
    
    # ==========================================================
    # SPDE hyperparameters (two components, e.g. encounter vs positive)
    # ==========================================================
    # Range parameters (log scale). Often shared across spatial + spatiotemporal fields
    # within a component; the C++ defines how they enter Q.
    ln_range_1 = rep(0.0, as.integer(use_any_spde_1)),
    # Marginal SDs on log scale:
    ln_sigma_0_1 = rep(0.0, as.integer(has_component1 && pop_spatial_on)),    # spatial (omega) SD
    ln_sigma_t_1 = rep(0.0, as.integer(has_component1 && pop_spatiotemporal_on)),  # spatiotemporal (epsilon) SD
    transf_rho_1 = rep(0.0, as.integer(has_component1 && pop_spatiotemporal_on)),  # spatiotemporal AR1 correlation (working scale)
    
    ln_range_2 = rep(0.0, as.integer(use_any_spde_2)),
    ln_sigma_0_2 = rep(0.0, as.integer(pop_spatial_on)),
    ln_sigma_t_2 = rep(0.0, as.integer(pop_spatiotemporal_on)),
    transf_rho_2 = rep(0.0, as.integer(pop_spatiotemporal_on)),

    # ==========================================================
    # Vessel random effects (core catchability component)
    # ==========================================================
    ves_v_1 = if (has_component1 && n_v_component1 > 1L) rep(0.0, n_v_component1) else numeric(0),
    ves_v_2 = rep(0.0, n_v),
    ves_ln_std_dev_1 = if (has_component1 && n_v_component1 > 1L) 0.0 else numeric(0),
    ves_ln_std_dev_2 = 0.0,

    # ==========================================================
    # Default fixed time effects (core population component)
    # ==========================================================
    yq_t_1 = comp1_vec(n_t),
    yq_t_2 = rep(0.0, n_t),
    
    # ==========================================================
    # Spatial + spatiotemporal random fields (SPDE vertices)
    # ==========================================================
    # omega_s_* : spatial random field at mesh vertices (length n_s).
    omega_s_1 = if (pop_spatial_on) comp1_vec(n_s) else numeric(0),
    omega_s_2 = if (pop_spatial_on) rep(0.0, n_s) else numeric(0),
    
    # epsilon_st_* : spatiotemporal random field at mesh vertices across time.
    # Dimension convention here is (n_s x n_t). This must match C++.
    epsilon_st_1 = if (pop_spatiotemporal_on) comp1_mat(n_s, n_t) else matrix(0.0, nrow = 0L, ncol = 0L),
    epsilon_st_2 = if (pop_spatiotemporal_on) empty_mat(n_s, n_t) else matrix(0.0, nrow = 0L, ncol = 0L),
    
    # ==========================================================
    # Flag systematic differences (relative to baseline flag = 0)
    # ==========================================================
    # flag_f_* : one coefficient per non-baseline flag (length n_f-1).
    # If n_f == 1, store as length-0 numeric to keep shapes consistent.
    flag_f_1 = if (has_component1 && n_f_component1 > 0L) rep(0.0, n_f_component1) else numeric(0),
    flag_f_2 = if (n_f > 1L) rep(0.0, n_f - 1L) else numeric(0),
    
    # Log SD for flag systematic differences (if modeled as random effects).
    flag_ln_std_dev_1 = rep(0.0, as.integer(has_component1 && n_f_component1 > 0L)),
    flag_ln_std_dev_2 = 0.0,
    
    # ==========================================================
    # Flag temporal differences (relative to baseline flag = 0)
    # ==========================================================
    # flag_t_* : free coefficients only for observed, non-reference flag-time cells.
    # Missing cells and the first observed time per flag are reconstructed in C++ as 0.
    flag_t_1 = zero_or_vec(has_component1 && use_q_diffs_time, n_flag_t_free_1),
    flag_t_2 = zero_or_vec(use_q_diffs_time, n_flag_t_free_2),
    
    # Log SD for flag temporal differences (per component).
    # If no flags are estimable, fit() may map these to NA to fix variance at 0.
    flag_t_ln_std_dev_1 = flag_t_sd(has_component1 && isTRUE(use_q_diffs_time) && n_flag_t_free_1 > 0L),
    flag_t_ln_std_dev_2 = flag_t_sd(isTRUE(use_q_diffs_time) && n_flag_t_free_2 > 0L),
    
    # ==========================================================
    # Flag spatial differences (relative to baseline flag = 0)
    # ==========================================================
    # flag_s_* : (n_s x (n_f-1)) matrix. Each column is a flag-specific spatial field.
    # SD handled via ln_sigma_flag_* above (component-specific).
    flag_s_1 = zero_or_mat(has_component1 && q_diffs_spatial_on, n_s, n_f_component1),
    flag_s_2 = zero_or_mat(q_diffs_spatial_on, n_s, max(0L, n_f - 1L)),
    
    # flag spatial-diff field SD
    ln_sigma_flag_1 = if (has_component1 && q_diffs_spatial_on && n_f_component1 > 0L) 0.0 else numeric(0),
    ln_sigma_flag_2 = if (q_diffs_spatial_on && n_f > 1L) 0.0 else numeric(0),
    
    # ==========================================================
    # Catchability smooth terms (mgcv::s())
    # ==========================================================
    bs_catch = empty_mat(K_smooth_catch, n_comp),
    b_smooth_catch = empty_mat(sum_k_catch, n_comp),
    ln_smooth_sigma_catch = empty_mat(n_smooth_catch, n_comp),
    bs_catch_1 = if (split_component1_design) rep(0.0, K_smooth_catch_1) else numeric(0),
    b_smooth_catch_1 = if (split_component1_design) rep(0.0, sum_k_catch_1) else numeric(0),
    ln_smooth_sigma_catch_1 = if (split_component1_design) rep(0.0, n_smooth_catch_1) else numeric(0),

    # ==========================================================
    # Catchability factor terms
    # ==========================================================
    bf_catch = empty_mat(K_factor_catch, n_comp),
    b_factor_catch = empty_mat(sum_k_factor_catch, n_comp),
    ln_factor_sigma_catch = empty_mat(n_factor_random_catch, n_comp),
    bf_catch_1 = if (split_component1_design) rep(0.0, K_factor_catch_1) else numeric(0),
    b_factor_catch_1 = if (split_component1_design) rep(0.0, sum_k_factor_catch_1) else numeric(0),
    ln_factor_sigma_catch_1 = if (split_component1_design) rep(0.0, n_factor_random_catch_1) else numeric(0),

    # ==========================================================
    # Population smooth terms (mgcv::s())
    # ==========================================================
    bs_pop = empty_mat(K_smooth_pop, n_comp),
    b_smooth_pop = empty_mat(sum_k_pop, n_comp),
    ln_smooth_sigma_pop = empty_mat(n_smooth_pop, n_comp),
    bs_pop_1 = if (split_component1_design) rep(0.0, K_smooth_pop_1) else numeric(0),
    b_smooth_pop_1 = if (split_component1_design) rep(0.0, sum_k_pop_1) else numeric(0),
    ln_smooth_sigma_pop_1 = if (split_component1_design) rep(0.0, n_smooth_pop_1) else numeric(0),

    # ==========================================================
    # Population factor terms
    # ==========================================================
    pop_intercept = rep(0.0, n_pop_comp),
    bf_pop = empty_mat(K_factor_pop, n_comp),
    b_factor_pop = empty_mat(sum_k_factor_pop, n_comp),
    ln_factor_sigma_pop = empty_mat(n_factor_random_pop, n_comp),
    bf_pop_1 = if (split_component1_design) rep(0.0, K_factor_pop_1) else numeric(0),
    b_factor_pop_1 = if (split_component1_design) rep(0.0, sum_k_factor_pop_1) else numeric(0),
    ln_factor_sigma_pop_1 = if (split_component1_design) rep(0.0, n_factor_random_pop_1) else numeric(0),
    
    # ==========================================================
    # Optional "epsilon trick" / stabilization term
    # ==========================================================
    # eps_index is currently empty; if used, it should be length n_t or similar
    # and must match the way C++ consumes it.
    eps_index = numeric(0)
  )

  # These names are controlled by compile-time TMB wrappers.  When a wrapper
  # disables a module, the corresponding PARAMETER declarations do not exist in
  # C++, so the R parameter list must omit them rather than pass length-0
  # placeholders.
  drop_params <- character(0)
  if (!pop_spatial_on && !pop_spatiotemporal_on && !q_diffs_spatial_on) {
    drop_params <- c(drop_params, "ln_H_input", "ln_range_1", "ln_range_2")
  }
  if (!pop_spatial_on) {
    drop_params <- c(
      drop_params,
      "ln_sigma_0_1", "ln_sigma_0_2",
      "omega_s_1", "omega_s_2"
    )
  }
  if (!pop_spatiotemporal_on) {
    drop_params <- c(
      drop_params,
      "ln_sigma_t_1", "ln_sigma_t_2",
      "transf_rho_1", "transf_rho_2",
      "epsilon_st_1", "epsilon_st_2"
    )
  }
  if (!q_diffs_spatial_on) {
    drop_params <- c(
      drop_params,
      "flag_s_1", "flag_s_2",
      "ln_sigma_flag_1", "ln_sigma_flag_2"
    )
  }
  params[unique(drop_params)] <- NULL

  params
}
