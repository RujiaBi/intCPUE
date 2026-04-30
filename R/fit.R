#' @useDynLib intCPUE, .registration = TRUE
NULL

#' Fit an integrated spatiotemporal CPUE standardization model
#'
#' Fits an integrated spatiotemporal model for CPUE standardization that jointly
#' models observation and sampling processes across one or more fisheries or surveys.
#' The framework is designed to reduce bias caused by preferential sampling,
#' targeting behavior, and heterogeneous effort distributions, and to estimate
#' a coherent latent abundance index.
#'
#' The model is implemented using Template Model Builder (TMB). Spatial and
#' spatiotemporal random fields are represented using the SPDE (stochastic partial
#' differential equation) approach for computational efficiency. See the model
#' description vignette for details:
#' https://github.com/RujiaBi/intCPUE/blob/main/vignettes/intCPUE-intro.Rmd
#'
#' @param data_utm A data.frame with required columns.
#' @param mesh An `intCPUEmesh` or a bare fmesher mesh.
#' @param formula_catchability Optional one-sided formula defining
#'   smooth and factor terms that affect catchability only.
#' @param formula_population Optional one-sided formula defining
#'   smooth and factor terms that affect the latent population surface and
#'   therefore also enter projection.
#' @param projection_data Optional projection-grid data for population smooth
#'   covariates. Passed to [make_data()]. For static population covariates,
#'   supply one row per extrapolation grid cell with columns `utm_x_scale`,
#'   `utm_y_scale`, and the covariate names used in `formula_population`. For
#'   time-varying population covariates, additionally include `tid` and supply
#'   one row per grid cell-time combination. If `formula_population` mixes
#'   static and time-varying covariates, use the grid cell-time format and
#'   repeat static covariate values across `tid`.
#' @param data_type Observation-data structure. `"auto"` detects the structure
#'   by flag/type; `"delta"` requires every flag/type to have both zero and
#'   positive catches; `"positive"` requires all observations to be positive;
#'   `"mixed"` requires some delta flags and some positive-only flags. For
#'   mixed data, the reference flag (`flagid = 0`) must be a delta flag.
#' @param delta_link `"Poisson"` or `"traditional"`. `"Poisson"` is only valid
#'   for pure `data_type = "delta"`. `"traditional"` is required for
#'   `data_type = "mixed"` and uses a logit encounter probability with a
#'   lognormal positive-CPUE component. For `data_type = "positive"`,
#'   `delta_link` is ignored because only the positive lognormal component is
#'   fitted.
#' @param pop_spatial `"on"` or `"off"`. Controls the time-constant spatial
#'   population field (`omega`). If `"off"`, the corresponding parameters have
#'   length 0 and the C++ template skips the module.
#' @param pop_spatiotemporal `"on"` or `"off"`. Controls the spatiotemporal
#'   population field (`epsilon`). If `"off"`, the corresponding parameters
#'   have length 0 and the C++ template skips the module.
#' @param pop_spatiotemporal_type `"rw"` or `"ar1"`. Controls whether the
#'   spatiotemporal population field follows a random-walk or AR1 evolution
#'   over time.
#' @param q_diffs_system_type `"random"` or `"fixed"`. Controls only the
#'   flag/system main effect (`flag_f`). Flag temporal and spatial deviations
#'   remain random effects when enabled.
#' @param q_diffs_time "on" or "off". Controls whether flag-specific temporal
#'   deviations are included. Free flag-time cells are reduced in `make_data()`
#'   and passed as compact random-effect vectors.
#' @param q_diffs_spatial "on" or "off". Controls whether flag-specific spatial
#'   deviations are included. If off, the corresponding parameters have length
#'   0 and the C++ template skips the module.
#' @param obs_sd `"shared"` or `"flag"`. If `"flag"`, the positive-catch
#'   lognormal observation SD is estimated separately for each flag.
#' @param control Control list passed to [stats::nlminb()].
#' @param ncores Optional integer. If provided, sets the number of OpenMP threads. Passed to [TMB::openmp()].
#' @param ... Passed to [intCPUE::make_data()] (for example, `area_scale`).
#' @param silent Logical. Passed to [TMB::MakeADFun()].
#' @param restart_max Non-negative integer. Maximum number of automatic
#'   `nlminb()` restarts attempted from the current best parameter vector when
#'   convergence or gradient checks remain unsatisfactory.
#' @param newton_max Non-negative integer. Maximum number of
#'   [TMB::newton()] refinement attempts after `nlminb()`/restart stages.
#' @param coord_max Non-negative integer. Maximum number of coordinate-polish
#'   iterations attempted after the restart/Newton stages. Set
#'   `restart_max = 0`, `newton_max = 0`, and `coord_max = 0` to run only a
#'   single `nlminb()` pass.
#'
#' @return An object of class `intCPUE` with elements `obj`, `opt`, `rep`, `prep`, etc.
#'
#' @details
#' The baseline population-density linear predictor is
#' `temporal effect (tid) + omega + epsilon + population factor effects + population smooth effects`.
#' The temporal effect (`tid`) is always present; the spatial random field
#' (`omega`) and spatiotemporal random field (`epsilon`) are controlled by
#' `pop_spatial` and `pop_spatiotemporal`.
#'
#' By default, `tid` is included as a fixed effect. If `formula_population`
#' contains `f(tid)`, the default fixed `tid` effect is replaced by
#' `intercept + random tid effect`.
#'
#' Vessel effects are part of the core catchability model and are always included.
#' If `n_v = 1`, the vessel effect is fixed to 0.
#'
#' Population-layer smooths and factor effects enter both the observation model
#' and projection. Catchability-layer smooths and factor effects enter the
#' observation model only.
#'
#' For `data_type = "mixed"`, component 1 is estimated only from delta flags.
#' Positive-only flags have encounter probability fixed at 1 and contribute only
#' the positive-catch likelihood. Component-1 vessel, flag, smooth, and factor
#' structures are compressed in [make_data()] so positive-only flags do not
#' introduce component-1 vessel effects, flag effects, flag temporal/spatial
#' effects, or unsupported component-1 smooth/factor terms. Component 2 keeps
#' the full positive-catch design.
#' @author Rujia Bi \email{rbi@@iattc.org}
#' @export
intCPUE <- function(
    data_utm,
    mesh,
    formula_catchability = NULL,
    formula_population = NULL,
    projection_data = NULL,
    data_type = c("auto", "delta", "positive", "mixed"),
    delta_link = c("Poisson", "traditional"),
    pop_spatial = c("on", "off"),
    pop_spatiotemporal = c("on", "off"),
    pop_spatiotemporal_type = c("rw", "ar1"),
    q_diffs_system_type = c("random", "fixed"),
    q_diffs_time = c("on", "off"),
    q_diffs_spatial = c("on", "off"),
    obs_sd = c("shared", "flag"),
    control = list(eval.max = 1e5, iter.max = 1e5),
    ncores = NULL,
    ...,
    silent = FALSE,
    restart_max = 1L,
    newton_max = 2L,
    coord_max = 5L
) {
  delta_link_was_missing <- missing(delta_link)
  data_type <- match.arg(data_type)
  delta_link <- match.arg(delta_link)
  q_diffs_system_type <- match.arg(q_diffs_system_type)
  pop_spatial <- match.arg(pop_spatial)
  pop_spatiotemporal <- match.arg(pop_spatiotemporal)
  q_diffs_time <- match.arg(q_diffs_time)
  q_diffs_spatial <- match.arg(q_diffs_spatial)
  pop_spatiotemporal_type <- match.arg(pop_spatiotemporal_type)
  obs_sd <- match.arg(obs_sd)
  restart_max <- .validate_nonneg_count(restart_max, "restart_max")
  newton_max <- .validate_nonneg_count(newton_max, "newton_max")
  coord_max <- .validate_nonneg_count(coord_max, "coord_max")
  pop_spatial_on <- identical(pop_spatial, "on")
  pop_spatiotemporal_on <- identical(pop_spatiotemporal, "on")
  
  data_utm <- as.data.frame(data_utm)
  data_type_info <- .detect_data_type_intCPUE(data_utm, data_type = data_type)
  data_type <- data_type_info$data_type
  if (delta_link_was_missing && !identical(data_type, "delta")) {
    delta_link <- "traditional"
  }
  delta_link <- .validate_delta_link_intCPUE(delta_link, data_type = data_type)
  
  # ---- 1) Data prep (single source of truth) ----
  prep <- make_data(
    data_utm = data_utm,
    mesh = mesh,
    formula_catchability = formula_catchability,
    formula_population = formula_population,
    projection_data = projection_data,
    flag_uses_component1 = data_type_info$flag_uses_component1,
    ...
  )
  
  data_tmb <- prep$data
  
  # ---- 2) Defensive checks (catch mismatches early) ----
  n_t <- data_tmb$n_t
  n_v <- data_tmb$n_v
  n_f <- data_tmb$n_f
  n_s <- data_tmb$spde$n_s
  use_random_tid_effect <- .factor_has_random_var(prep$smooth_info$population$factor_terms, "tid")
  has_component1 <- isTRUE(data_tmb$has_component1 == 1L)
  q_diffs_spatial_active <- q_diffs_spatial == "on" && n_f > 1L

  if (!identical(data_type_info$flag_id, seq.int(0L, n_f - 1L))) {
    stop("Internal data-type flag IDs do not match recoded `flagid`; re-run with contiguous 0-based flag IDs.", call. = FALSE)
  }
  data_tmb$use_poisson_link <- as.integer(identical(delta_link, "Poisson"))
  data_tmb$use_flag_f_random <- as.integer(identical(q_diffs_system_type, "random"))
  data_tmb$use_q_diffs_time <- as.integer(q_diffs_time == "on" && n_f > 1L)
  data_tmb$use_pop_spatiotemporal_rw <- as.integer(pop_spatiotemporal_type == "rw")
  data_tmb$use_pop_spatiotemporal_ar1 <- as.integer(pop_spatiotemporal_type == "ar1")
  data_tmb$use_flag_sd <- as.integer(obs_sd == "flag" && n_f > 0L)

  split_component1_design <- isTRUE(data_tmb$use_split_component1_design == 1L)
  n_v_component1 <- sum(data_tmb$v_component1_keep == 1L)
  n_f_component1 <- sum(data_tmb$flag_component1_keep == 1L)
  has_tf <- if (q_diffs_time == "on" && n_f > 1L) data_tmb$has_tf > 0L else NULL
  
  # Smooth dims
  has_smooths_catch <- isTRUE(data_tmb$has_smooths_catch == 1L)
  K_smooth_catch <- if (has_smooths_catch) ncol(data_tmb$Xs_catch) else 0L
  n_smooth_catch <- if (has_smooths_catch) length(data_tmb$Zs_catch) else 0L
  sum_k_catch <- if (has_smooths_catch && n_smooth_catch > 0L) sum(vapply(data_tmb$Zs_catch, ncol, 0L)) else 0L

  has_smooths_pop <- isTRUE(data_tmb$has_smooths_pop == 1L)
  K_smooth_pop <- if (has_smooths_pop) ncol(data_tmb$Xs_pop_i) else 0L
  n_smooth_pop <- if (has_smooths_pop) length(data_tmb$Zs_pop_i) else 0L
  sum_k_pop <- if (has_smooths_pop && n_smooth_pop > 0L) sum(vapply(data_tmb$Zs_pop_i, ncol, 0L)) else 0L

  has_fixed_factors_catch <- isTRUE(data_tmb$has_fixed_factors_catch == 1L)
  K_factor_catch <- if (has_fixed_factors_catch) ncol(data_tmb$Xf_catch) else 0L
  has_random_factors_catch <- isTRUE(data_tmb$has_random_factors_catch == 1L)
  n_factor_random_catch <- if (has_random_factors_catch) length(data_tmb$Zf_catch) else 0L
  sum_k_factor_catch <- if (has_random_factors_catch && n_factor_random_catch > 0L) sum(vapply(data_tmb$Zf_catch, ncol, 0L)) else 0L

  has_fixed_factors_pop <- isTRUE(data_tmb$has_fixed_factors_pop == 1L)
  K_factor_pop <- if (has_fixed_factors_pop) ncol(data_tmb$Xf_pop_i) else 0L
  has_random_factors_pop <- isTRUE(data_tmb$has_random_factors_pop == 1L)
  n_factor_random_pop <- if (has_random_factors_pop) length(data_tmb$Zf_pop_i) else 0L
  sum_k_factor_pop <- if (has_random_factors_pop && n_factor_random_pop > 0L) sum(vapply(data_tmb$Zf_pop_i, ncol, 0L)) else 0L

  K_smooth_catch_1 <- if (split_component1_design && isTRUE(data_tmb$has_smooths_catch_1 == 1L)) ncol(data_tmb$Xs_catch_1) else 0L
  n_smooth_catch_1 <- if (split_component1_design && isTRUE(data_tmb$has_smooths_catch_1 == 1L)) length(data_tmb$Zs_catch_1) else 0L
  sum_k_catch_1 <- if (n_smooth_catch_1 > 0L) sum(vapply(data_tmb$Zs_catch_1, ncol, 0L)) else 0L
  K_factor_catch_1 <- if (split_component1_design && isTRUE(data_tmb$has_fixed_factors_catch_1 == 1L)) ncol(data_tmb$Xf_catch_1) else 0L
  n_factor_random_catch_1 <- if (split_component1_design && isTRUE(data_tmb$has_random_factors_catch_1 == 1L)) length(data_tmb$Zf_catch_1) else 0L
  sum_k_factor_catch_1 <- if (n_factor_random_catch_1 > 0L) sum(vapply(data_tmb$Zf_catch_1, ncol, 0L)) else 0L
  K_smooth_pop_1 <- if (split_component1_design && isTRUE(data_tmb$has_smooths_pop_1 == 1L)) ncol(data_tmb$Xs_pop_i_1) else 0L
  n_smooth_pop_1 <- if (split_component1_design && isTRUE(data_tmb$has_smooths_pop_1 == 1L)) length(data_tmb$Zs_pop_i_1) else 0L
  sum_k_pop_1 <- if (n_smooth_pop_1 > 0L) sum(vapply(data_tmb$Zs_pop_i_1, ncol, 0L)) else 0L
  K_factor_pop_1 <- if (split_component1_design && isTRUE(data_tmb$has_fixed_factors_pop_1 == 1L)) ncol(data_tmb$Xf_pop_i_1) else 0L
  n_factor_random_pop_1 <- if (split_component1_design && isTRUE(data_tmb$has_random_factors_pop_1 == 1L)) length(data_tmb$Zf_pop_i_1) else 0L
  sum_k_factor_pop_1 <- if (n_factor_random_pop_1 > 0L) sum(vapply(data_tmb$Zf_pop_i_1, ncol, 0L)) else 0L
  
  # ---- 3) Initial parameters (must match cpp) ----
  parameters <- .make_parameters_intCPUE(
    n_t = n_t, n_v = n_v, n_f = n_f, n_s = n_s,
    data_type = data_type,
    n_v_component1 = n_v_component1,
    n_f_component1 = n_f_component1,
    n_flag_t_free_1 = if (q_diffs_time == "on" && n_f > 1L) data_tmb$n_flag_t_free_1 else 0L,
    n_flag_t_free_2 = if (q_diffs_time == "on" && n_f > 1L) data_tmb$n_flag_t_free_2 else 0L,
    pop_spatial_on = pop_spatial_on,
    pop_spatiotemporal_on = pop_spatiotemporal_on,
    use_q_diffs_time = (q_diffs_time == "on" && n_f > 1L),
    q_diffs_spatial_on = q_diffs_spatial_active,
    K_smooth_catch = K_smooth_catch,
    n_smooth_catch = n_smooth_catch,
    sum_k_catch = sum_k_catch,
    K_smooth_pop = K_smooth_pop,
    n_smooth_pop = n_smooth_pop,
    sum_k_pop = sum_k_pop,
    K_factor_catch = K_factor_catch,
    n_factor_random_catch = n_factor_random_catch,
    sum_k_factor_catch = sum_k_factor_catch,
    K_factor_pop = K_factor_pop,
    n_factor_random_pop = n_factor_random_pop,
    sum_k_factor_pop = sum_k_factor_pop,
    split_component1_design = split_component1_design,
    K_smooth_catch_1 = K_smooth_catch_1,
    n_smooth_catch_1 = n_smooth_catch_1,
    sum_k_catch_1 = sum_k_catch_1,
    K_factor_catch_1 = K_factor_catch_1,
    n_factor_random_catch_1 = n_factor_random_catch_1,
    sum_k_factor_catch_1 = sum_k_factor_catch_1,
    K_smooth_pop_1 = K_smooth_pop_1,
    n_smooth_pop_1 = n_smooth_pop_1,
    sum_k_pop_1 = sum_k_pop_1,
    K_factor_pop_1 = K_factor_pop_1,
    n_factor_random_pop_1 = n_factor_random_pop_1,
    sum_k_factor_pop_1 = sum_k_factor_pop_1
  )

  # ---- 4) MAP (fix non-identifiable scalars; structural modules use length-0 parameters) ----
  map <- .make_map_intCPUE(
    parameters = parameters,
    n_f = n_f,
    n_v = n_v,
    obs_sd = obs_sd,
    pop_spatiotemporal_type = pop_spatiotemporal_type,
    use_random_tid_effect = use_random_tid_effect,
    q_diffs_system_type = q_diffs_system_type,
    q_diffs_time = q_diffs_time,
    has_tf = has_tf,
    estimable_flag = if (q_diffs_time == "on" && n_f > 1L) prep$component1_info$flag_t_constraint_2$estimable_flag else NULL
  )
  
  # ---- 5) Random effects list ----
  random <- character(0)

  if (n_v > 1L) {
    if (has_component1 && length(parameters$ves_v_1) > 0L) random <- c(random, "ves_v_1")
    random <- c(random, "ves_v_2")
  }
  if (pop_spatial_on) {
    if (has_component1 && length(parameters$omega_s_1) > 0L) random <- c(random, "omega_s_1")
    if (length(parameters$omega_s_2) > 0L) random <- c(random, "omega_s_2")
  }
  if (pop_spatiotemporal_on) {
    if (has_component1 && length(parameters$epsilon_st_1) > 0L) random <- c(random, "epsilon_st_1")
    if (length(parameters$epsilon_st_2) > 0L) random <- c(random, "epsilon_st_2")
  }

  if (n_f > 1L && identical(q_diffs_system_type, "random")) {
    if (has_component1 && length(parameters$flag_f_1) > 0L) random <- c(random, "flag_f_1")
    random <- c(random, "flag_f_2")
  }
  if (length(parameters$flag_t_1) > 0L) random <- c(random, "flag_t_1")
  if (length(parameters$flag_t_2) > 0L) random <- c(random, "flag_t_2")
  if (q_diffs_spatial_active) {
    if (has_component1 && length(parameters$flag_s_1) > 0L) random <- c(random, "flag_s_1")
    if (length(parameters$flag_s_2) > 0L) random <- c(random, "flag_s_2")
  }
  
  # Smooth random coeffs
  if (has_smooths_catch && sum_k_catch > 0L) random <- c(random, "b_smooth_catch")
  if (split_component1_design && sum_k_catch_1 > 0L) random <- c(random, "b_smooth_catch_1")
  if (has_smooths_pop && sum_k_pop > 0L) random <- c(random, "b_smooth_pop")
  if (split_component1_design && sum_k_pop_1 > 0L) random <- c(random, "b_smooth_pop_1")
  if (has_random_factors_catch && sum_k_factor_catch > 0L) random <- c(random, "b_factor_catch")
  if (split_component1_design && sum_k_factor_catch_1 > 0L) random <- c(random, "b_factor_catch_1")
  if (has_random_factors_pop && sum_k_factor_pop > 0L) random <- c(random, "b_factor_pop")
  if (split_component1_design && sum_k_factor_pop_1 > 0L) random <- c(random, "b_factor_pop_1")
  
  random <- unique(random)
  
  # ---- 6) Build & optimize ----
  .check_fit_inputs_intCPUE(data_tmb)
  
  DLL <- .ensure_template_dll_intCPUE(
    q_diffs_spatial = if (q_diffs_spatial_active) "on" else "off",
    pop_spatial = pop_spatial,
    pop_spatiotemporal = pop_spatiotemporal
  )
  
  if (!is.null(ncores)) {
    ncores <- as.integer(ncores)
    if (is.na(ncores) || ncores < 1L) {
      stop("`ncores` must be a positive integer.", call. = FALSE)
    }
    TMB::openmp(ncores, autopar = ncores > 1L, DLL = DLL)
  }
  
  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = parameters,
    map = map,
    random = random,
    DLL = DLL,
    silent = silent
  )
  
  opt <- .safe_optimize(
    obj = obj,
    control = control,
    restart_max = restart_max,
    newton_max = newton_max,
    coord_max = coord_max
  )
  par_structured <- try(obj$env$parList(opt$par), silent = TRUE)
  if (!inherits(par_structured, "try-error")) {
    opt$par_list <- par_structured
  }
  
  rep_info <- .safe_sdreport(obj, opt)
  rep <- rep_info$rep
  
  out <- list(
    obj = obj,
    opt = opt,
    rep = rep,
    prep = prep,
    data_tmb = data_tmb,
    map = map,
    random = random,
    control = control,
    settings = list(
      formula_catchability = formula_catchability,
      formula_population = formula_population,
      data_type = data_type,
      data_type_detected = data_type_info$detected,
      flag_type = stats::setNames(data_type_info$flag_type, data_type_info$flag_id),
      delta_link = delta_link,
      pop_spatial = pop_spatial,
      pop_spatiotemporal = pop_spatiotemporal,
      pop_spatiotemporal_type = pop_spatiotemporal_type,
      q_diffs_system_type = q_diffs_system_type,
      q_diffs_time = q_diffs_time,
      q_diffs_spatial = q_diffs_spatial,
      obs_sd = obs_sd,
      DLL = DLL,
      ncores = ncores,
      restart_max = restart_max,
      newton_max = newton_max,
      coord_max = coord_max
    ),
    diagnostics = list(
      convergence = opt$convergence,
      message = opt$message,
      max_grad = if (!is.null(opt$max_grad)) opt$max_grad else max(abs(obj$gr(opt$par))),
      newton_used = if (!is.null(opt$newton_used)) opt$newton_used else NA_integer_,
      coord_used = if (!is.null(opt$coord_used)) opt$coord_used else NA_integer_,
      sdreport_method = rep_info$method,
      sdreport_error = rep_info$error
    )
  )
  class(out) <- "intCPUE"
  out
}
