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
#' @param formula Legacy one-sided or two-sided model formula with optional
#'   mgcv::s() smooth terms. If supplied, it is treated as
#'   `formula_catchability`.
#' @param data_utm A data.frame with required columns.
#' @param mesh An `intCPUEmesh` or a bare fmesher mesh.
#' @param formula_catchability Optional one-sided or two-sided formula defining
#'   smooth terms that affect catchability only.
#' @param formula_population Optional one-sided or two-sided formula defining
#'   smooth terms that affect the latent population surface and therefore also
#'   enter projection.
#' @param projection_data Optional projection-grid data for population smooth
#'   covariates. Passed to [make_data()]. For static population covariates,
#'   supply one row per extrapolation grid cell with columns `utm_x_scale`,
#'   `utm_y_scale`, and the covariate names used in `formula_population`. For
#'   time-varying population covariates, additionally include `tid` and supply
#'   one row per grid cell-time combination. If `formula_population` mixes
#'   static and time-varying covariates, use the grid cell-time format and
#'   repeat static covariate values across `tid`.
#' @param q_diffs_time "on" or "off". Controls whether flag-specific temporal
#'   deviations are included. Implemented via dedicated TMB templates rather
#'   than `map`.
#' @param q_diffs_spatial "on" or "off". Controls whether flag-specific spatial
#'   deviations are included. Implemented via dedicated TMB templates rather
#'   than `map`.
#' @param pop_spatiotemporal_type `"rw"` or `"ar1"`. Controls whether the
#'   always-on spatiotemporal population field follows a random-walk or AR1
#'   evolution over time.
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
#' @author Rujia Bi \email{rbi@@iattc.org}
#' @export
intCPUE <- function(
    formula = NULL,
    data_utm,
    mesh,
    formula_catchability = NULL,
    formula_population = NULL,
    projection_data = NULL,
    q_diffs_time = c("on", "off"),
    q_diffs_spatial = c("on", "off"),
    pop_spatiotemporal_type = c("rw", "ar1"),
    obs_sd = c("shared", "flag"),
    control = list(eval.max = 1e5, iter.max = 1e5),
    ncores = NULL,
    ...,
    silent = FALSE,
    restart_max = 1L,
    newton_max = 2L,
    coord_max = 5L
) {
  q_diffs_time <- match.arg(q_diffs_time)
  q_diffs_spatial <- match.arg(q_diffs_spatial)
  pop_spatiotemporal_type <- match.arg(pop_spatiotemporal_type)
  obs_sd <- match.arg(obs_sd)
  restart_max <- .validate_nonneg_count(restart_max, "restart_max")
  newton_max <- .validate_nonneg_count(newton_max, "newton_max")
  coord_max <- .validate_nonneg_count(coord_max, "coord_max")
  
  data_utm <- as.data.frame(data_utm)
  
  # ---- 1) Data prep (single source of truth) ----
  prep <- make_data(
    formula = formula,
    data_utm = data_utm,
    mesh = mesh,
    formula_catchability = formula_catchability,
    formula_population = formula_population,
    projection_data = projection_data,
    ...
  )
  
  data_tmb <- prep$data
  
  # ---- 2) Defensive checks (catch mismatches early) ----
  n_t <- data_tmb$n_t
  n_v <- data_tmb$n_v
  n_f <- data_tmb$n_f
  n_s <- data_tmb$spde$n_s

  data_tmb$use_q_diffs_time <- as.integer(q_diffs_time == "on" && n_f > 1L)
  data_tmb$use_q_diffs_spatial <- as.integer(q_diffs_spatial == "on" && n_f > 1L)
  data_tmb$use_pop_spatiotemporal_rw <- as.integer(pop_spatiotemporal_type == "rw")
  data_tmb$use_pop_spatiotemporal_ar1 <- as.integer(pop_spatiotemporal_type == "ar1")
  data_tmb$use_flag_sd <- as.integer(obs_sd == "flag" && n_f > 0L)

  has_tf <- NULL
  flag_t_constraint <- NULL
  if (q_diffs_time == "on" && n_f > 1L) {
    has_tf <- data_tmb$has_tf > 0L
    flag_t_constraint <- .build_flag_t_constraint(has_tf)
    data_tmb$flag_t_index <- matrix(
      as.integer(flag_t_constraint$flag_t_index),
      nrow = nrow(flag_t_constraint$flag_t_index),
      ncol = ncol(flag_t_constraint$flag_t_index)
    )
  } else {
    data_tmb$flag_t_index <- matrix(integer(0), nrow = 0L, ncol = max(0L, n_f - 1L))
  }
  
  # Smooth dims
  has_smooths_catch <- isTRUE(data_tmb$has_smooths_catch == 1L)
  K_smooth_catch <- if (has_smooths_catch) ncol(data_tmb$Xs_catch) else 0L
  n_smooth_catch <- if (has_smooths_catch) length(data_tmb$Zs_catch) else 0L
  sum_k_catch <- if (has_smooths_catch && n_smooth_catch > 0L) sum(vapply(data_tmb$Zs_catch, ncol, 0L)) else 0L

  has_smooths_pop <- isTRUE(data_tmb$has_smooths_pop == 1L)
  K_smooth_pop <- if (has_smooths_pop) ncol(data_tmb$Xs_pop_i) else 0L
  n_smooth_pop <- if (has_smooths_pop) length(data_tmb$Zs_pop_i) else 0L
  sum_k_pop <- if (has_smooths_pop && n_smooth_pop > 0L) sum(vapply(data_tmb$Zs_pop_i, ncol, 0L)) else 0L
  
  # ---- 3) Initial parameters (must match cpp) ----
  parameters <- .make_parameters_intCPUE(
    n_t = n_t, n_v = n_v, n_f = n_f, n_s = n_s,
    n_flag_t_free = if (is.null(flag_t_constraint)) 0L else flag_t_constraint$n_free,
    use_q_diffs_time = (q_diffs_time == "on" && n_f > 1L),
    use_q_diffs_spatial = (q_diffs_spatial == "on" && n_f > 1L),
    K_smooth_catch = K_smooth_catch,
    n_smooth_catch = n_smooth_catch,
    sum_k_catch = sum_k_catch,
    K_smooth_pop = K_smooth_pop,
    n_smooth_pop = n_smooth_pop,
    sum_k_pop = sum_k_pop
  )

  # Match the active TMB template exactly: q-time/q-spatial parameters do not
  # exist at all in the corresponding off-templates.
  if (q_diffs_time == "off" || n_f <= 1L) {
    parameters$flag_t_1 <- NULL
    parameters$flag_t_2 <- NULL
    parameters$flag_t_ln_std_dev_1 <- NULL
    parameters$flag_t_ln_std_dev_2 <- NULL
  }
  if (q_diffs_spatial == "off" || n_f <= 1L) {
    parameters$flag_s_1 <- NULL
    parameters$flag_s_2 <- NULL
    parameters$ln_sigma_flag_1 <- NULL
    parameters$ln_sigma_flag_2 <- NULL
  }
  
  # ---- 4) MAP (turn on/off components without touching cpp) ----
  map <- .make_map_intCPUE(
    parameters = parameters,
    n_f = n_f,
    obs_sd = obs_sd,
    pop_spatiotemporal_type = pop_spatiotemporal_type,
    q_diffs_time = q_diffs_time,
    has_tf = has_tf,
    estimable_flag = if (is.null(flag_t_constraint)) NULL else flag_t_constraint$estimable_flag
  )
  
  # ---- 5) Random effects list ----
  random <- character(0)

  random <- c(random, "omega_s_1", "omega_s_2")
  random <- c(random, "epsilon_st_1", "epsilon_st_2")

  if (n_v > 0L) {
    random <- c(random, "ves_v_1", "ves_v_2")
  }
  
  if (n_f > 1L) {
    random <- c(random, "flag_f_1", "flag_f_2")
  }
  if (n_f > 1L && q_diffs_time == "on") {
    random <- c(random, "flag_t_1", "flag_t_2")
  }
  if (n_f > 1L && q_diffs_spatial == "on") {
    random <- c(random, "flag_s_1", "flag_s_2")
  }
  
  # Smooth random coeffs
  if (has_smooths_catch && sum_k_catch > 0L) random <- c(random, "b_smooth_catch")
  if (has_smooths_pop && sum_k_pop > 0L) random <- c(random, "b_smooth_pop")
  
  random <- unique(random)
  
  # ---- 6) Build & optimize ----
  .check_fit_inputs_intCPUE(data_tmb)
  
  DLL <- .ensure_template_dll_intCPUE(
    q_diffs_time = q_diffs_time,
    q_diffs_spatial = q_diffs_spatial
  )
  
  if (!is.null(ncores)) {
    ncores <- as.integer(ncores)
    if (is.na(ncores) || ncores < 1L) {
      stop("`ncores` must be a positive integer.", call. = FALSE)
    }
    if (ncores > 1L) {
      TMB::openmp(ncores, autopar = TRUE, DLL = DLL)
    }
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
      formula = formula,
      formula_catchability = if (is.null(formula_catchability)) formula else formula_catchability,
      formula_population = formula_population,
      pop_spatiotemporal_type = pop_spatiotemporal_type,
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
