#' @useDynLib intCPUE, .registration = TRUE
NULL

#' Fits a spatiotemporal generalized linear mixed-effects model (GLMM)
#' using the TMB (Template Model Builder) R package. Spatial and spatiotemporal
#' random fields are approximated using the SPDE (stochastic partial differential
#' equation) approach, allowing efficient modeling of data that are correlated
#' in space and/or time.
#'
#' @param formula A model formula with optional mgcv::s() smooth terms.
#' @param data_utm A data.frame with required columns.
#' @param mesh An `intCPUEmesh` or a bare fmesher mesh.
#' @param vessel_effect "on" or "off". If "off", vessel RE is fixed at 0 via map.
#' @param q_diffs_system "on" or "off". If "off", flag systematic differences are fixed.
#' @param q_diffs_time "on" or "off". If "off", flag-specific temporal effects are fixed.
#' @param q_diffs_spatial "on" or "off". If "off", flag-specific spatial fields are fixed.
#' @param control Control list passed to [stats::nlminb()].
#' @param ncores Optional integer. If provided, sets the number of OpenMP threads. Passed to [TMB::openmp()].
#' @param silent Logical. Passed to [TMB::MakeADFun()].
#' @param ... Passed to [intCPUE::make_data()] (e.g., utm_zone, coord_scale, area_scale).
#'
#' @return An object of class `intCPUE` with elements `obj`, `opt`, `rep`, `prep`, etc.
#' @author Rujia Bi \email{rbi@@iattc.org}
#' @export
intCPUE <- function(
    formula,
    data_utm,
    mesh,
    vessel_effect = c("on", "off"),
    q_diffs_system = c("on", "off"),
    q_diffs_time = c("on", "off"),
    q_diffs_spatial = c("on", "off"),
    control = list(eval.max = 1e5, iter.max = 1e5),
    ncores = NULL,
    ...,
    silent = FALSE
) {
  vessel_effect   <- match.arg(vessel_effect)
  q_diffs_system  <- match.arg(q_diffs_system)
  q_diffs_time    <- match.arg(q_diffs_time)
  q_diffs_spatial <- match.arg(q_diffs_spatial)
  
  data_utm <- as.data.frame(data_utm)
  
  # ---- 1) Data prep (single source of truth) ----
  prep <- make_data(
    formula = formula,
    data_utm = data_utm,
    mesh = mesh,
    ...
  )
  
  data_tmb <- prep$data
  
  # ---- 2) Defensive checks (catch mismatches early) ----
  n_t <- data_tmb$n_t
  n_v <- data_tmb$n_v
  n_f <- data_tmb$n_f
  n_s <- data_tmb$spde$n_s
  
  # Smooth dims
  has_smooths <- isTRUE(data_tmb$has_smooths == 1L)
  K_smooth <- if (has_smooths) ncol(data_tmb$Xs) else 0L
  n_smooth <- if (has_smooths) length(data_tmb$Zs) else 0L
  sum_k <- if (has_smooths && n_smooth > 0L) sum(vapply(data_tmb$Zs, ncol, 0L)) else 0L
  
  # ---- 3) Initial parameters (must match cpp) ----
  parameters <- .make_parameters_intCPUE(
    n_t = n_t, n_v = n_v, n_f = n_f, n_s = n_s,
    K_smooth = K_smooth, n_smooth = n_smooth, sum_k = sum_k
  )
  
  # ---- 4) MAP (turn on/off components without touching cpp) ----
  has_tf <- NULL
  if (!is.null(data_tmb$has_tf)) {
    # has_tf is integer matrix 0/1 in your make_data()
    has_tf <- (data_tmb$has_tf > 0L)
  }
  
  map <- .make_map_intCPUE(
    parameters = parameters,
    n_f = n_f,
    vessel_effect = vessel_effect,
    q_diffs_system = q_diffs_system,
    q_diffs_time = q_diffs_time,
    q_diffs_spatial = q_diffs_spatial,
    has_tf = has_tf
  )
  
  # ---- 5) Random effects list ----
  random <- c(
    "omega_s_1", "epsilon_st_1",
    "omega_s_2", "epsilon_st_2"
  )
  
  if (vessel_effect == "on") {
    random <- c(random, "ves_v_1", "ves_v_2")
  }
  
  # Only include these if flags exist:
  if (n_f > 1L && q_diffs_system == "on") {
    random <- c(random, "flag_f_1", "flag_f_2")
  }
  if (n_f > 1L && q_diffs_time == "on") {
    random <- c(random, "flag_t_1", "flag_t_2")
  }
  if (n_f > 1L && q_diffs_spatial == "on") {
    random <- c(random, "flag_s_1", "flag_s_2")
  }
  
  # Smooth random coeffs
  if (has_smooths && sum_k > 0L) random <- c(random, "b_smooth")
  
  random <- unique(random)
  
  # ---- 6) Build & optimize ----
  .check_fit_inputs_intCPUE(data_tmb)
  
  DLL <- "intCPUE"
  
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
  
  opt <- .safe_optimize(obj, control)
  
  rep <- TMB::sdreport(obj)
  
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
      vessel_effect = vessel_effect,
      q_diffs_system = q_diffs_system,
      q_diffs_time = q_diffs_time,
      q_diffs_spatial = q_diffs_spatial,
      DLL = DLL,
      ncores = ncores
    ),
    diagnostics = list(
      convergence = opt$convergence,
      message = opt$message,
      max_grad = max(abs(obj$gr(opt$par)))
    )
  )
  class(out) <- "intCPUE"
  out
}
