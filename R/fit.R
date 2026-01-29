#' Fit intCPUE model (encounter + positive)
#'
#' @param formula A model formula with optional mgcv::s() smooth terms.
#'   (For now, formula is used only to build smooth objects; other fixed effects
#'   are not yet parsed in this draft because your cpp currently only uses
#'   vessel/year/flag + spatial fields + smooths.)
#' @param data A data.frame with required columns.
#' @param DLL Name of compiled TMB DLL (without extension).
#' @param vessel_effect "on" or "off". If "off", vessel RE is fixed at 0 via map.
#' @param q_diffs_spatial "on" or "off". If "off", flag_s_1/2 spatial fields are fixed at 0 via map.
#' @param control List. e.g. control = list(eval.max=1e5, iter.max=1e5, profile=FALSE)
#' @param ... Passed to make_data() (mesh/scale settings).
#'
#' @return A list of class "intCPUE" containing obj, opt, rep, prep, maps, etc.
#' @author Rujia Bi \email{rbi@@iattc.org}
#' @export
intCPUE <- function(
    formula,
    data,
    DLL = "intCPUE",
    vessel_effect = c("on", "off"),
    q_diffs_time = c("on", "off"),
    q_diffs_spatial = c("on", "off"),
    control = list(eval.max = 1e5, iter.max = 1e5),
    ...,
    silent = FALSE
) {
  vessel_effect   <- match.arg(vessel_effect)
  q_diffs_spatial <- match.arg(q_diffs_spatial)
  
  data <- as.data.frame(data)
  
  # ---- 1) Data prep (your function) ----
  prep <- make_data(formula = formula, data_input = data, ...)
  
  # We will work on a copy and enforce 0-based contiguous indices:
  dat <- prep$combined_data
  t0  <- .recode_0based(dat$tid)
  v0  <- .recode_0based(dat$vesid)
  f0  <- .recode_0based(dat$flagid)  # IMPORTANT: reference level should be 0 (cpp assumes fid==0 is baseline)
  
  # Overwrite TMB data list safely:
  data_tmb <- prep$data
  data_tmb$n_t <- t0$n
  data_tmb$n_v <- v0$n
  data_tmb$n_f <- f0$n
  
  data_tmb$t_i <- t0$idx
  data_tmb$v_i <- v0$idx
  data_tmb$f_i <- f0$idx
  
  # also ensure encounter is integer 0/1:
  data_tmb$e_i <- as.integer(data_tmb$e_i)
  
  # ---- 2) Dimension bookkeeping ----
  n_i <- data_tmb$n_i
  n_t <- data_tmb$n_t
  n_v <- data_tmb$n_v
  n_f <- data_tmb$n_f
  n_s <- data_tmb$spde$n_s
  
  # smooth dims from data:
  has_smooths <- isTRUE(as.integer(data_tmb$has_smooths) == 1L)
  K_smooth <- if (has_smooths) ncol(data_tmb$Xs) else 0L
  n_smooth <- if (has_smooths) length(data_tmb$Zs) else 0L
  sum_k <- if (has_smooths && n_smooth > 0L) sum(vapply(data_tmb$Zs, ncol, 0L)) else 0L
  
  # ---- 3) Initial parameters (match your cpp names) ----
  parameters <- .make_parameters_intCPUE(
    n_t = n_t, n_v = n_v, n_f = n_f, n_s = n_s,
    K_smooth = K_smooth, n_smooth = n_smooth, sum_k = sum_k
  )
  
  # ---- 4) Build MAP to turn on/off parts without user touching cpp ----
  map <- .make_map_intCPUE(
    parameters = parameters,
    n_f = n_f,
    vessel_effect = vessel_effect,
    q_diffs_spatial = q_diffs_spatial
  )
  
  # ---- 5) Random effects list ----
  # NOTE:
  # - omega/epsilon/flag_s are GMRF random fields.
  random <- c(
    "omega_s_1", "epsilon_st_1",
    "omega_s_2", "epsilon_st_2"
  )
  if (vessel_effect == "on") {
    random <- c(random, "ves_v_1", "ves_v_2")
  }
  if (q_diffs_spatial == "on" && n_f > 1L) {
    random <- c(random, "flag_s_1", "flag_s_2")
  }
  # Smooth penalized coefficients:
  if (has_smooths && sum_k > 0L) random <- c(random, "b_smooth")
  
  # ---- 6) MakeADFun ----
  .load_intcpue_dll(DLL)
  
  obj <- TMB::MakeADFun(
    data = data_tmb,
    parameters = parameters,
    map = map,
    random = unique(random),
    DLL = DLL,
    silent = silent
  )
  
  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = control
  )
  
  rep <- TMB::sdreport(obj)
  
  out <- list(
    obj = obj,
    opt = opt,
    rep = rep,
    prep = prep,
    data_tmb = data_tmb,
    map = map,
    control = control,
    settings = list(
      vessel_effect = vessel_effect,
      q_diffs_spatial = q_diffs_spatial,
      DLL = DLL
    )
  )
  class(out) <- "intCPUE"
  out
}

# ---- helpers ---------------------------------------------------------------

.recode_0based <- function(x) {
  # robust: works for numeric IDs that are not consecutive
  f <- factor(x, levels = sort(unique(x)))
  idx <- as.integer(f) - 1L
  list(idx = idx, n = nlevels(f), levels = levels(f))
}

.make_parameters_intCPUE <- function(n_t, n_v, n_f, n_s, K_smooth, n_smooth, sum_k) {
  
  # convenience for "empty matrices"
  empty_mat <- function(nr, nc) {
    matrix(0.0, nrow = nr, ncol = nc)
  }
  
  list(
    # obs likelihood
    ln_sd = 0.0,
    
    # anisotropy params: length 2
    ln_H_input = c(0.0, 0.0),
    
    # SPDE hypers (two components)
    ln_range_1 = 0.0,
    ln_sigma_0_1 = 0.0,
    ln_sigma_t_1 = 0.0,
    ln_sigma_flag_1 = 0.0,
    
    ln_range_2 = 0.0,
    ln_sigma_0_2 = 0.0,
    ln_sigma_t_2 = 0.0,
    ln_sigma_flag_2 = 0.0,
    
    # vessel/year/flag fixed effects
    ves_v_1 = rep(0.0, n_v),
    yq_t_1  = rep(0.0, n_t),
    flag_f_1 = if (n_f > 1L) rep(0.0, n_f - 1L) else numeric(0),
    
    ves_ln_std_dev_1 = 0.0,
    flag_ln_std_dev_1 = 0.0,
    
    ves_v_2 = rep(0.0, n_v),
    yq_t_2  = rep(0.0, n_t),
    flag_f_2 = if (n_f > 1L) rep(0.0, n_f - 1L) else numeric(0),
    
    ves_ln_std_dev_2 = 0.0,
    flag_ln_std_dev_2 = 0.0,
    
    # spatial fields
    omega_s_1 = rep(0.0, n_s),
    epsilon_st_1 = empty_mat(n_s, n_t),
    
    omega_s_2 = rep(0.0, n_s),
    epsilon_st_2 = empty_mat(n_s, n_t),
    
    # spatial flag differences (only meaningful if n_f>1)
    flag_s_1 = empty_mat(n_s, max(0L, n_f - 1L)),
    flag_s_2 = empty_mat(n_s, max(0L, n_f - 1L)),
    
    # smoothers: matrices with 2 columns (encounter col0, positive col1)
    bs = empty_mat(K_smooth, 2L),
    b_smooth = empty_mat(sum_k, 2L),
    ln_smooth_sigma = empty_mat(n_smooth, 2L),
    
    # optional epsilon trick
    eps_index = numeric(0)
  )
}

.make_map_intCPUE <- function(parameters, n_f,
                              vessel_effect = c("on","off"),
                              q_diffs_spatial = c("on","off")) {
  
  vessel_effect <- match.arg(vessel_effect)
  q_diffs_spatial <- match.arg(q_diffs_spatial)
  
  map <- list()
  
  # ---- Vessel RE off: fix ves_v_* and ves_ln_std_dev_* ----
  if (vessel_effect == "off") {
    map$ves_v_1 <- factor(rep(NA, length(parameters$ves_v_1)))
    map$ves_v_2 <- factor(rep(NA, length(parameters$ves_v_2)))
    map$ves_ln_std_dev_1 <- factor(NA)
    map$ves_ln_std_dev_2 <- factor(NA)
    
    # also optional: if you want strictly no vessel variance, you can fix the SD too
    # parameters$ves_ln_std_dev_* will remain in par list but fixed by map.
  }
  
  # ---- Spatial flag differences off: fix flag_s_*, ln_sigma_flag_* ----
  if (q_diffs_spatial == "off" || n_f <= 1L) {
    # if n_f==1, these matrices are n_s x 0 anyway; mapping is harmless.
    map$flag_s_1 <- .map_matrix_NA(parameters$flag_s_1)
    map$flag_s_2 <- .map_matrix_NA(parameters$flag_s_2)
    map$ln_sigma_flag_1 <- factor(NA)
    map$ln_sigma_flag_2 <- factor(NA)
  }

  map
}

.load_intcpue_dll <- function(DLL = "intCPUE") {
  if (!DLL %in% names(getLoadedDLLs())) {
    dynlib <- TMB::dynlib(DLL)
    if (!file.exists(dynlib)) {
      stop("Cannot find compiled TMB DLL: ", dynlib, call. = FALSE)
    }
    TMB::dyn.load(dynlib)
  }
  invisible(TRUE)
}

.map_matrix_NA <- function(x) {
  f <- factor(rep(NA, length(x)))
  dim(f) <- dim(x)
  f
}