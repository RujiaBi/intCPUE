#' Get bias-corrected index and uncertainty from an intCPUE fit
#'
#' Computes a bias-corrected index on the original scale using the "epsilon trick"
#' (via eps_index) and returns log-scale SE from sdreport (ADREPORT(link_total)).
#'
#' @param object An object of class `intCPUE` returned by [intCPUE::intCPUE()].
#' @param level Confidence level for intervals. Default 0.95.
#' @param inner.control List passed to TMB::MakeADFun(inner.control=...).
#'   Default uses sparse + lowrank for memory efficiency.
#' @param silent Logical; passed to TMB::MakeADFun() for the bias-correction step.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item time: 1:n_t
#'     \item index: bias-corrected index (original scale)
#'     \item log_index: log(index)
#'     \item cv: SE on log scale (from sdreport for ADREPORT(link_total))
#'     \item lwr, upr: CI on original scale using lognormal approximation
#'   }
#' @export
get_index <- function(
    object,
    level = 0.95,
    inner.control = list(sparse = TRUE, lowrank = TRUE, trace = FALSE),
    silent = TRUE
) {
  if (!inherits(object, "intCPUE")) {
    stop("`object` must be an `intCPUE` fit from intCPUE().", call. = FALSE)
  }
  
  # --- pull what we need from fit ---
  obj   <- object$obj
  opt   <- object$opt
  rep   <- object$rep
  data  <- object$data_tmb
  DLL   <- object$settings$DLL
  random <- object$random
  ncores <- object$settings$ncores
  
  if (is.null(data$n_t)) {
    stop("`object$data_tmb$n_t` is missing. Cannot determine index length.", call. = FALSE)
  }
  n_t <- as.integer(data$n_t)
  
  # --- 1) log-scale SE from sdreport (ADREPORT(link_total)) ---
  ssr <- summary(rep, "report")
  if (!("link_total" %in% rownames(ssr))) {
    stop("`link_total` not found in sdreport(summary(rep, 'report')). Did you ADREPORT(link_total) in C++?",
         call. = FALSE)
  }
  log_se <- ssr[rownames(ssr) == "link_total", "Std. Error"]
  if (length(log_se) != n_t) {
    stop("Unexpected length of log SE for link_total: got ", length(log_se),
         ", expected n_t=", n_t, ".", call. = FALSE)
  }
  
  # --- 2) bias correction step using eps_index gradient ---
  parhat <- obj$env$parList(opt$par)
  parhat[["eps_index"]] <- rep(0, n_t)
  
  if (!is.null(ncores)) {
    ncores <- as.integer(ncores)
    if (is.na(ncores) || ncores < 1L) {
      stop("`ncores` must be a positive integer.", call. = FALSE)
    }
    if (ncores > 1L) {
      TMB::openmp(ncores, autopar = TRUE, DLL = DLL)
    }
  }
  
  new_obj <- TMB::MakeADFun(
    data = data,
    parameters = parhat,
    map = object$map,
    random = random,
    DLL = DLL,
    silent = silent,
    intern = FALSE,
    inner.control = inner.control
  )
  
  grad <- new_obj$gr()
  nm   <- names(new_obj$par)

  if (is.null(nm) || !any(nm == "eps_index")) {
    stop(
      "Could not find `eps_index` in new_obj$par. Ensure eps_index is a PARAMETER_VECTOR in C++ and is not mapped out.",
      call. = FALSE
    )
  }
  
  index_bc <- as.numeric(grad[nm == "eps_index"])
  if (length(index_bc) != n_t) {
    stop(
      "Bias-corrected index length mismatch: got ", length(index_bc),
      ", expected n_t=", n_t, ".",
      call. = FALSE
    )
  }
  if (any(!is.finite(index_bc)) || any(index_bc <= 0)) {
    stop(
      "Bias-corrected index contains non-finite or non-positive values; cannot take logs. ",
      "Check model fit / bias-correction step.",
      call. = FALSE
    )
  }
  
  # --- 3) assemble output on original scale + log scale ---
  log_index_bc <- log(index_bc)
  z <- stats::qnorm(1 - (1 - level)/2)
  
  out <- data.frame(
    time = seq_len(n_t),
    index = index_bc,
    log_index = log_index_bc,
    cv = log_se,
    lwr = exp(log_index_bc - z * log_se),
    upr = exp(log_index_bc + z * log_se)
  )
  
  out
}

