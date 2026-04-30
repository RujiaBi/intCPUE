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
#' @param regions `"total"` returns only the total-region index, preserving the
#'   historical output shape. `"all"` returns region-specific indices plus the
#'   total index. Bias correction is done jointly either way.
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
    silent = TRUE,
    regions = c("total", "all")
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
  regions <- match.arg(regions)
  n_t <- as.integer(data$n_t)
  n_region <- if (!is.null(data$n_region)) as.integer(data$n_region) else 0L
  n_index <- n_region + 1L
  region_labels <- object$prep$region_labels
  if (is.null(region_labels) || length(region_labels) != n_region) {
    region_labels <- as.character(seq_len(n_region) - 1L)
  }
  index_labels <- c(region_labels, "total")
  
  # --- 1) log-scale SE from sdreport (ADREPORT(link_index)) ---
  ssr <- summary(rep, "report")
  if ("link_index" %in% rownames(ssr)) {
    log_se <- ssr[rownames(ssr) == "link_index", "Std. Error"]
    if (length(log_se) != n_index * n_t) {
      stop("Unexpected length of log SE for link_index: got ", length(log_se),
           ", expected ", n_index * n_t, ".", call. = FALSE)
    }
    log_se_mat <- matrix(log_se, nrow = n_index, ncol = n_t)
  } else if ("link_total" %in% rownames(ssr)) {
    log_se_total <- ssr[rownames(ssr) == "link_total", "Std. Error"]
    if (length(log_se_total) != n_t) {
      stop("Unexpected length of log SE for link_total: got ", length(log_se_total),
           ", expected n_t=", n_t, ".", call. = FALSE)
    }
    log_se_mat <- matrix(NA_real_, nrow = n_index, ncol = n_t)
    log_se_mat[n_index, ] <- log_se_total
  } else {
    stop("Neither `link_index` nor `link_total` was found in sdreport(summary(rep, 'report')).",
         call. = FALSE)
  }
  
  # --- 2) joint bias correction step using eps_index gradient ---
  parhat <- obj$env$parList(opt$par)
  parhat[["eps_index"]] <- rep(0, n_index * n_t)
  
  if (!is.null(ncores)) {
    ncores <- as.integer(ncores)
    if (is.na(ncores) || ncores < 1L) {
      stop("`ncores` must be a positive integer.", call. = FALSE)
    }
    TMB::openmp(ncores, autopar = ncores > 1L, DLL = DLL)
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
  if (length(index_bc) != n_index * n_t) {
    stop(
      "Bias-corrected index length mismatch: got ", length(index_bc),
      ", expected ", n_index * n_t, ".",
      call. = FALSE
    )
  }
  index_bc_mat <- matrix(index_bc, nrow = n_index, ncol = n_t)
  if (any(!is.finite(index_bc_mat)) || any(index_bc_mat <= 0)) {
    stop(
      "Bias-corrected index contains non-finite or non-positive values; cannot take logs. ",
      "Check model fit / bias-correction step.",
      call. = FALSE
    )
  }
  
  # --- 3) assemble output on original scale + log scale ---
  keep <- if (identical(regions, "total")) n_index else seq_len(n_index)
  index_bc_vec <- as.vector(index_bc_mat[keep, , drop = FALSE])
  log_index_bc <- log(index_bc_vec)
  log_se_vec <- as.vector(log_se_mat[keep, , drop = FALSE])
  z <- stats::qnorm(1 - (1 - level)/2)
  
  out <- data.frame(
    region = rep(index_labels[keep], times = n_t),
    time = rep(seq_len(n_t), each = length(keep)),
    index = index_bc_vec,
    log_index = log_index_bc,
    cv = log_se_vec,
    lwr = exp(log_index_bc - z * log_se_vec),
    upr = exp(log_index_bc + z * log_se_vec)
  )
  if (identical(regions, "total")) {
    out$region <- NULL
  }
  
  out
}

