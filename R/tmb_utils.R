# tmb-utils.R ---------------------------------------------------------------
# Internal helpers for intCPUE TMB workflows.
# - ID recoding / sanity checks for 0-based consecutive indices
# - map helpers for fixing parameters via TMB map
# - small wrappers for optimization / debugging
#
# These are intentionally small, side-effect-free utilities used by fit().

# ---- 1) Index helpers / input checks --------------------------------------
.check_0based_contiguous <- function(x, name, allow_empty = FALSE) {
  if (anyNA(x)) {
    stop(sprintf("`%s` must not contain NA.", name), call. = FALSE)
  }
  
  # allow factor/character but require clean integer after coercion
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  if (any(!is.finite(x_num))) {
    stop(sprintf("`%s` must be numeric/integer (0-based contiguous).", name), call. = FALSE)
  }
  
  # require integer-valued
  if (any(abs(x_num - round(x_num)) > 1e-9)) {
    stop(sprintf("`%s` must be integer-valued (0-based contiguous).", name), call. = FALSE)
  }
  x_int <- as.integer(round(x_num))
  
  if (!allow_empty && length(x_int) == 0L) {
    stop(sprintf("`%s` is empty.", name), call. = FALSE)
  }
  
  if (min(x_int) != 0L) {
    stop(sprintf("`%s` must be 0-based: min(%s) must be 0.", name, name), call. = FALSE)
  }
  
  mx <- max(x_int)
  # must be exactly all integers 0..mx
  tab <- tabulate(x_int + 1L, nbins = mx + 1L)
  if (any(tab == 0L)) {
    miss <- which(tab == 0L) - 1L
    stop(sprintf(
      "`%s` must be contiguous (0..max). Missing values: %s",
      name, paste(miss, collapse = ", ")
    ), call. = FALSE)
  }
  
  invisible(list(x = x_int, n = mx + 1L))
}

.check_fit_inputs_intCPUE <- function(data_tmb) {
  # NOTE:
  # This is the single "gatekeeper" for fit().
  # If it passes, indices and core objects are consistent with what the C++ code expects.
  # If it fails, fix in make_data() (NOT by patching indices manually).
  
  # ---- required fields ----
  req <- c(
    "n_i","n_t","n_v","n_f",
    "b_i","e_i","t_i","v_i","f_i",
    "A_is","A_gs","spde",
    "has_smooths","Xs","Zs"
  )
  miss <- setdiff(req, names(data_tmb))
  if (length(miss)) {
    stop(
      "make_data() returned `data` missing: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  
  # ---- basic type checks ----
  if (!is.numeric(data_tmb$b_i))
    stop("`b_i` must be numeric.", call. = FALSE)
  
  if (anyNA(data_tmb$e_i))
    stop("`e_i` contains NA; encounter must be 0/1 with no NA.", call. = FALSE)
  
  if (!all(data_tmb$e_i %in% c(0L, 1L)))
    stop("`e_i` must contain only 0/1.", call. = FALSE)
  
  # ---- index sanity checks ----
  check_index <- function(x, n, name) {
    if (!is.integer(x)) x <- as.integer(x)
    
    if (length(x) == 0L)
      stop("`", name, "` has length 0; check make_data().", call. = FALSE)
    
    if (anyNA(x))
      stop("`", name, "` contains NA; do not modify indices manually.", call. = FALSE)
    
    # range check (0-based)
    if (min(x) < 0L || max(x) > n - 1L) {
      stop(
        sprintf(
          "`%s` must be 0..%d (0-based). Found range [%d, %d]. Re-run make_data().",
          name, n - 1L, min(x), max(x)
        ),
        call. = FALSE
      )
    }
    
    # consecutive check (CRITICAL):
    # must contain every value 0..n-1 at least once
    ux <- sort(unique(x))
    if (length(ux) != n || !identical(ux, seq.int(0L, n - 1L))) {
      stop(
        sprintf(
          "`%s` must contain consecutive values 0..%d. Found: [%s].\nDid you drop levels or modify indices after make_data()? Re-run make_data().",
          name, n - 1L, paste(utils::head(ux, 12), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    
    invisible(TRUE)
  }
  
  check_index(data_tmb$t_i, data_tmb$n_t, "t_i")
  check_index(data_tmb$v_i, data_tmb$n_v, "v_i")
  check_index(data_tmb$f_i, data_tmb$n_f, "f_i")
  
  invisible(TRUE)
}

# ---- 2) TMB map helpers ----------------------------------------------------

.map_matrix_NA <- function(x) {
  # NOTE:
  # In TMB's `map`, NA means "fixed at its initial value".
  # This helper creates a factor matrix of NAs with the same dimensions as x.
  f <- factor(rep(NA, length(x)))
  dim(f) <- dim(x)
  f
}

.map_matrix_partial_fix <- function(x, keep) {
  # NOTE:
  # Partial mapping for matrix parameters.
  # - keep == TRUE  -> estimate the parameter
  # - keep == FALSE -> fix at initial value (so ensure init value is what you want, e.g. 0)
  #
  # This is useful for sparse designs like flag-by-time effects:
  # you can estimate only cells that are supported by the data and fix the rest to 0.
  
  if (!is.matrix(x) || !is.matrix(keep))
    stop("x and keep must be matrices.", call. = FALSE)
  
  if (!all(dim(x) == dim(keep)))
    stop("x and keep must have same dim.", call. = FALSE)
  
  # Factor indices label the parameters to be estimated.
  # NA entries are fixed at initial value.
  f <- factor(seq_len(length(x)))
  f[!as.vector(keep)] <- NA
  dim(f) <- dim(x)
  f
}

# ---- 3) Optimization helpers ----------------------------------------------

.safe_optimize <- function(obj, control) {
  # NOTE:
  # Small wrapper around nlminb() that returns the result and warns on non-zero convergence.
  # You can expand later to do:
  # - optimizer restarts
  # - Newton steps
  # - gradient checks
  opt <- stats::nlminb(
    start     = obj$par,
    objective = obj$fn,
    gradient  = obj$gr,
    control   = control
  )
  
  if (!is.null(opt$convergence) && opt$convergence != 0) {
    msg <- if (!is.null(opt$message)) opt$message else "no message"
    warning(
      sprintf("Optimizer may not have converged (convergence=%s): %s", opt$convergence, msg),
      call. = FALSE
    )
  }
  
  opt
}

# ---- 4) Debug helpers ------------------------------------------------------

.get_parlist <- function(obj) {
  # NOTE:
  # Convenience function to view parameters in structured list form:
  # useful for checking dimensions, initial values, and whether mapping worked.
  obj$env$parList(obj$par)
}
