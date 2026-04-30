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
    "n_i","n_t","n_v","n_f","n_g","n_region",
    "has_component1","all_obs_use_component1","use_split_component1_design","use_poisson_link",
    "use_q_diffs_time",
    "use_flag_f_random",
    "obs_uses_component1","v_component1_keep","v_component1_col",
    "flag_component1_keep",
    "flag_component1_col",
    "flag_t_index_1","flag_t_index_2",
    "b_i","e_i","t_i","v_i","f_i","region_g",
    "A_is","A_gs","spde",
    "has_smooths_catch","Xs_catch","Zs_catch","b_smooth_start_catch",
    "has_fixed_factors_catch","has_random_factors_catch","Xf_catch","Zf_catch","b_factor_start_catch",
    "has_fixed_factors_catch_1","has_random_factors_catch_1","Xf_catch_1","Zf_catch_1","b_factor_start_catch_1",
    "has_smooths_catch_1","Xs_catch_1","Zs_catch_1","b_smooth_start_catch_1",
    "has_smooths_pop","Xs_pop_i","Zs_pop_i","Xs_pop_g","Zs_pop_g","b_smooth_start_pop",
    "has_fixed_factors_pop","has_random_factors_pop","has_pop_intercept","Xf_pop_i","Zf_pop_i","Xf_pop_g","Zf_pop_g","b_factor_start_pop",
    "has_smooths_pop_1","Xs_pop_i_1","Zs_pop_i_1","Xs_pop_g_1","Zs_pop_g_1","b_smooth_start_pop_1",
    "has_fixed_factors_pop_1","has_random_factors_pop_1","Xf_pop_i_1","Zf_pop_i_1","Xf_pop_g_1","Zf_pop_g_1","b_factor_start_pop_1"
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

  if (length(data_tmb$region_g) != data_tmb$n_g) {
    stop("`region_g` must have length `n_g`.", call. = FALSE)
  }
  if (data_tmb$n_region > 0L) {
    check_index(data_tmb$region_g, data_tmb$n_region, "region_g")
  } else if (!all(data_tmb$region_g == -1L)) {
    stop("`region_g` must be all -1 when `n_region = 0`.", call. = FALSE)
  }

  if (length(data_tmb$v_component1_keep) != data_tmb$n_v) {
    stop("`v_component1_keep` must have length `n_v`.", call. = FALSE)
  }
  if (!all(data_tmb$v_component1_keep %in% c(0L, 1L))) {
    stop("`v_component1_keep` must contain only 0/1.", call. = FALSE)
  }
  if (length(data_tmb$v_component1_col) != data_tmb$n_v) {
    stop("`v_component1_col` must have length `n_v`.", call. = FALSE)
  }
  if (!all(data_tmb$v_component1_col >= 0L)) {
    stop("`v_component1_col` must contain non-negative indices.", call. = FALSE)
  }
  n_v_component1 <- sum(data_tmb$v_component1_keep == 1L)
  if (n_v_component1 > 0L &&
      any(data_tmb$v_component1_col[data_tmb$v_component1_keep == 1L] >= n_v_component1)) {
    stop("`v_component1_col` is out of bounds for component-1 vessels.", call. = FALSE)
  }
  if (length(data_tmb$flag_component1_keep) != data_tmb$n_f) {
    stop("`flag_component1_keep` must have length `n_f`.", call. = FALSE)
  }
  if (!all(data_tmb$flag_component1_keep %in% c(0L, 1L))) {
    stop("`flag_component1_keep` must contain only 0/1.", call. = FALSE)
  }
  if (length(data_tmb$flag_component1_col) != data_tmb$n_f) {
    stop("`flag_component1_col` must have length `n_f`.", call. = FALSE)
  }
  if (!all(data_tmb$flag_component1_col >= 0L)) {
    stop("`flag_component1_col` must contain non-negative indices.", call. = FALSE)
  }
  n_f_component1 <- sum(data_tmb$flag_component1_keep == 1L)
  if (n_f_component1 > 0L &&
      any(data_tmb$flag_component1_col[data_tmb$flag_component1_keep == 1L] >= n_f_component1)) {
    stop("`flag_component1_col` is out of bounds for component-1 flags.", call. = FALSE)
  }

  if (!identical(nrow(data_tmb$Xs_catch), data_tmb$n_i))
    stop("`Xs_catch` must have `n_i` rows.", call. = FALSE)

  if (!identical(nrow(data_tmb$Xs_pop_i), data_tmb$n_i))
    stop("`Xs_pop_i` must have `n_i` rows.", call. = FALSE)

  if (!identical(nrow(data_tmb$Xs_pop_g), data_tmb$n_g * data_tmb$n_t))
    stop("`Xs_pop_g` must have `n_g * n_t` rows.", call. = FALSE)

  if (length(data_tmb$Zs_catch) > 0L &&
      any(vapply(data_tmb$Zs_catch, nrow, integer(1)) != data_tmb$n_i)) {
    stop("Each element of `Zs_catch` must have `n_i` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zs_pop_i) > 0L &&
      any(vapply(data_tmb$Zs_pop_i, nrow, integer(1)) != data_tmb$n_i)) {
    stop("Each element of `Zs_pop_i` must have `n_i` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zs_pop_g) > 0L &&
      any(vapply(data_tmb$Zs_pop_g, nrow, integer(1)) != data_tmb$n_g * data_tmb$n_t)) {
    stop("Each element of `Zs_pop_g` must have `n_g * n_t` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zs_pop_i) != length(data_tmb$Zs_pop_g)) {
    stop("`Zs_pop_i` and `Zs_pop_g` must have the same length.", call. = FALSE)
  }

  if (!identical(nrow(data_tmb$Xs_catch_1), data_tmb$n_i))
    stop("`Xs_catch_1` must have `n_i` rows.", call. = FALSE)

  if (!identical(nrow(data_tmb$Xs_pop_i_1), data_tmb$n_i))
    stop("`Xs_pop_i_1` must have `n_i` rows.", call. = FALSE)

  if (!identical(nrow(data_tmb$Xs_pop_g_1), data_tmb$n_g * data_tmb$n_t))
    stop("`Xs_pop_g_1` must have `n_g * n_t` rows.", call. = FALSE)

  if (length(data_tmb$Zs_catch_1) > 0L &&
      any(vapply(data_tmb$Zs_catch_1, nrow, integer(1)) != data_tmb$n_i)) {
    stop("Each element of `Zs_catch_1` must have `n_i` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zs_pop_i_1) > 0L &&
      any(vapply(data_tmb$Zs_pop_i_1, nrow, integer(1)) != data_tmb$n_i)) {
    stop("Each element of `Zs_pop_i_1` must have `n_i` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zs_pop_g_1) > 0L &&
      any(vapply(data_tmb$Zs_pop_g_1, nrow, integer(1)) != data_tmb$n_g * data_tmb$n_t)) {
    stop("Each element of `Zs_pop_g_1` must have `n_g * n_t` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zs_pop_i_1) != length(data_tmb$Zs_pop_g_1)) {
    stop("`Zs_pop_i_1` and `Zs_pop_g_1` must have the same length.", call. = FALSE)
  }

  if (!is.null(data_tmb$Xf_catch) &&
      !identical(nrow(data_tmb$Xf_catch), data_tmb$n_i)) {
    stop("`Xf_catch` must have `n_i` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zf_catch) > 0L &&
      any(vapply(data_tmb$Zf_catch, nrow, integer(1)) != data_tmb$n_i)) {
    stop("Each element of `Zf_catch` must have `n_i` rows.", call. = FALSE)
  }

  if (!is.null(data_tmb$Xf_pop_i) &&
      !identical(nrow(data_tmb$Xf_pop_i), data_tmb$n_i)) {
    stop("`Xf_pop_i` must have `n_i` rows.", call. = FALSE)
  }

  if (!is.null(data_tmb$Xf_pop_g) &&
      !identical(nrow(data_tmb$Xf_pop_g), data_tmb$n_g * data_tmb$n_t)) {
    stop("`Xf_pop_g` must have `n_g * n_t` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zf_pop_i) > 0L &&
      any(vapply(data_tmb$Zf_pop_i, nrow, integer(1)) != data_tmb$n_i)) {
    stop("Each element of `Zf_pop_i` must have `n_i` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zf_pop_g) > 0L &&
      any(vapply(data_tmb$Zf_pop_g, nrow, integer(1)) != data_tmb$n_g * data_tmb$n_t)) {
    stop("Each element of `Zf_pop_g` must have `n_g * n_t` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zf_pop_i) != length(data_tmb$Zf_pop_g)) {
    stop("`Zf_pop_i` and `Zf_pop_g` must have the same length.", call. = FALSE)
  }

  if (!is.null(data_tmb$Xf_catch_1) &&
      !identical(nrow(data_tmb$Xf_catch_1), data_tmb$n_i)) {
    stop("`Xf_catch_1` must have `n_i` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zf_catch_1) > 0L &&
      any(vapply(data_tmb$Zf_catch_1, nrow, integer(1)) != data_tmb$n_i)) {
    stop("Each element of `Zf_catch_1` must have `n_i` rows.", call. = FALSE)
  }

  if (!is.null(data_tmb$Xf_pop_i_1) &&
      !identical(nrow(data_tmb$Xf_pop_i_1), data_tmb$n_i)) {
    stop("`Xf_pop_i_1` must have `n_i` rows.", call. = FALSE)
  }

  if (!is.null(data_tmb$Xf_pop_g_1) &&
      !identical(nrow(data_tmb$Xf_pop_g_1), data_tmb$n_g * data_tmb$n_t)) {
    stop("`Xf_pop_g_1` must have `n_g * n_t` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zf_pop_i_1) > 0L &&
      any(vapply(data_tmb$Zf_pop_i_1, nrow, integer(1)) != data_tmb$n_i)) {
    stop("Each element of `Zf_pop_i_1` must have `n_i` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zf_pop_g_1) > 0L &&
      any(vapply(data_tmb$Zf_pop_g_1, nrow, integer(1)) != data_tmb$n_g * data_tmb$n_t)) {
    stop("Each element of `Zf_pop_g_1` must have `n_g * n_t` rows.", call. = FALSE)
  }

  if (length(data_tmb$Zf_pop_i_1) != length(data_tmb$Zf_pop_g_1)) {
    stop("`Zf_pop_i_1` and `Zf_pop_g_1` must have the same length.", call. = FALSE)
  }
  
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

.validate_nonneg_count <- function(x, name) {
  if (length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single non-negative integer.", name), call. = FALSE)
  }

  x_num <- suppressWarnings(as.numeric(x))
  if (!is.finite(x_num) || abs(x_num - round(x_num)) > 1e-9 || x_num < 0) {
    stop(sprintf("`%s` must be a single non-negative integer.", name), call. = FALSE)
  }

  as.integer(round(x_num))
}

.max_grad_opt <- function(obj, par) {
  if (length(par) == 0L) {
    return(0)
  }
  max(abs(obj$gr(par)))
}

.safe_fn_eval <- function(obj, par) {
  val <- try(obj$fn(par), silent = TRUE)
  if (inherits(val, "try-error") || !is.finite(val)) {
    return(Inf)
  }
  val
}

.safe_gr_eval <- function(obj, par) {
  gr <- try(obj$gr(par), silent = TRUE)
  if (inherits(gr, "try-error") || any(!is.finite(gr))) {
    return(NULL)
  }
  gr
}

.run_nlminb_once <- function(obj, start, control) {
  stats::nlminb(
    start     = start,
    objective = obj$fn,
    gradient  = obj$gr,
    control   = control
  )
}

.should_restart_optimize <- function(opt, grad, grad_tol) {
  bad_conv <- !is.null(opt$convergence) && opt$convergence != 0
  bad_grad <- is.finite(grad) && grad > grad_tol
  bad_conv || bad_grad
}

.try_newton_refine <- function(obj, start, grad_tol, trace = 0L) {
  if (!is.function(obj$fn) || !is.function(obj$gr) || !is.function(obj$he)) {
    return(list(ok = FALSE, result = NULL, error = "Missing fn/gr/he on TMB object."))
  }
  if (length(start) == 0L) {
    return(list(ok = FALSE, result = NULL, error = "Empty parameter vector."))
  }

  he_sparse <- function(par) {
    Matrix::Matrix(obj$he(par), sparse = TRUE)
  }

  out <- tryCatch(
    TMB::newton(
      par = start,
      fn = obj$fn,
      gr = obj$gr,
      he = he_sparse,
      trace = as.integer(trace),
      maxit = 50,
      tol = grad_tol,
      grad.tol = grad_tol,
      step.tol = 1e-10,
      tol10 = 1e-6,
      smartsearch = TRUE
    ),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    return(list(ok = FALSE, result = NULL, error = conditionMessage(out)))
  }

  list(ok = TRUE, result = out, error = NULL)
}

.try_coordinate_refine <- function(obj, start, grad_tol, maxit = 5L, top_k = 5L) {
  if (length(start) == 0L) {
    return(NULL)
  }

  par <- start
  best_value <- .safe_fn_eval(obj, par)
  grad_vec <- .safe_gr_eval(obj, par)
  if (!is.finite(best_value) || is.null(grad_vec)) {
    return(NULL)
  }

  best_grad <- max(abs(grad_vec))
  used <- 0L

  for (iter in seq_len(maxit)) {
    if (!is.finite(best_grad) || best_grad <= grad_tol) {
      break
    }

    ord <- order(abs(grad_vec), decreasing = TRUE)
    ord <- ord[seq_len(min(length(ord), max(1L, as.integer(top_k))))]
    improved <- FALSE

    for (idx in ord) {
      base <- par[idx]
      g_i <- grad_vec[idx]
      step0 <- max(
        1e-3,
        0.1 * max(1, abs(base)),
        0.05 * abs(g_i)
      )
      deltas <- step0 * c(-4, -2, -1, -0.5, 0.5, 1, 2, 4)

      cand_vals <- rep(Inf, length(deltas))
      cand_pars <- vector("list", length(deltas))
      for (k in seq_along(deltas)) {
        par_k <- par
        par_k[idx] <- base + deltas[k]
        cand_pars[[k]] <- par_k
        cand_vals[k] <- .safe_fn_eval(obj, par_k)
      }

      k_best <- which.min(cand_vals)
      if (length(k_best) == 0L || !is.finite(cand_vals[k_best])) {
        next
      }

      if (cand_vals[k_best] < best_value - 1e-8) {
        par <- cand_pars[[k_best]]
        best_value <- cand_vals[k_best]
        grad_vec_new <- .safe_gr_eval(obj, par)
        if (is.null(grad_vec_new)) {
          return(NULL)
        }
        grad_vec <- grad_vec_new
        best_grad <- max(abs(grad_vec))
        used <- used + 1L
        improved <- TRUE
        break
      }
    }

    if (!improved) {
      break
    }
  }

  list(
    par = par,
    value = best_value,
    gradient = grad_vec,
    max_grad = best_grad,
    used = used
  )
}

.safe_optimize <- function(obj, control, restart_max = 1L, grad_tol = 1e-3, newton_max = 2L, coord_max = 5L) {
  # NOTE:
  # Small wrapper around nlminb() with a couple of automatic restarts.
  # This helps when nlminb reports false convergence while gradients remain non-trivial.
  restart_max <- .validate_nonneg_count(restart_max, "restart_max")
  newton_max <- .validate_nonneg_count(newton_max, "newton_max")
  coord_max <- .validate_nonneg_count(coord_max, "coord_max")

  opt <- .run_nlminb_once(obj = obj, start = obj$par, control = control)
  best_opt <- opt
  best_grad <- .max_grad_opt(obj, best_opt$par)
  newton_used <- 0L
  coord_used <- 0L
  newton_attempted <- 0L
  coord_attempted <- 0L
  newton_error <- NULL
  coord_error <- NULL

  if (restart_max > 0L && .should_restart_optimize(best_opt, best_grad, grad_tol)) {
    for (iter in seq_len(restart_max)) {
      restart_opt <- .run_nlminb_once(
        obj = obj,
        start = best_opt$par,
        control = control
      )
      restart_grad <- .max_grad_opt(obj, restart_opt$par)

      better_value <- is.finite(restart_opt$objective) &&
        (!is.finite(best_opt$objective) || restart_opt$objective < best_opt$objective - 1e-8)
      better_grad <- is.finite(restart_grad) &&
        (!is.finite(best_grad) || restart_grad < best_grad - 1e-8)
      better_conv <- !is.null(best_opt$convergence) &&
        !is.null(restart_opt$convergence) &&
        restart_opt$convergence < best_opt$convergence

      if (better_value || better_grad || better_conv) {
        best_opt <- restart_opt
        best_grad <- restart_grad
      } else {
        break
      }

      if (!.should_restart_optimize(best_opt, best_grad, grad_tol)) {
        break
      }
    }
  }

  if (.should_restart_optimize(best_opt, best_grad, grad_tol) && newton_max > 0L) {
    for (iter in seq_len(newton_max)) {
      newton_attempted <- newton_attempted + 1L
      newton_opt <- .try_newton_refine(
        obj = obj,
        start = best_opt$par,
        grad_tol = grad_tol,
        trace = 0L
      )

      if (!isTRUE(newton_opt$ok)) {
        newton_error <- newton_opt$error
        break
      }
      newton_res <- newton_opt$result

      newton_grad <- max(abs(newton_res$gradient))
      better_value <- is.finite(newton_res$value) &&
        (!is.finite(best_opt$objective) || newton_res$value < best_opt$objective - 1e-8)
      better_grad <- is.finite(newton_grad) &&
        (!is.finite(best_grad) || newton_grad < best_grad - 1e-8)

      if (!(better_value || better_grad)) {
        break
      }

      best_opt$par <- newton_res$par
      best_opt$objective <- newton_res$value
      best_opt$max_grad <- newton_grad
      best_opt$convergence <- if (newton_grad <= grad_tol) 0L else best_opt$convergence
      best_grad <- newton_grad
      newton_used <- iter

      if (!.should_restart_optimize(best_opt, best_grad, grad_tol)) {
        break
      }
    }

    if (newton_used > 0L) {
      msg <- if (!is.null(best_opt$message)) best_opt$message else "no message"
      best_opt$message <- paste0(msg, "; refined with TMB::newton x", newton_used)
    }
  }

  if (.should_restart_optimize(best_opt, best_grad, grad_tol) && coord_max > 0L) {
    coord_attempted <- 1L
    coord_opt <- .try_coordinate_refine(
      obj = obj,
      start = best_opt$par,
      grad_tol = grad_tol,
      maxit = coord_max,
      top_k = 5L
    )

    if (!is.null(coord_opt) && coord_opt$used > 0L) {
      coord_grad <- coord_opt$max_grad
      better_value <- is.finite(coord_opt$value) &&
        (!is.finite(best_opt$objective) || coord_opt$value < best_opt$objective - 1e-8)
      better_grad <- is.finite(coord_grad) &&
        (!is.finite(best_grad) || coord_grad < best_grad - 1e-8)

      if (better_value || better_grad) {
        best_opt$par <- coord_opt$par
        best_opt$objective <- coord_opt$value
        best_opt$max_grad <- coord_grad
        best_opt$convergence <- if (coord_grad <= grad_tol) 0L else best_opt$convergence
        best_grad <- coord_grad
        coord_used <- coord_opt$used
        msg <- if (!is.null(best_opt$message)) best_opt$message else "no message"
        best_opt$message <- paste0(msg, "; refined with coordinate-polish x", coord_used)
      }
    } else if (is.null(coord_opt)) {
      coord_error <- "Coordinate polish failed to obtain a usable objective/gradient path."
    }
  }

  best_opt$max_grad <- best_grad
  best_opt$newton_used <- newton_used
  best_opt$coord_used <- coord_used
  best_opt$newton_attempted <- newton_attempted
  best_opt$coord_attempted <- coord_attempted
  best_opt$newton_error <- newton_error
  best_opt$coord_error <- coord_error

  if ((!is.null(best_opt$convergence) && best_opt$convergence != 0) || best_grad > grad_tol) {
    msg <- if (!is.null(best_opt$message)) best_opt$message else "no message"
    if (!is.null(best_opt$convergence) && best_opt$convergence == 0 && best_grad > grad_tol) {
      warning(
        sprintf(
          paste0(
            "Optimizer stopped with convergence=0, but max_grad=%.6g exceeds ",
            "the intCPUE tolerance %.6g; treating this as not fully converged. ",
            "nlminb message: %s"
          ),
          best_grad,
          grad_tol,
          msg
        ),
        call. = FALSE
      )
    } else {
      warning(
        sprintf(
          "Optimizer may not have converged (convergence=%s, max_grad=%.6g): %s",
          best_opt$convergence,
          best_grad,
          msg
        ),
        call. = FALSE
      )
    }
  }

  best_opt
}

.safe_sdreport <- function(obj, opt) {
  rep_default <- try(TMB::sdreport(obj), silent = TRUE)
  if (!inherits(rep_default, "try-error")) {
    return(list(
      rep = rep_default,
      method = "default",
      error = NULL
    ))
  }

  default_msg <- conditionMessage(attr(rep_default, "condition"))

  list(
    rep = NULL,
    method = "failed",
    error = paste0("Default sdreport failed: ", default_msg)
  )
}

# ---- 4) Debug helpers ------------------------------------------------------

.get_parlist <- function(obj) {
  # NOTE:
  # Convenience function to view parameters in structured list form:
  # useful for checking dimensions, initial values, and whether mapping worked.
  obj$env$parList(obj$par)
}
