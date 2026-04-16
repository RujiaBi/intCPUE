# Internal: build TMB map to switch model components on/off
#
# The map controls which parameters are estimated vs fixed.
# In TMB, map entries set to NA are fixed at their initial value.
#
# Important statistical principle:
# We automatically prevent estimation of variance parameters that are not identifiable
# (e.g., flag temporal random effects with < 2 time points).
#
.build_flag_t_constraint <- function(has_tf) {
  if (is.null(has_tf)) {
    stop("`has_tf` is required to build temporal flag constraints.", call. = FALSE)
  }
  if (!is.matrix(has_tf)) {
    stop("`has_tf` must be a matrix.", call. = FALSE)
  }

  keep_tf <- has_tf > 0
  estimable_flag <- colSums(keep_tf) >= 2L
  flag_t_index <- matrix(-1L, nrow = nrow(keep_tf), ncol = ncol(keep_tf))
  next_idx <- 0L

  if (ncol(keep_tf) > 0L) {
    for (j in seq_len(ncol(keep_tf))) {
      if (!estimable_flag[j]) {
        keep_tf[, j] <- FALSE
        next
      }
      rows_j <- which(keep_tf[, j])
      # Use the first observed time as the reference level and fix it to zero.
      anchor_row <- rows_j[1L]
      keep_tf[anchor_row, j] <- FALSE

      free_rows <- rows_j[-1L]
      if (length(free_rows) > 0L) {
        idx_j <- seq.int(from = next_idx, length.out = length(free_rows))
        flag_t_index[free_rows, j] <- as.integer(idx_j)
        next_idx <- next_idx + length(free_rows)
      }
    }
  }

  list(
    estimable_flag = estimable_flag,
    flag_t_index = flag_t_index,
    n_free = next_idx
  )
}

.make_map_intCPUE <- function(parameters, n_f,
                              obs_sd = c("shared","flag"),
                              pop_spatiotemporal_type = c("rw", "ar1"),
                              q_diffs_time = c("on","off"),
                              has_tf = NULL,
                              estimable_flag = NULL) {
  q_diffs_time <- match.arg(q_diffs_time)
  obs_sd <- match.arg(obs_sd)
  pop_spatiotemporal_type <- match.arg(pop_spatiotemporal_type)
  
  map <- list()
  has_param <- function(name) !is.null(parameters[[name]])

  # ---------------------------------------------------------
  # Observation SD
  # ---------------------------------------------------------
  if (obs_sd == "shared") {
    if (has_param("ln_sd_flag")) {
      map$ln_sd_flag <- factor(rep(NA, length(parameters$ln_sd_flag)))
    }
  } else {
    if (has_param("ln_sd")) {
      map$ln_sd <- factor(NA)
    }
  }

  # Fixed core architecture:
  # spatial + spatiotemporal + vessel + q_diffs_system are always on.
  # q_diffs_time / q_diffs_spatial are handled by separate TMB templates, not map.
  # The AR1 correlation parameters are only active for the AR1 branch; under RW
  # they must be fixed at the neutral working-scale value to avoid singular Hessians.
  if (pop_spatiotemporal_type == "rw") {
    if (has_param("transf_rho_1")) map$transf_rho_1 <- factor(NA)
    if (has_param("transf_rho_2")) map$transf_rho_2 <- factor(NA)
  }

  if (has_param("flag_t_ln_std_dev_1")) {
    if (q_diffs_time == "off" || n_f <= 1L) {
      map$flag_t_ln_std_dev_1 <- factor(NA)
      if (has_param("flag_t_ln_std_dev_2")) map$flag_t_ln_std_dev_2 <- factor(NA)
    } else {
      if (is.null(has_tf)) {
        stop("q_diffs_time='on' requires `has_tf` to construct temporal flag effects.", call. = FALSE)
      }
      if (is.null(estimable_flag)) {
        tf_constraint <- .build_flag_t_constraint(has_tf)
        estimable_flag <- tf_constraint$estimable_flag
      }
      if (!any(estimable_flag)) {
        map$flag_t_ln_std_dev_1 <- factor(NA)
        if (has_param("flag_t_ln_std_dev_2")) map$flag_t_ln_std_dev_2 <- factor(NA)
      }
    }
  }

  # If there is only one flag in the data, the flag-system terms are not identifiable.
  if (n_f <= 1L) {
    if (has_param("flag_f_1")) map$flag_f_1 <- factor(rep(NA, length(parameters$flag_f_1)))
    if (has_param("flag_f_2")) map$flag_f_2 <- factor(rep(NA, length(parameters$flag_f_2)))
    if (has_param("flag_ln_std_dev_1")) map$flag_ln_std_dev_1 <- factor(NA)
    if (has_param("flag_ln_std_dev_2")) map$flag_ln_std_dev_2 <- factor(NA)
  }
  
  map
}
