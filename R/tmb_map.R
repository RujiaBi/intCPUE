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

.make_map_intCPUE <- function(parameters, n_f, n_v,
                              obs_sd = c("shared","flag"),
                              pop_spatiotemporal_type = c("rw", "ar1"),
                              use_random_tid_effect = FALSE,
                              q_diffs_system_type = c("random", "fixed"),
                              q_diffs_time = c("on","off"),
                              has_tf = NULL,
                              estimable_flag = NULL) {
  q_diffs_time <- match.arg(q_diffs_time)
  q_diffs_system_type <- match.arg(q_diffs_system_type)
  obs_sd <- match.arg(obs_sd)
  pop_spatiotemporal_type <- match.arg(pop_spatiotemporal_type)
  
  map <- list()
  has_param <- function(name) !is.null(parameters[[name]]) && length(parameters[[name]]) > 0L
  fix_all <- function(name) factor(rep(NA, length(parameters[[name]])))

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

  # Structural modules such as omega, epsilon, and flag_s are controlled by
  # wrapper templates and omitted from `parameters` when disabled. This map only
  # fixes parameters that remain in the template but are not identifiable for a
  # given data/model configuration.
  # The AR1 correlation parameters are only active for the AR1 branch; under RW
  # they must be fixed at the neutral working-scale value to avoid singular Hessians.
  if (pop_spatiotemporal_type == "rw") {
    if (has_param("transf_rho_1")) map$transf_rho_1 <- fix_all("transf_rho_1")
    if (has_param("transf_rho_2")) map$transf_rho_2 <- fix_all("transf_rho_2")
  }

  if (n_v <= 1L) {
    if (has_param("ves_v_1")) map$ves_v_1 <- factor(rep(NA, length(parameters$ves_v_1)))
    if (has_param("ves_v_2")) map$ves_v_2 <- factor(rep(NA, length(parameters$ves_v_2)))
    if (has_param("ves_ln_std_dev_1")) map$ves_ln_std_dev_1 <- fix_all("ves_ln_std_dev_1")
    if (has_param("ves_ln_std_dev_2")) map$ves_ln_std_dev_2 <- factor(NA)
  }

  if (isTRUE(use_random_tid_effect)) {
    if (has_param("yq_t_1")) map$yq_t_1 <- factor(rep(NA, length(parameters$yq_t_1)))
    if (has_param("yq_t_2")) map$yq_t_2 <- factor(rep(NA, length(parameters$yq_t_2)))
  } else {
    if (has_param("pop_intercept")) map$pop_intercept <- factor(rep(NA, length(parameters$pop_intercept)))
  }

  if (has_param("flag_t_ln_std_dev_1") || has_param("flag_t_ln_std_dev_2")) {
    if (is.null(has_tf)) {
      stop("q_diffs_time='on' requires `has_tf` to construct temporal flag effects.", call. = FALSE)
    }
    if (is.null(estimable_flag)) {
      tf_constraint <- .build_flag_t_constraint(has_tf)
      estimable_flag <- tf_constraint$estimable_flag
    }
    if (!any(estimable_flag)) {
      if (has_param("flag_t_ln_std_dev_1")) map$flag_t_ln_std_dev_1 <- fix_all("flag_t_ln_std_dev_1")
      if (has_param("flag_t_ln_std_dev_2")) map$flag_t_ln_std_dev_2 <- fix_all("flag_t_ln_std_dev_2")
    }
  }

  # If there is only one flag in the data, or the flag main effect is fixed,
  # the flag-system variance terms are not identifiable.
  if (n_f <= 1L || q_diffs_system_type == "fixed") {
    if (has_param("flag_ln_std_dev_1")) map$flag_ln_std_dev_1 <- fix_all("flag_ln_std_dev_1")
    if (has_param("flag_ln_std_dev_2")) map$flag_ln_std_dev_2 <- factor(NA)
  }

  map
}
