.normalize_smoother_output <- function(sm, n_rows) {
  if (!isTRUE(sm$has_smooths)) {
    return(list(
      has_smooths = FALSE,
      Xs = matrix(0, nrow = n_rows, ncol = 0L),
      Zs = list(),
      basis_out = list(),
      labels = list(),
      classes = list(),
      sm_dims = integer(0),
      b_smooth_start = integer(0)
    ))
  }

  Xs <- sm$Xs
  Zs <- sm$Zs
  sm_dims <- as.integer(sm$sm_dims)
  b_smooth_start <- as.integer(sm$b_smooth_start)

  if (!identical(nrow(Xs), n_rows)) {
    stop("parse_smoothers() returned Xs with nrow != expected nrow.", call. = FALSE)
  }

  if (length(Zs)) {
    for (k in seq_along(Zs)) {
      if (!identical(nrow(Zs[[k]]), n_rows)) {
        stop("parse_smoothers() returned Zs[[k]] with nrow != expected nrow.", call. = FALSE)
      }
      if (!inherits(Zs[[k]], "sparseMatrix")) {
        Zs[[k]] <- Matrix::Matrix(Zs[[k]], sparse = TRUE)
      }
    }
  }

  list(
    has_smooths = TRUE,
    Xs = Xs,
    Zs = Zs,
    basis_out = sm$basis_out,
    labels = sm$labels,
    classes = sm$classes,
    sm_dims = sm_dims,
    b_smooth_start = b_smooth_start
  )
}

.mask_formula_vars <- function(data, formula, keep_rows) {
  if (is.null(formula)) {
    return(data)
  }
  vars <- intersect(all.vars(formula), names(data))
  if (!length(vars)) {
    return(data)
  }
  out <- data
  out[!keep_rows, vars] <- NA
  out
}

.mask_unseen_factor_levels <- function(newdata, basis_out) {
  if (!length(basis_out)) {
    return(newdata)
  }
  out <- newdata
  for (basis_i in basis_out) {
    var_i <- basis_i$var
    if (!is.null(var_i) && var_i %in% names(out) && !is.null(basis_i$levels)) {
      vals <- as.character(out[[var_i]])
      vals[!is.na(vals) & !vals %in% basis_i$levels] <- NA
      out[[var_i]] <- vals
    }
  }
  out
}

.component1_supported_formula <- function(formula, data, keep_rows, min_smooth_rows = 5L) {
  if (is.null(formula)) {
    return(list(formula = NULL, dropped = character(0)))
  }
  terms <- all_terms(formula)
  if (!length(terms)) {
    return(list(formula = NULL, dropped = character(0)))
  }

  keep_term <- logical(length(terms))
  for (i in seq_along(terms)) {
    term_i <- terms[i]
    if (grepl("s\\(", term_i)) {
      expr <- str2expression(term_i)[[1]]
      eval_env <- new.env(parent = baseenv())
      eval_env$s <- mgcv::s
      sm <- eval(expr, envir = eval_env)
      keep_i <- keep_rows & .complete_rows_for_sm(sm, data)
      keep_term[i] <- sum(keep_i) >= min_smooth_rows
    } else if (grepl("^f\\(", term_i)) {
      var_i <- .parse_factor_var(term_i, random = TRUE)
      vals_i <- if (var_i %in% names(data)) data[[var_i]][keep_rows & !is.na(data[[var_i]])] else NULL
      keep_term[i] <- length(unique(as.character(vals_i))) >= 2L
    } else {
      var_i <- .parse_factor_var(term_i, random = FALSE)
      vals_i <- if (var_i %in% names(data)) data[[var_i]][keep_rows & !is.na(data[[var_i]])] else NULL
      keep_term[i] <- length(unique(as.character(vals_i))) >= 2L
    }
  }

  kept <- terms[keep_term]
  list(
    formula = if (length(kept)) stats::reformulate(kept) else NULL,
    dropped = terms[!keep_term]
  )
}

.normalize_factor_output <- function(fac, n_rows) {
  Xf <- fac$Xf
  Zf <- fac$Zf
  rand_dims <- as.integer(fac$rand_dims)
  rand_start <- as.integer(fac$rand_start)

  if (!identical(nrow(Xf), n_rows)) {
    stop("parse_factors() returned Xf with nrow != expected nrow.", call. = FALSE)
  }

  if (length(Zf)) {
    for (k in seq_along(Zf)) {
      if (!identical(nrow(Zf[[k]]), n_rows)) {
        stop("parse_factors() returned Zf[[k]] with nrow != expected nrow.", call. = FALSE)
      }
      if (!inherits(Zf[[k]], "sparseMatrix")) {
        Zf[[k]] <- Matrix::Matrix(Zf[[k]], sparse = TRUE)
      }
    }
  }

  list(
    Xf = Xf,
    Zf = Zf,
    has_factor_fixed = isTRUE(fac$has_factor_fixed),
    has_factor_random = isTRUE(fac$has_factor_random),
    basis_out = fac$basis_out,
    rand_dims = rand_dims,
    rand_start = rand_start
  )
}

.smooth_vars_from_basis <- function(basis_out) {
  .keep_var_names <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x) & x != "NA"]
    x
  }

  needed <- character(0)

  if (!length(basis_out)) {
    return(needed)
  }

  for (basis_i in basis_out) {
    if (!length(basis_i)) next
    for (sm in basis_i) {
      term_i <- sm$term
      if (length(term_i)) {
        term_i <- .keep_var_names(term_i)
        needed <- c(needed, term_i)
      }
      if (!is.null(sm$by) && is.character(sm$by) && length(sm$by) == 1L &&
          length(.keep_var_names(sm$by)) == 1L) {
        needed <- c(needed, sm$by)
      }
    }
  }

  unique(.keep_var_names(needed))
}

.expand_projection_over_time <- function(projection_data, tid_values) {
  projection_data <- as.data.frame(projection_data)
  tid_values <- as.integer(tid_values)

  out <- projection_data[rep(seq_len(nrow(projection_data)), times = length(tid_values)), , drop = FALSE]
  out$tid <- rep(tid_values, each = nrow(projection_data))
  rownames(out) <- NULL
  out
}

.warn_projection_na <- function(projection_data, needed_vars) {
  if (!length(needed_vars)) {
    return(invisible(NULL))
  }

  na_counts <- vapply(
    needed_vars,
    function(nm) sum(is.na(projection_data[[nm]])),
    integer(1)
  )
  na_counts <- na_counts[na_counts > 0L]

  if (!length(na_counts)) {
    return(invisible(NULL))
  }

  warning(
    paste0(
      "`projection_data` contains NA in population covariates: ",
      paste(sprintf("%s (%d)", names(na_counts), as.integer(na_counts)), collapse = ", "),
      ". For projection rows with NA, the affected population smooth term(s) will contribute 0."
    ),
    call. = FALSE
  )

  invisible(NULL)
}

.align_projection_data_to_key <- function(projection_data, key, needed_vars, tid_values) {
  projection_data <- as.data.frame(projection_data)
  tid_values <- as.integer(tid_values)

  coord_cols <- c("utm_x_scale", "utm_y_scale")
  has_tid <- "tid" %in% names(projection_data)
  req <- unique(c(coord_cols, if (has_tid) "tid", needed_vars))
  miss <- setdiff(req, names(projection_data))
  if (length(miss)) {
    stop(
      "`projection_data` is missing required columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  if (has_tid) {
    if (anyNA(projection_data$tid)) {
      stop("`projection_data$tid` must not contain NA.", call. = FALSE)
    }

    tid_num <- suppressWarnings(as.numeric(as.character(projection_data$tid)))
    if (any(!is.finite(tid_num)) || any(abs(tid_num - round(tid_num)) > 1e-9)) {
      stop("`projection_data$tid` must be integer-valued and use the same coding as `data_utm$tid`.", call. = FALSE)
    }
    projection_data$tid <- as.integer(round(tid_num))

    bad_tid <- setdiff(unique(projection_data$tid), tid_values)
    if (length(bad_tid)) {
      stop(
        "`projection_data$tid` contains values not used by the model: ",
        paste(sort(bad_tid), collapse = ", "),
        call. = FALSE
      )
    }

    proj_id <- paste(projection_data$utm_x_scale, projection_data$utm_y_scale, projection_data$tid, sep = "\r")
    key_gt <- .expand_projection_over_time(key[, coord_cols, drop = FALSE], tid_values)
    key_id <- paste(key_gt$utm_x_scale, key_gt$utm_y_scale, key_gt$tid, sep = "\r")

    if (anyDuplicated(proj_id)) {
      stop("`projection_data` must contain at most one row per extrapolation cell-time combination.", call. = FALSE)
    }

    idx <- match(key_id, proj_id)
    if (anyNA(idx)) {
      stop(
        "`projection_data` must contain one row for every extrapolation grid cell-time combination ",
        "used by the model. Missing combinations detected.",
        call. = FALSE
      )
    }

    out <- projection_data[idx, req, drop = FALSE]
  } else {
    proj_id <- paste(projection_data$utm_x_scale, projection_data$utm_y_scale, sep = "\r")
    key_id <- paste(key$utm_x_scale, key$utm_y_scale, sep = "\r")

    if (anyDuplicated(proj_id)) {
      stop("`projection_data` must contain at most one row per extrapolation cell.", call. = FALSE)
    }

    idx <- match(key_id, proj_id)
    if (anyNA(idx)) {
      stop(
        "`projection_data` must contain one row for every extrapolation grid cell ",
        "used by the model. Missing cells detected.",
        call. = FALSE
      )
    }

    out <- .expand_projection_over_time(
      projection_data[idx, unique(c(coord_cols, needed_vars)), drop = FALSE],
      tid_values = tid_values
    )
  }

  rownames(out) <- NULL
  out
}

.default_population_projection_data <- function(data_utm, key, basis_out, tid_values) {
  needed_vars <- .smooth_vars_from_basis(basis_out)
  coord_cols <- c("utm_x_scale", "utm_y_scale")

  if (!length(needed_vars)) {
    out <- .expand_projection_over_time(key[, coord_cols, drop = FALSE], tid_values = tid_values)
    rownames(out) <- NULL
    return(out)
  }

  miss <- setdiff(needed_vars, names(data_utm))
  if (length(miss)) {
    stop(
      "Population smooth terms require covariates not found in `data_utm`: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  cell_id <- paste(data_utm$utm_x_scale, data_utm$utm_y_scale, sep = "\r")
  key_id <- paste(key$utm_x_scale, key$utm_y_scale, sep = "\r")

  out_g <- key[, coord_cols, drop = FALSE]
  for (nm in needed_vars) {
    proto <- data_utm[[nm]]
    out_g[[nm]] <- proto[rep(NA_integer_, nrow(key))]
  }

  for (g in seq_len(nrow(key))) {
    rows <- which(cell_id == key_id[g])
    if (!length(rows)) {
      stop("Internal error: extrapolation cell not found in `data_utm`.", call. = FALSE)
    }

    for (nm in needed_vars) {
      vals <- data_utm[[nm]][rows]
      vals_non_na <- vals[!is.na(vals)]
      uniq <- unique(vals_non_na)
      if (length(uniq) > 1L) {
        stop(
          "Population smooth covariate `", nm, "` varies within an extrapolation grid cell. ",
          "Please supply `projection_data` explicitly. Use a `tid` column there ",
          "if the covariate varies over time.",
          call. = FALSE
        )
      }
      if (length(uniq) == 1L) {
        out_g[[nm]][g] <- uniq
      }
    }
  }

  out <- .expand_projection_over_time(out_g, tid_values = tid_values)
  rownames(out) <- NULL
  out
}

.derive_region_g <- function(data_utm, key) {
  coord_cols <- c("utm_x_scale", "utm_y_scale")

  if (!"region" %in% names(data_utm)) {
    return(list(region_g = rep(-1L, nrow(key)), region_labels = character(0)))
  }

  source <- as.data.frame(data_utm)
  miss <- setdiff(c(coord_cols, "region"), names(source))
  if (length(miss)) {
    stop("`data_utm` is missing required region columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (anyNA(source$region)) {
    stop("`region` must not contain NA.", call. = FALSE)
  }

  cell_id <- paste(source$utm_x_scale, source$utm_y_scale, sep = "\r")
  key_id <- paste(key$utm_x_scale, key$utm_y_scale, sep = "\r")
  region_chr <- as.character(source$region)

  region_by_key <- character(nrow(key))
  for (g in seq_len(nrow(key))) {
    vals <- unique(region_chr[cell_id == key_id[g]])
    vals <- vals[nzchar(vals)]
    if (!length(vals)) {
      stop("`region` must be supplied for every extrapolation grid cell.", call. = FALSE)
    }
    if (length(vals) > 1L) {
      stop("Each extrapolation grid cell must belong to exactly one `region`.", call. = FALSE)
    }
    region_by_key[g] <- vals
  }

  labels <- sort(unique(region_by_key))
  region_g <- match(region_by_key, labels) - 1L
  list(region_g = as.integer(region_g), region_labels = labels)
}

#' Prepare data objects and mesh for intCPUE workflows
#'
#' Main data-prep function for intCPUE.
#' - mesh/SPDE/A matrices
#' - extrapolation key grid + areas
#' - parse catchability and population smooths (mgcv `s()`) into Xs/Zs
#'
#' @param data_utm A data.frame containing required columns (with utm_x/y_scale)
#' @param mesh intCPUEmesh built from make_mesh(), or a custom mesh
#' @param formula_catchability Optional one-sided formula defining
#'   smooth and factor terms that affect catchability only.
#' @param formula_population Optional one-sided formula defining
#'   smooth and factor terms that affect both the observation model and the
#'   projected population density surface.
#' @param projection_data Optional data.frame giving projection-grid covariates
#'   for `formula_population`. It must contain `utm_x_scale` and `utm_y_scale`
#'   matching the extrapolation grid. If it also contains `tid`, it is treated
#'   as time-varying projection data with one row per grid cell-time
#'   combination. If omitted, `make_data()` attempts to derive a static
#'   projection table from `data_utm` and replicate it across time; this only
#'   works when each population-smooth covariate is uniquely defined within
#'   each grid cell. In summary:
#'   (1) static population covariates: one row per grid cell;
#'   (2) time-varying population covariates: one row per grid cell-time
#'   combination, with `tid` coded the same way as `data_utm$tid`;
#'   (3) if both static and time-varying covariates are used together, provide
#'   one row per grid cell-time combination and repeat the static covariate
#'   values across `tid`.
#' @param flag_uses_component1 Optional logical vector with one value per
#'   0-based flag level, indicating which flags contribute to the encounter
#'   component. Used internally for mixed delta/positive data. If `NULL`, all
#'   flags are treated as component-1 flags.
#' @param area_scale Numeric or "auto". Scaling factor for area_km2.
#'
#' @return A list with elements mesh, data, key, scales, smooth_basis, smooth_info.
#' @author Rujia Bi \email{rbi@@iattc.org}
#' @export
make_data <- function(
    data_utm,
    mesh,
    formula_catchability = NULL,
    formula_population = NULL,
    projection_data = NULL,
    flag_uses_component1 = NULL,
    area_scale = "auto"
) {
  data_utm <- as.data.frame(data_utm)

  .check_required_cols(data_utm, c("cpue", "encounter", "lon", "lat", "vesid", "tid", "flagid", "utm_x_scale", "utm_y_scale"))
  .check_numeric(data_utm, c("cpue", "encounter", "lon", "lat", "vesid", "tid", "flagid", "utm_x_scale", "utm_y_scale"))

  if (anyNA(data_utm$lon) || anyNA(data_utm$lat)) {
    stop("`lon`/`lat` must not contain NA.", call. = FALSE)
  }

  if (anyNA(data_utm$utm_x_scale) || anyNA(data_utm$utm_y_scale)) {
    stop("`utm_x_scale`/`utm_y_scale` must not contain NA.", call. = FALSE)
  }

  # ---- SPDE + A matrix (handle intCPUEmesh or bare mesh) ----
  loc_xy <- as.matrix(data_utm[, c("utm_x_scale", "utm_y_scale"), drop = FALSE])

  mesh_in <- mesh
  mesh_obj <- .as_intCPUEmesh(
    mesh = mesh_in,
    loc_xy = loc_xy,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    recompute_A = "auto"
  )

  if (!is.null(mesh_obj$loc_xy)) {
    r1 <- range(mesh_obj$loc_xy[, 1])
    r2 <- range(loc_xy[, 1])
    if (is.finite(r1[1]) && is.finite(r2[1])) {
      if (abs(diff(r1) - diff(r2)) / max(1e-12, diff(r2)) > 0.5) {
        warning("Mesh coordinate scale may not match `utm_x_scale/utm_y_scale`. Check scaling.", call. = FALSE)
      }
    }
  }

  mesh <- mesh_obj$mesh
  spde <- mesh_obj$spde

  # ---- A matrices ----
  A_is <- mesh_obj$A
  A_isT <- methods::as(A_is, "TsparseMatrix")
  Ais_ij <- cbind(A_isT@i, A_isT@j)
  Ais_x <- A_isT@x

  # ---- key/extrapolation grid ----
  key_out <- .prep_key_area(data_utm, mesh, area_scale = area_scale)
  key <- key_out$key
  A_gs <- key_out$A_gs
  region_info <- .derive_region_g(data_utm = data_utm, key = key)

  n_i <- nrow(data_utm)
  n_g <- nrow(key)

  # ---- user-supplied 0-based indices: validate only ----
  t_chk <- .check_0based_contiguous(data_utm$tid, "tid")
  v_chk <- .check_0based_contiguous(data_utm$vesid, "vesid")
  f_chk <- .check_0based_contiguous(data_utm$flagid, "flagid")

  t_i <- t_chk$x
  v_i <- v_chk$x
  f_i <- f_chk$x

  n_t <- t_chk$n
  n_v <- v_chk$n
  n_f <- f_chk$n
  tid_values <- seq.int(0L, n_t - 1L)

  if (is.null(flag_uses_component1)) {
    flag_uses_component1 <- rep(TRUE, n_f)
  } else {
    flag_uses_component1 <- as.logical(flag_uses_component1)
    if (!identical(length(flag_uses_component1), n_f)) {
      stop("`flag_uses_component1` must have one value per recoded flag level.", call. = FALSE)
    }
    if (anyNA(flag_uses_component1)) {
      stop("`flag_uses_component1` must not contain NA.", call. = FALSE)
    }
  }

  has_component1 <- any(flag_uses_component1)
  obs_uses_component1 <- flag_uses_component1[f_i + 1L]
  obs_component1_rows <- if (has_component1) obs_uses_component1 else rep(FALSE, n_i)
  split_component1_design <- has_component1 && !all(obs_uses_component1)

  vessel_uses_component1 <- tabulate(v_i[obs_uses_component1] + 1L, nbins = n_v) > 0L
  v_component1_keep <- as.integer(vessel_uses_component1)
  if (sum(v_component1_keep) <= 1L) v_component1_keep[] <- 0L
  v_component1_col <- integer(n_v)
  component1_vessels <- which(v_component1_keep == 1L)
  if (length(component1_vessels) > 0L) {
    v_component1_col[component1_vessels] <- seq_along(component1_vessels) - 1L
  }

  flag_component1_keep <- as.integer(flag_uses_component1 & seq.int(0L, n_f - 1L) > 0L)
  flag_component1_col <- integer(n_f)
  component1_nonref_flags <- which(flag_component1_keep[-1L] == 1L)
  if (length(component1_nonref_flags) > 0L) {
    flag_component1_col[component1_nonref_flags + 1L] <- seq_along(component1_nonref_flags) - 1L
  }

  # ---- smooth parsing: catchability layer ----
  sm_catch <- .normalize_smoother_output(
    parse_smoothers(
      formula = formula_catchability,
      data = data_utm,
      knots = NULL,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )

  # ---- factor parsing: catchability layer ----
  fac_catch <- .normalize_factor_output(
    parse_factors(
      formula = formula_catchability,
      data = data_utm,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )
  if (.factor_has_var(fac_catch$basis_out, "vesid")) {
    stop(
      "`vesid` should not be supplied in `formula_catchability`; vessel effects are included by default.",
      call. = FALSE
    )
  }

  # ---- smooth parsing: population layer (obs + projection) ----
  sm_pop_obs <- .normalize_smoother_output(
    parse_smoothers(
      formula = formula_population,
      data = data_utm,
      knots = NULL,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )

  fac_pop_obs <- .normalize_factor_output(
    parse_factors(
      formula = formula_population,
      data = data_utm,
      newdata = NULL,
      basis_prev = NULL
    ),
    n_rows = n_i
  )
  if (.factor_has_var(fac_pop_obs$basis_out, "vesid")) {
    stop(
      "`vesid` should not be supplied in `formula_population`; vessel effects belong to catchability and are included by default.",
      call. = FALSE
    )
  }
  if (.factor_has_var(fac_pop_obs$basis_out, "tid", kind = "fixed")) {
    stop(
      "`tid` is already included by default as a fixed time effect. Remove `tid` from `formula_population`, or use `f(tid)` to request a random time effect.",
      call. = FALSE
    )
  }
  has_pop_intercept <- as.integer(.factor_has_random_var(fac_pop_obs$basis_out, "tid"))

  if (isTRUE(sm_pop_obs$has_smooths) ||
      isTRUE(fac_pop_obs$has_factor_fixed) ||
      isTRUE(fac_pop_obs$has_factor_random)) {
    needed_vars_pop <- unique(c(
      .smooth_vars_from_basis(sm_pop_obs$basis_out),
      .factor_vars_from_basis(fac_pop_obs$basis_out)
    ))
    projection_data_use <- if (is.null(projection_data)) {
      .default_population_projection_data(
        data_utm = data_utm,
        key = key,
        basis_out = sm_pop_obs$basis_out,
        tid_values = tid_values
      )
    } else {
      .align_projection_data_to_key(
        projection_data = projection_data,
        key = key,
        needed_vars = needed_vars_pop,
        tid_values = tid_values
      )
    }

    .warn_projection_na(
      projection_data = projection_data_use,
      needed_vars = needed_vars_pop
    )

    sm_pop_proj <- .normalize_smoother_output(
      parse_smoothers(
        formula = formula_population,
        data = data_utm,
        knots = NULL,
        newdata = projection_data_use,
        basis_prev = sm_pop_obs$basis_out
      ),
      n_rows = n_g * n_t
    )

    fac_pop_proj <- .normalize_factor_output(
      parse_factors(
        formula = formula_population,
        data = data_utm,
        newdata = projection_data_use,
        basis_prev = fac_pop_obs$basis_out
      ),
      n_rows = n_g * n_t
    )
  } else {
    projection_data_use <- .expand_projection_over_time(
      key[, c("utm_x_scale", "utm_y_scale"), drop = FALSE],
      tid_values = tid_values
    )
    sm_pop_proj <- .normalize_smoother_output(
      parse_smoothers(
        formula = NULL,
        data = data_utm,
        knots = NULL,
        newdata = projection_data_use,
        basis_prev = NULL
      ),
      n_rows = n_g * n_t
    )
    fac_pop_proj <- .normalize_factor_output(
      parse_factors(
        formula = NULL,
        data = data_utm,
        newdata = projection_data_use,
        basis_prev = NULL
      ),
      n_rows = n_g * n_t
    )
  }

  if (split_component1_design) {
    formula_catchability_1_info <- .component1_supported_formula(
      formula_catchability, data_utm, obs_component1_rows
    )
    formula_population_1_info <- .component1_supported_formula(
      formula_population, data_utm, obs_component1_rows
    )
    dropped_c1_terms <- unique(c(
      formula_catchability_1_info$dropped,
      formula_population_1_info$dropped
    ))
    if (length(dropped_c1_terms) > 0L) {
      warning(
        "Dropped component-1 formula term(s) without enough delta/component-1 support: ",
        paste(dropped_c1_terms, collapse = ", "),
        call. = FALSE
      )
    }

    formula_catchability_1 <- formula_catchability_1_info$formula
    formula_population_1 <- formula_population_1_info$formula

    data_c1_catch <- .mask_formula_vars(data_utm, formula_catchability_1, obs_component1_rows)
    sm_catch_1 <- .normalize_smoother_output(
      parse_smoothers(
        formula = formula_catchability_1,
        data = data_c1_catch,
        knots = NULL,
        newdata = NULL,
        basis_prev = NULL
      ),
      n_rows = n_i
    )
    fac_catch_1 <- .normalize_factor_output(
      parse_factors(
        formula = formula_catchability_1,
        data = data_c1_catch,
        newdata = NULL,
        basis_prev = NULL
      ),
      n_rows = n_i
    )

    data_c1_pop <- .mask_formula_vars(data_utm, formula_population_1, obs_component1_rows)
    sm_pop_i_1 <- .normalize_smoother_output(
      parse_smoothers(
        formula = formula_population_1,
        data = data_c1_pop,
        knots = NULL,
        newdata = NULL,
        basis_prev = NULL
      ),
      n_rows = n_i
    )
    fac_pop_i_1 <- .normalize_factor_output(
      parse_factors(
        formula = formula_population_1,
        data = data_c1_pop,
        newdata = NULL,
        basis_prev = NULL
      ),
      n_rows = n_i
    )
    sm_pop_g_1 <- .normalize_smoother_output(
      parse_smoothers(
        formula = formula_population_1,
        data = data_c1_pop,
        knots = NULL,
        newdata = projection_data_use,
        basis_prev = sm_pop_i_1$basis_out
      ),
      n_rows = n_g * n_t
    )
    fac_pop_projection_1 <- .mask_unseen_factor_levels(projection_data_use, fac_pop_i_1$basis_out)
    fac_pop_g_1 <- .normalize_factor_output(
      parse_factors(
        formula = formula_population_1,
        data = data_c1_pop,
        newdata = fac_pop_projection_1,
        basis_prev = fac_pop_i_1$basis_out
      ),
      n_rows = n_g * n_t
    )
  } else {
    sm_catch_1 <- .normalize_smoother_output(parse_smoothers(NULL, data_utm), n_rows = n_i)
    fac_catch_1 <- .normalize_factor_output(parse_factors(NULL, data_utm), n_rows = n_i)
    sm_pop_i_1 <- .normalize_smoother_output(parse_smoothers(NULL, data_utm), n_rows = n_i)
    sm_pop_g_1 <- .normalize_smoother_output(parse_smoothers(NULL, data_utm, newdata = projection_data_use), n_rows = n_g * n_t)
    fac_pop_i_1 <- .normalize_factor_output(parse_factors(NULL, data_utm), n_rows = n_i)
    fac_pop_g_1 <- .normalize_factor_output(parse_factors(NULL, data_utm, newdata = projection_data_use), n_rows = n_g * n_t)
  }

  # ---- build has_tf: n_t x (n_f-1), for flag-specific time effects ----
  has_tf <- matrix(FALSE, nrow = n_t, ncol = max(0L, n_f - 1L))
  if (n_f > 1L) {
    ii <- which(f_i > 0L)
    if (length(ii) > 0L) {
      tt <- t_i[ii] + 1L
      ff <- f_i[ii]
      has_tf[cbind(tt, ff)] <- TRUE
    }
  }

  has_tf_1 <- has_tf
  if (ncol(has_tf_1) > 0L) {
    for (j in seq_len(ncol(has_tf_1))) {
      if (!flag_uses_component1[j + 1L]) has_tf_1[, j] <- FALSE
    }
  }
  flag_t_constraint_1 <- .build_flag_t_constraint(has_tf_1)
  flag_t_constraint_2 <- .build_flag_t_constraint(has_tf)

  data <- list(
    n_i = n_i,
    n_t = n_t,
    n_v = n_v,
    n_f = n_f,
    n_g = n_g,
    n_region = length(region_info$region_labels),

    b_i = data_utm$cpue,
    e_i = data_utm$encounter,
    t_i = t_i,
    v_i = v_i,
    f_i = f_i,

    has_component1 = as.integer(has_component1),
    all_obs_use_component1 = as.integer(has_component1 && all(obs_uses_component1)),
    use_split_component1_design = as.integer(split_component1_design),
    obs_uses_component1 = as.integer(obs_uses_component1),
    v_component1_keep = as.integer(v_component1_keep),
    v_component1_col = as.integer(v_component1_col),
    flag_component1_keep = as.integer(flag_component1_keep),
    flag_component1_col = as.integer(flag_component1_col),

    has_tf = has_tf * 1L,
    flag_t_index_1 = matrix(
      as.integer(flag_t_constraint_1$flag_t_index),
      nrow = nrow(flag_t_constraint_1$flag_t_index),
      ncol = ncol(flag_t_constraint_1$flag_t_index)
    ),
    flag_t_index_2 = matrix(
      as.integer(flag_t_constraint_2$flag_t_index),
      nrow = nrow(flag_t_constraint_2$flag_t_index),
      ncol = ncol(flag_t_constraint_2$flag_t_index)
    ),
    n_flag_t_free_1 = as.integer(flag_t_constraint_1$n_free),
    n_flag_t_free_2 = as.integer(flag_t_constraint_2$n_free),

    area_g = key$area_km2_scaled,
    region_g = region_info$region_g,

    A_is = A_is,
    A_gs = A_gs,
    Ais_ij = Ais_ij,
    Ais_x = Ais_x,

    matern_range = diff(range(mesh$loc[, 1])) / 5,
    range_prob = 0.5,
    matern_sigma_0 = 1,
    matern_sigma_t = 1,
    matern_sigma_flag = 1,
    sigma_prob = 0.05,

    has_smooths_catch = as.integer(isTRUE(sm_catch$has_smooths)),
    Xs_catch = sm_catch$Xs,
    Zs_catch = sm_catch$Zs,
    b_smooth_start_catch = as.integer(sm_catch$b_smooth_start),
    has_fixed_factors_catch = as.integer(isTRUE(fac_catch$has_factor_fixed)),
    has_random_factors_catch = as.integer(isTRUE(fac_catch$has_factor_random)),
    Xf_catch = fac_catch$Xf,
    Zf_catch = fac_catch$Zf,
    b_factor_start_catch = as.integer(fac_catch$rand_start),

    has_smooths_catch_1 = as.integer(split_component1_design && isTRUE(sm_catch_1$has_smooths)),
    Xs_catch_1 = sm_catch_1$Xs,
    Zs_catch_1 = sm_catch_1$Zs,
    b_smooth_start_catch_1 = as.integer(sm_catch_1$b_smooth_start),
    has_fixed_factors_catch_1 = as.integer(split_component1_design && isTRUE(fac_catch_1$has_factor_fixed)),
    has_random_factors_catch_1 = as.integer(split_component1_design && isTRUE(fac_catch_1$has_factor_random)),
    Xf_catch_1 = fac_catch_1$Xf,
    Zf_catch_1 = fac_catch_1$Zf,
    b_factor_start_catch_1 = as.integer(fac_catch_1$rand_start),

    has_smooths_pop = as.integer(isTRUE(sm_pop_obs$has_smooths)),
    Xs_pop_i = sm_pop_obs$Xs,
    Zs_pop_i = sm_pop_obs$Zs,
    Xs_pop_g = sm_pop_proj$Xs,
    Zs_pop_g = sm_pop_proj$Zs,
    b_smooth_start_pop = as.integer(sm_pop_obs$b_smooth_start),
    has_fixed_factors_pop = as.integer(isTRUE(fac_pop_obs$has_factor_fixed)),
    has_random_factors_pop = as.integer(isTRUE(fac_pop_obs$has_factor_random)),
    has_pop_intercept = has_pop_intercept,
    Xf_pop_i = fac_pop_obs$Xf,
    Zf_pop_i = fac_pop_obs$Zf,
    Xf_pop_g = fac_pop_proj$Xf,
    Zf_pop_g = fac_pop_proj$Zf,
    b_factor_start_pop = as.integer(fac_pop_obs$rand_start),

    has_smooths_pop_1 = as.integer(split_component1_design && isTRUE(sm_pop_i_1$has_smooths)),
    Xs_pop_i_1 = sm_pop_i_1$Xs,
    Zs_pop_i_1 = sm_pop_i_1$Zs,
    Xs_pop_g_1 = sm_pop_g_1$Xs,
    Zs_pop_g_1 = sm_pop_g_1$Zs,
    b_smooth_start_pop_1 = as.integer(sm_pop_i_1$b_smooth_start),
    has_fixed_factors_pop_1 = as.integer(split_component1_design && isTRUE(fac_pop_i_1$has_factor_fixed)),
    has_random_factors_pop_1 = as.integer(split_component1_design && isTRUE(fac_pop_i_1$has_factor_random)),
    Xf_pop_i_1 = fac_pop_i_1$Xf,
    Zf_pop_i_1 = fac_pop_i_1$Zf,
    Xf_pop_g_1 = fac_pop_g_1$Xf,
    Zf_pop_g_1 = fac_pop_g_1$Zf,
    b_factor_start_pop_1 = as.integer(fac_pop_i_1$rand_start)
  )

  data$spde <- .prep_anisotropy(mesh = mesh, spde = spde)

  list(
    data = data,
    key = key,
    scales = list(area_scale = key_out$area_scale_val),
    region_labels = region_info$region_labels,
    projection_data = projection_data_use,
    component1_info = list(
      flag_uses_component1 = flag_uses_component1,
      split_component1_design = split_component1_design,
      n_v_component1 = length(component1_vessels),
      n_f_component1 = length(component1_nonref_flags),
      flag_t_constraint_1 = flag_t_constraint_1,
      flag_t_constraint_2 = flag_t_constraint_2
    ),
    smooth_basis = list(
      catchability = sm_catch$basis_out,
      population = sm_pop_obs$basis_out
    ),
    smooth_info = list(
      catchability = list(
        labels = sm_catch$labels,
        classes = sm_catch$classes,
        sm_dims = sm_catch$sm_dims,
        b_smooth_start = sm_catch$b_smooth_start,
        K_smooth = ncol(sm_catch$Xs),
        n_smooth = length(sm_catch$Zs),
        factor_terms = fac_catch$basis_out,
        K_factor = ncol(fac_catch$Xf),
        n_factor_random = length(fac_catch$Zf)
      ),
      population = list(
        labels = sm_pop_obs$labels,
        classes = sm_pop_obs$classes,
        sm_dims = sm_pop_obs$sm_dims,
        b_smooth_start = sm_pop_obs$b_smooth_start,
        K_smooth = ncol(sm_pop_obs$Xs),
        n_smooth = length(sm_pop_obs$Zs),
        factor_terms = fac_pop_obs$basis_out,
        K_factor = ncol(fac_pop_obs$Xf),
        n_factor_random = length(fac_pop_obs$Zf)
      )
    )
  )
}
