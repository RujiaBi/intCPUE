.expand_back_sparse <- function(M_cc, keep, n_all) {
  if (is.null(M_cc)) {
    return(NULL)
  }
  if (!inherits(M_cc, "sparseMatrix")) {
    M_cc <- Matrix::Matrix(M_cc, sparse = TRUE)
  }
  M_all <- Matrix::Matrix(0, nrow = n_all, ncol = ncol(M_cc), sparse = TRUE)
  if (nrow(M_cc) > 0L && any(keep)) {
    M_all[keep, ] <- M_cc
  }
  M_all
}

.extract_formula_terms <- function(formula) {
  terms <- all_terms(formula)
  if (!length(terms)) {
    return(list(
      all = character(0),
      smooth = character(0),
      factor_random = character(0),
      factor_fixed = character(0)
    ))
  }

  smooth_i <- get_smooth_terms(terms)
  smooth_terms <- terms[smooth_i]
  non_smooth <- if (length(smooth_i)) terms[-smooth_i] else terms

  random_i <- grepl("^f\\(", non_smooth)
  factor_random <- non_smooth[random_i]
  factor_fixed <- non_smooth[!random_i]

  list(
    all = terms,
    smooth = smooth_terms,
    factor_random = factor_random,
    factor_fixed = factor_fixed
  )
}

.parse_factor_var <- function(term, random = FALSE) {
  if (random) {
    m <- regexec("^f\\(([[:alnum:]_.]+)\\)$", term)
    reg <- regmatches(term, m)[[1]]
    if (length(reg) != 2L) {
      stop(
        "Random factor terms must use simple syntax like f(var). Unsupported term: ",
        term,
        call. = FALSE
      )
    }
    return(reg[2])
  }

  if (!grepl("^[[:alnum:]_.]+$", term)) {
    stop(
      "Fixed factor terms must be simple variable names like `tid`. Unsupported term: ",
      term,
      call. = FALSE
    )
  }
  term
}

.normalize_factor_levels <- function(x, term) {
  x_chr <- as.character(x)
  x_chr <- x_chr[!is.na(x_chr)]
  if (!length(x_chr)) {
    stop("Factor term `", term, "` has no non-missing values.", call. = FALSE)
  }

  if (is.factor(x)) {
    lvls <- levels(x)
    lvls[lvls %in% unique(x_chr)]
  } else {
    sort(unique(x_chr))
  }
}

.build_fixed_factor_matrix <- function(values, levels, term) {
  fac <- factor(as.character(values), levels = levels)
  mm <- stats::model.matrix(~ fac - 1)
  colnames(mm) <- paste0(term, "::", levels)
  mm
}

.build_random_factor_matrix <- function(values, levels, term) {
  fac <- factor(as.character(values), levels = levels)
  mm <- Matrix::sparse.model.matrix(~ fac - 1)
  colnames(mm) <- paste0(term, "::", levels)
  mm
}

.factor_vars_from_basis <- function(basis_out) {
  if (!length(basis_out)) {
    return(character(0))
  }
  vars <- vapply(basis_out, `[[`, "", "var")
  vars <- vars[!is.na(vars) & nzchar(vars)]
  unique(vars)
}

parse_factors <- function(formula, data, newdata = NULL, basis_prev = NULL) {
  parts <- .extract_formula_terms(formula)
  eval_data <- if (is.null(newdata)) data else newdata
  n_all <- nrow(eval_data)

  Xf_list <- list()
  Zf_list <- list()
  basis_out <- list()

  if (!length(parts$factor_fixed) && !length(parts$factor_random)) {
    return(list(
      Xf = matrix(0, nrow = n_all, ncol = 0L),
      Zf = list(),
      has_factor_fixed = FALSE,
      has_factor_random = FALSE,
      basis_out = list(),
      rand_dims = integer(0),
      rand_start = integer(0)
    ))
  }

  term_counter <- 0L

  for (term in parts$factor_fixed) {
    term_counter <- term_counter + 1L
    var <- .parse_factor_var(term, random = FALSE)
    if (!var %in% names(eval_data)) {
      stop("Fixed factor term `", term, "` requires variable `", var, "` in data.", call. = FALSE)
    }
    keep <- !is.na(eval_data[[var]])

    levels_i <- if (is.null(newdata)) {
      .normalize_factor_levels(eval_data[[var]], term)
    } else {
      if (is.null(basis_prev) || length(basis_prev) < term_counter) {
        stop("basis_prev is missing factor basis information for term `", term, "`.", call. = FALSE)
      }
      basis_prev[[term_counter]]$levels
    }

    vals_cc <- eval_data[[var]][keep]
    unseen <- setdiff(unique(as.character(vals_cc)), levels_i)
    if (length(unseen)) {
      stop(
        "Factor term `", term, "` in prediction data contains unseen level(s): ",
        paste(unseen, collapse = ", "),
        call. = FALSE
      )
    }

    X_cc <- .build_fixed_factor_matrix(vals_cc, levels_i, term)
    Xf_list[[length(Xf_list) + 1L]] <- .expand_back(X_cc, keep, n_all)
    basis_out[[term_counter]] <- list(term = term, var = var, levels = levels_i, kind = "fixed")
  }

  for (term in parts$factor_random) {
    term_counter <- term_counter + 1L
    var <- .parse_factor_var(term, random = TRUE)
    if (!var %in% names(eval_data)) {
      stop("Random factor term `", term, "` requires variable `", var, "` in data.", call. = FALSE)
    }
    keep <- !is.na(eval_data[[var]])

    levels_i <- if (is.null(newdata)) {
      .normalize_factor_levels(eval_data[[var]], term)
    } else {
      if (is.null(basis_prev) || length(basis_prev) < term_counter) {
        stop("basis_prev is missing factor basis information for term `", term, "`.", call. = FALSE)
      }
      basis_prev[[term_counter]]$levels
    }

    vals_cc <- eval_data[[var]][keep]
    unseen <- setdiff(unique(as.character(vals_cc)), levels_i)
    if (length(unseen)) {
      stop(
        "Factor term `", term, "` in prediction data contains unseen level(s): ",
        paste(unseen, collapse = ", "),
        call. = FALSE
      )
    }

    Z_cc <- .build_random_factor_matrix(vals_cc, levels_i, term)
    Zf_list[[length(Zf_list) + 1L]] <- .expand_back_sparse(Z_cc, keep, n_all)
    basis_out[[term_counter]] <- list(term = term, var = var, levels = levels_i, kind = "random")
  }

  Xf <- if (length(Xf_list)) do.call(cbind, Xf_list) else matrix(0, nrow = n_all, ncol = 0L)
  rand_dims <- if (length(Zf_list)) as.integer(vapply(Zf_list, ncol, 0L)) else integer(0)
  rand_start <- if (length(rand_dims)) c(0L, cumsum(rand_dims)[-length(rand_dims)]) else integer(0)

  list(
    Xf = Xf,
    Zf = Zf_list,
    has_factor_fixed = ncol(Xf) > 0L,
    has_factor_random = length(Zf_list) > 0L,
    basis_out = basis_out,
    rand_dims = rand_dims,
    rand_start = rand_start
  )
}

.factor_has_random_var <- function(basis_out, var) {
  if (!length(basis_out)) {
    return(FALSE)
  }
  any(vapply(
    basis_out,
    function(x) identical(x$kind, "random") && identical(x$var, var),
    logical(1)
  ))
}

.factor_has_var <- function(basis_out, var, kind = NULL) {
  if (!length(basis_out)) {
    return(FALSE)
  }
  any(vapply(
    basis_out,
    function(x) {
      same_var <- identical(x$var, var)
      same_kind <- is.null(kind) || identical(x$kind, kind)
      same_var && same_kind
    },
    logical(1)
  ))
}

.add_factor_intercept <- function(fac, n_rows, label = "(Intercept)") {
  if (!identical(nrow(fac$Xf), n_rows)) {
    stop("Factor intercept requested with incompatible row count.", call. = FALSE)
  }

  intercept_col <- matrix(1.0, nrow = n_rows, ncol = 1L)
  colnames(intercept_col) <- label
  fac$Xf <- cbind(intercept_col, fac$Xf)
  fac$has_factor_fixed <- TRUE
  fac
}
