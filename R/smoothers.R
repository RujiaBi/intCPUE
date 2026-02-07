# ------------------------------------------------------------------------------
# Adapted from sdmTMB smoother parsing workflow (GPL-3)
# https://github.com/sdmTMB/sdmTMB
#
# Includes helper utilities based on:
#   brms internal functions (GPL-3)
#   mgcv smooth2random documentation example (GPL-2)
#
# Modified in intCPUE to support NA-safe smooth evaluation.
# ------------------------------------------------------------------------------

# ---- helpers: complete-case build + expand back to full rows -----------------

# Expand a matrix back to full n rows, fill missing rows with 0
.expand_back <- function(M_cc, keep, n_all) {
  if (is.null(M_cc)) return(NULL)
  if (is.vector(M_cc)) M_cc <- matrix(M_cc, ncol = 1)
  M_all <- matrix(0, nrow = n_all, ncol = ncol(M_cc))
  if (nrow(M_cc) > 0L && any(keep)) M_all[keep, ] <- M_cc
  M_all
}

# For a given smooth object sm (or mgcv smooth), find complete rows required for evaluation
.complete_rows_for_sm <- function(sm, data) {
  needed <- sm$term
  if (!is.null(sm$by) && is.character(sm$by) && length(sm$by) == 1L &&
      nzchar(sm$by) && sm$by %in% names(data)) {
    needed <- unique(c(needed, sm$by))
  }
  keep <- rep(TRUE, nrow(data))
  for (nm in needed) {
    if (!nm %in% names(data)) {
      stop("Smooth term requires variable '", nm, "' but it is not in data.", call. = FALSE)
    }
    keep <- keep & !is.na(data[[nm]])
  }
  keep
}

# ---- brms helpers ------------------------------------------------------------

# from brms:::rm_wsp()
rm_wsp <- function (x) {
  out <- gsub("[ \t\r\n]+", "", x, perl = TRUE)
  dim(out) <- dim(x)
  out
}

# from brms:::all_terms()
all_terms <- function (x) {
  if (!length(x)) {
    return(character(0))
  }
  if (!inherits(x, "terms")) {
    x <- stats::terms(stats::as.formula(x))
  }
  rm_wsp(attr(x, "term.labels"))
}

get_smooth_terms <- function(terms) {
  grep("s\\(", terms)
}

# ---- mgcv prediction helper --------------------------------------------------

# from mgcv docs ?mgcv::smooth2random
s2rPred <- function(sm, re, data) {
  ## Function to aid prediction from smooths represented as type==2
  ## random effects. re must be the result of smooth2random(sm,...,type=2).
  needed <- sm$term
  if (!is.null(sm$by) && is.character(sm$by) && length(sm$by) == 1L &&
      nzchar(sm$by) && sm$by %in% colnames(data)) {
    needed <- unique(c(needed, sm$by))
  }
  if (!all(needed %in% colnames(data))) {
    miss <- needed[!needed %in% colnames(data)]
    stop(paste("A smoother term is missing from 'newdata':",
               paste(miss, collapse = ", ")),
         call. = FALSE)
  }
  X <- mgcv::PredictMat(sm, data) ## get prediction matrix for new data
  ## transform to r.e. parameterization
  if (!is.null(re$trans.U)) {
    X <- X %*% re$trans.U
  }
  X <- t(t(X) * re$trans.D)
  ## re-order columns according to random effect re-ordering...
  X[, re$rind] <- X[, re$pen.ind != 0]
  ## re-order penalization index in same way
  pen.ind <- re$pen.ind
  pen.ind[re$rind] <- pen.ind[pen.ind > 0]
  ## start return object...
  r <- list(rand = list(), Xf = X[, which(re$pen.ind == 0), drop = FALSE])
  for (i in seq_along(re$rand)) { ## loop over random effect matrices
    r$rand[[i]] <- X[, which(pen.ind == i), drop = FALSE]
    attr(r$rand[[i]], "s.label") <- attr(re$rand[[i]], "s.label")
  }
  names(r$rand) <- names(re$rand)
  r
}

# ---- Core parser: only s(), complete-case basis build + NA rows -> 0 ---------
# Returns:
#   Xs: combined fixed-effect smoother design
#   Zs: list of random-effect (penalized) basis matrices
#   has_smooths, labels, classes, basis_out, sm_dims, b_smooth_start
#
# Notes:
# - keep mgcv::smoothCon(... absorb.cons=TRUE, modCon=3, diagonal.penalty=FALSE).
# - NA handling:
#   * For training: build each smooth basis ONLY on complete rows for that smooth.
#     Then expand Xf/rand back to full nrow(data) with missing rows = 0.
#   * For prediction: compute Xf/rand ONLY for complete rows in newdata (per smooth),
#     then expand back with missing rows = 0.
parse_smoothers <- function(formula, data, knots = NULL,
                            newdata = NULL, basis_prev = NULL) {
  
  terms <- all_terms(formula)
  smooth_i <- get_smooth_terms(terms)
  
  basis <- list()
  basis_out <- list()
  Zs <- list()
  Xs <- list()
  labels <- list()
  classes <- list()
  
  # which dataset we're evaluating on (training vs prediction)
  eval_data <- if (is.null(newdata)) data else newdata
  n_all <- nrow(eval_data)
  
  if (length(smooth_i) > 0) {
    has_smooths <- TRUE
    smterms <- terms[smooth_i]
    ns <- 0
    ns_Xf <- 0
    
    for (i in seq_along(smterms)) {
      if (grepl('bs\\=\\"re', smterms[i])) {
        stop("bs = 're' is not currently supported for smooths", call. = FALSE)
      }
      if (grepl('fx\\=T', smterms[i])) {
        stop("fx = TRUE is not currently supported for smooths", call. = FALSE)
      }
      
      expr <- str2expression(smterms[i])[[1]]
      eval_env <- new.env(parent = baseenv())
      eval_env$s <- mgcv::s
      obj <- eval(expr, envir = eval_env)
      
      labels[[i]] <- obj$label
      classes[[i]] <- attr(obj, "class")
      
      # complete rows needed for this smooth term (based on obj$term)
      keep_i <- .complete_rows_for_sm(obj, eval_data)
      
      if (is.null(newdata)) {
        # TRAINING: build basis ONLY on complete rows for this smooth
        data_cc <- eval_data[keep_i, , drop = FALSE]
        
        # guard against too few rows
        if (nrow(data_cc) < 5L) {
          stop("Too few non-missing rows (", nrow(data_cc),
               ") to build smooth term: ", labels[[i]], call. = FALSE)
        }
        
        basis_out[[i]] <- basis[[i]] <- mgcv::smoothCon(
          object = obj, data = data_cc,
          knots = knots, absorb.cons = TRUE, modCon = 3,
          diagonal.penalty = FALSE
        )
      } else {
        # PREDICTION: use basis from training
        if (is.null(basis_prev) || length(basis_prev) < length(smterms)) {
          stop("basis_prev must be provided and have at least one entry per smooth term.", call. = FALSE)
        }
        
        basis[[i]] <- basis_prev[[i]]
      }
      
      # For each element (multiple for by-terms)
      for (j in seq_along(basis[[i]])) {
        ns_Xf <- ns_Xf + 1
        
        if (is.null(newdata)) {
          # TRAINING: smooth2random on cc rows, then expand back
          data_cc <- eval_data[keep_i, , drop = FALSE]
          
          rasm_cc <- mgcv::smooth2random(basis[[i]][[j]], names(data), type = 2)
          
          Xf <- .expand_back(rasm_cc$Xf, keep_i, n_all)
          rand_list <- lapply(rasm_cc$rand, .expand_back, keep = keep_i, n_all = n_all)
          
        } else {
          # PREDICTION: predict only on cc rows of newdata, then expand back
          new_cc <- eval_data[keep_i, , drop = FALSE]
          
          # baseline re-parameterization object; must be created using training columns
          # names(data) here refers to training data column names expected by mgcv
          rasm0 <- mgcv::smooth2random(basis[[i]][[j]], names(data), type = 2)
          
          # get matrices for complete rows only
          rasm_cc <- s2rPred(basis[[i]][[j]], rasm0, new_cc)
          
          Xf <- .expand_back(rasm_cc$Xf, keep_i, n_all)
          rand_list <- lapply(rasm_cc$rand, .expand_back, keep = keep_i, n_all = n_all)
        }
        
        # Collect rand matrices (elements > 1 with if s(x, y))
        for (k in seq_along(rand_list)) {
          if (is.null(rand_list[[k]])) next
          ns <- ns + 1
          Zs[[ns]] <- rand_list[[k]]
        }
        Xs[[ns_Xf]] <- Xf
      }
    }
    
    sm_dims <- unlist(lapply(Zs, ncol))
    Xs <- do.call(cbind, Xs)
    
    # robust start index (handle case with zero smooth random cols)
    if (length(sm_dims) > 0L) {
      b_smooth_start <- c(0, cumsum(sm_dims)[-length(sm_dims)])
    } else {
      b_smooth_start <- 0L
    }
    
  } else {
    has_smooths <- FALSE
    sm_dims <- 0L
    b_smooth_start <- 0L
    Xs <- matrix(nrow = 0L, ncol = 0L)
  }
  
  list(
    Xs = Xs,
    Zs = Zs,
    has_smooths = has_smooths,
    labels = labels,
    classes = classes,
    basis_out = basis_out,
    sm_dims = sm_dims,
    b_smooth_start = b_smooth_start
  )
}
