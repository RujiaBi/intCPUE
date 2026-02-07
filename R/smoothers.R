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

# from mgcv docs ?mgcv::smooth2random
s2rPred <- function(sm, re, data) {
  ## Function to aid prediction from smooths represented as type==2
  ## random effects. re must be the result of smooth2random(sm,...,type=2).
  if (!all(sm$term %in% colnames(data))) {
    stop(paste("A smoother term is missing from 'newdata':",
               sm$term[!sm$term %in% colnames(data)]),
         call. = FALSE
    )
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

# ---- NA->0 wrapper around PredictMat + smooth2random output -----------------
# Goal: if any smooth covariate involved in sm is NA in a row, then
# that row contributes 0 to BOTH fixed part (Xf) and random part (rand) matrices.
.zero_rows_for_na <- function(sm, data, Xf, rand_list) {
  # sm$term contains all columns needed to evaluate the smoother (incl by vars)
  # Identify rows where ANY required term is NA:
  miss <- rep(FALSE, nrow(data))
  for (nm in sm$term) {
    if (!nm %in% names(data)) {
      stop("Internal error: smoother term '", nm, "' not found in data.", call. = FALSE)
    }
    miss <- miss | is.na(data[[nm]])
  }
  if (any(miss)) {
    if (!is.null(Xf) && nrow(Xf) > 0L) Xf[miss, ] <- 0
    if (length(rand_list)) {
      for (k in seq_along(rand_list)) {
        if (!is.null(rand_list[[k]]) && nrow(rand_list[[k]]) > 0L) {
          rand_list[[k]][miss, ] <- 0
        }
      }
    }
  }
  list(Xf = Xf, rand = rand_list, miss = miss)
}

# ---- Core parser: only s(), plus NA->0 -------------------
# Returns:
#   Xs: combined fixed-effect smoother design
#   Zs: list of random-effect (penalized) basis matrices
#   has_smooths, labels, classes, basis_out, sm_dims, b_smooth_start
#
# Notes:
# - keep mgcv::smoothCon(... absorb.cons=TRUE, modCon=3, diagonal.penalty=FALSE).
# - NA handling:
#   * For training data: we must feed mgcv something without NAs, so we
#     temporarily replace NAs with safe values before calling smoothCon/smooth2random,
#     then zero-out the corresponding rows in Xf/rand afterwards.
#   * For newdata prediction: same idea, but we must have basis_prev and we call s2rPred
#     which uses PredictMat internally; we again zero-out rows after.
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
  
  if (length(smooth_i) > 0) {
    has_smooths <- TRUE
    smterms <- terms[smooth_i]
    ns <- 0
    ns_Xf <- 0
    
    for (i in seq_along(smterms)) {
      if (grepl('bs\\=\\"re', smterms[i])) stop("bs = 're' is not currently supported for smooths", call. = FALSE)
      if (grepl('fx\\=T', smterms[i])) stop("fx = TRUE is not currently supported for smooths", call. = FALSE)
      
      expr <- str2expression(smterms[i])[[1]]
      eval_env <- new.env(parent = baseenv())
      eval_env$s <- mgcv::s
      obj <- eval(expr, envir = eval_env)
      labels[[i]] <- obj$label
      classes[[i]] <- attr(obj, "class")
      
      
      # Choose which data we're evaluating the matrices on:
      eval_data <- if (is.null(newdata)) data else newdata
      
      # ---- NA pre-fill to allow mgcv to build matrices ----
      # We must ensure eval_data has no NA in variables needed for the smooth term.
      # Strategy: for any variable used by obj (obj$term), replace NA with median (numeric)
      # or first level (factor/character). This is ONLY to let PredictMat run;
      # we will then zero-out affected rows.
      eval_data_safe <- eval_data
      needed <- obj$term
      for (nm in needed) {
        if (!nm %in% names(eval_data_safe)) {
          stop("Smooth term requires variable '", nm, "' but it is not in data.", call. = FALSE)
        }
        if (anyNA(eval_data_safe[[nm]])) {
          x <- eval_data_safe[[nm]]
          if (is.numeric(x)) {
            fill <- stats::median(x, na.rm = TRUE)
            if (!is.finite(fill)) fill <- 0
            x[is.na(x)] <- fill
            eval_data_safe[[nm]] <- x
          } else if (is.factor(x)) {
            # choose first existing level that appears (or first level)
            lev <- levels(x)
            if (length(lev) == 0) stop("Factor '", nm, "' has no levels.", call. = FALSE)
            fill <- lev[1]
            x <- as.character(x)
            x[is.na(x)] <- fill
            eval_data_safe[[nm]] <- factor(x, levels = lev)
          } else {
            # character / logical / etc.
            fill <- x[which(!is.na(x))[1]]
            if (is.na(fill)) fill <- ""
            x[is.na(x)] <- fill
            eval_data_safe[[nm]] <- x
          }
        }
      }
      
      if (is.null(newdata)) {
        # training: build basis from training data (safe-filled)
        basis_out[[i]] <- basis[[i]] <- mgcv::smoothCon(
          object = obj, data = eval_data_safe,
          knots = knots, absorb.cons = TRUE, modCon = 3,
          diagonal.penalty = FALSE
        )
      } else {
        # prediction: basis must be provided from training
        if (is.null(basis_prev)) stop("basis_prev must be provided when newdata is not NULL.", call. = FALSE)
        basis[[i]] <- basis_prev[[i]]
      }
      
      # For each element (multiple for by-terms)
      for (j in seq_along(basis[[i]])) {
        ns_Xf <- ns_Xf + 1
        
        rasm <- mgcv::smooth2random(basis[[i]][[j]], names(data), type = 2)
        
        if (!is.null(newdata)) {
          # build prediction matrices for newdata (safe-filled) then zero out NA rows
          rasm <- s2rPred(basis[[i]][[j]], rasm, eval_data_safe)
        } else {
          # For training we already built basis on safe data; rasm$Xf and rasm$rand correspond to eval_data_safe
          # (which is data_safe). That's fine because we will zero NA rows next.
        }
        
        # Apply NA->0 row logic using the ORIGINAL (possibly NA) eval_data
        tmp <- .zero_rows_for_na(basis[[i]][[j]], eval_data, rasm$Xf, rasm$rand)
        Xf <- tmp$Xf
        rand_list <- tmp$rand
        
        # Collect rand matrices (elements > 1 with if s(x, y))
        for (k in seq_along(rand_list)) {
          ns <- ns + 1
          Zs[[ns]] <- rand_list[[k]]
        }
        Xs[[ns_Xf]] <- Xf
      }
    }
    
    sm_dims <- unlist(lapply(Zs, ncol))
    Xs <- do.call(cbind, Xs)
    b_smooth_start <- c(0, cumsum(sm_dims)[-length(sm_dims)])
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