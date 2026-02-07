#' Prepare data objects and mesh for intCPUE workflows
#'
#' Main data-prep function for intCPUE.
#' - mesh/SPDE/A matrices
#' - extrapolation key grid + areas
#' - parse smoothers (mgcv s()) into Xs/Zs
#'
#' @param formula A model formula. Smooth terms must use s().
#' @param data_utm A data.frame containing required columns (with utm_x/y_scale)
#' @param mesh intCPUEmesh built from make_mesh(), or a custom mesh
#' @param area_scale Numeric or "auto". Scaling factor for area_km2.
#'
#' @return A list with elements mesh, data, key, scales, smooth_basis, smooth_info.
#' @author Rujia Bi \email{rbi@@iattc.org}
#' @export
make_data <- function(
    formula,
    data_utm,
    mesh,
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
  
  # mesh can be:
  #  - intCPUEmesh from make_mesh()
  #  - a bare fmesher mesh object (mesh$mesh)
  mesh_in <- mesh
  mesh_obj <- .as_intCPUEmesh(
    mesh = mesh_in,
    loc_xy = loc_xy,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    recompute_A = "auto"
  )
  
  if (!is.null(mesh_obj$loc_xy)) {
    r1 <- range(mesh_obj$loc_xy[,1]); r2 <- range(loc_xy[,1])
    if (is.finite(r1[1]) && is.finite(r2[1])) {
      if (abs(diff(r1) - diff(r2)) / max(1e-12, diff(r2)) > 0.5) {
        warning("Mesh coordinate scale may not match `utm_x_scale/utm_y_scale`. Check scaling.", call. = FALSE)
      }
    }
  }
  
  mesh <- mesh_obj$mesh
  
  # ---- SPDE (INLA) ----
  spde <- mesh_obj$spde
  
  # ---- A matrices ----
  A_is <- mesh_obj$A
  A_isT <- methods::as(A_is, "TsparseMatrix")
  Ais_ij <- cbind(A_isT@i, A_isT@j)
  Ais_x  <- A_is@x
  
  # ---- key/extrapolation grid ----
  key_out <- .prep_key_area(data_utm, mesh, area_scale = area_scale)
  key <- key_out$key
  A_gs <- key_out$A_gs
  
  # ---- smooth parsing ----
  # parse_smoothers() should:
  # - keep nrow(data) unchanged (na.pass internally or your own NA->0 logic)
  # - return $Xs (matrix), $Zs (list of sparse matrices), $sm_dims, $b_smooth_start
  sm <- parse_smoothers(
    formula = formula,
    data    = data_utm,
    knots   = NULL,
    newdata = NULL,
    basis_prev = NULL
  )
  
  n_i <- nrow(data_utm)
  
  if (!isTRUE(sm$has_smooths)) {
    Xs <- matrix(0, nrow = n_i, ncol = 0L)
    Zs <- list()
    sm_dims <- integer(0)
    b_smooth_start <- integer(0)
  } else {
    Xs <- sm$Xs
    Zs <- sm$Zs
    sm_dims <- as.integer(sm$sm_dims)
    b_smooth_start <- as.integer(sm$b_smooth_start)
    
    # defensive: ensure correct row count
    if (!identical(nrow(Xs), n_i)) {
      stop("parse_smoothers() returned Xs with nrow != nrow(data).", call. = FALSE)
    }
    if (length(Zs)) {
      for (k in seq_along(Zs)) {
        if (!identical(nrow(Zs[[k]]), n_i)) {
          stop("parse_smoothers() returned Zs[[k]] with nrow != nrow(data).", call. = FALSE)
        }
        # ensure sparse matrix class TMB likes
        if (!inherits(Zs[[k]], "sparseMatrix")) {
          Zs[[k]] <- Matrix::Matrix(Zs[[k]], sparse = TRUE)
        }
      }
    }
  }
  
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
  
  # ---- build has_tf: n_t x (n_f-1), for flag-specific time effects ----
  has_tf <- matrix(FALSE, nrow = n_t, ncol = max(0L, n_f - 1L))
  if (n_f > 1L) {
    ii <- which(f_i > 0L)  # exclude baseline 0
    if (length(ii) > 0L) {
      tt <- t_i[ii] + 1L   # R matrix rows are 1-based
      ff <- f_i[ii]        # 1..(n_f-1), columns are 1-based already
      has_tf[cbind(tt, ff)] <- TRUE
    }
  }
  
  # ---- assemble data list for TMB ----
  data <- list(
    n_i = n_i,
    n_t = n_t,
    n_v = n_v,
    n_f = n_f,
    n_g = nrow(key),
    
    b_i = data_utm$cpue,
    e_i = data_utm$encounter,
    t_i = t_i,
    v_i = v_i,
    f_i = f_i,
    
    has_tf = has_tf * 1L,   # logical -> integer (for DATA_IMATRIX)
    
    area_g = key$area_km2_scaled,
    
    A_is   = A_is,
    A_gs   = A_gs,
    Ais_ij = Ais_ij,
    Ais_x  = Ais_x,
    
    # PC priors
    matern_range   = diff(range(mesh$loc[, 1])) / 5,
    range_prob     = 0.5,
    matern_sigma_0 = 1,
    matern_sigma_t = 1,
    matern_sigma_flag = 1,
    sigma_prob     = 0.05,
    
    # smoothers
    has_smooths    = as.integer(isTRUE(sm$has_smooths)),
    Xs             = Xs,
    Zs             = Zs,
    sm_dims        = sm_dims,
    b_smooth_start = b_smooth_start
  )
  
  data$spde <- .prep_anisotropy(mesh = mesh, spde = spde)
  
  list(
    data = data,
    key = key,
    scales = list(area_scale = key_out$area_scale_val),
    
    # smooth outputs for plotting/prediction
    smooth_basis = sm$basis_out,
    smooth_info = list(
      labels = sm$labels,
      classes = sm$classes,
      sm_dims = sm_dims,
      b_smooth_start = b_smooth_start,
      K_smooth = ncol(Xs),
      n_smooth = length(Zs)
    )
  )
}
