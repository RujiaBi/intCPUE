#' Prepare data objects and mesh for intCPUE workflows
#'
#' Main data-prep function for intCPUE.
#' - UTM transform + shared scaling
#' - mesh/SPDE/A matrices
#' - extrapolation key grid + areas
#' - parse smoothers (mgcv s()) into Xs/Zs
#'
#' @param formula A model formula. Smooth terms must use s().
#' @param data_input A data.frame containing required columns.
#' @param utm_zone Optional integer. If NULL, uses the mode of per-record UTM zones.
#' @param coord_scale Numeric or "auto". Shared scaling factor for UTM x and y.
#' @param area_scale Numeric or "auto". Scaling factor for area_km2.
#' @param concave,convex,resx,resy,max_edge,bound_outer,inner_mesh_length,outer_mesh_length,cutoff_set
#'   Passed to make_mesh().
#'
#' @return A list with elements mesh, data, combined_data, key, scales, smooth_basis, smooth_info.
#' @author Rujia Bi \email{rbi@@iattc.org}
#' @export
make_data <- function(
    formula,
    data_input,
    utm_zone = NULL,
    coord_scale = "auto",
    area_scale = "auto",
    concave = -0.05,
    convex  = -0.05,
    resx = 100,
    resy = 100,
    max_edge = 0.5,
    bound_outer = 1,
    inner_mesh_length = 1,
    outer_mesh_length = 3,
    cutoff_set = 2
) {
  data_input <- as.data.frame(data_input)
  
  .check_required_cols(data_input, c("cpue", "encounter", "lon", "lat", "vesid", "tid", "flagid"))
  .check_numeric(data_input, c("cpue", "encounter", "lon", "lat", "vesid", "tid", "flagid"))
  
  if (anyNA(data_input$lon) || anyNA(data_input$lat)) {
    stop("`lon`/`lat` must not contain NA.", call. = FALSE)
  }
  
  # ---- coords -> UTM + scaling ----
  utm <- .prep_utm_scaled(data_input, utm_zone = utm_zone, coord_scale = coord_scale)
  combined_data <- utm$combined_data
  
  
  # ---- mesh ----
  mesh <- make_mesh(
    combined_data$utm_x_scale,
    combined_data$utm_y_scale,
    concave = concave,
    convex  = convex,
    resx = resx,
    resy = resy,
    max_edge = max_edge,
    bound_outer = bound_outer,
    inner_mesh_length = inner_mesh_length,
    outer_mesh_length = outer_mesh_length,
    cutoff_set = cutoff_set
  )
  
  # ---- SPDE (INLA) ----
  inla_spde <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = c(diff(range(mesh$loc[, 1])) / 5, 0.05),
    prior.sigma = c(1, 0.05)
  )
  
  # ---- A matrices ----
  A_is  <- fmesher::fm_basis(mesh, loc = as.matrix(combined_data[, c("utm_x_scale", "utm_y_scale")]))
  A_isT <- methods::as(A_is, "TsparseMatrix")
  Ais_ij <- cbind(A_isT@i, A_isT@j)
  Ais_x  <- A_is@x
  
  # ---- key/extrapolation grid ----
  key_out <- .prep_key_area(combined_data, mesh, area_scale = area_scale)
  key <- key_out$key
  A_gs <- key_out$A_gs
  
  # ---- smooth parsing (sdmTMB style: Xs + Zs list) ----
  # parse_smoothers() should:
  # - keep nrow(data) unchanged (na.pass internally or your own NA->0 logic)
  # - return $Xs (matrix), $Zs (list of sparse matrices), $sm_dims, $b_smooth_start
  sm <- parse_smoothers(
    formula = formula,
    data    = combined_data,
    knots   = NULL,
    newdata = NULL,
    basis_prev = NULL
  )
  
  n_i <- nrow(combined_data)
  
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
    
    # NA handling: rows with any NA in smooth covariates -> 0 effect
    # If the parse_smoothers already did this, this is harmless.
    if (!is.null(sm$na_rows) && length(sm$na_rows)) {
      na_rows <- sm$na_rows
      if (length(na_rows)) {
        if (ncol(Xs) > 0) Xs[na_rows, ] <- 0
        if (length(Zs)) {
          for (k in seq_along(Zs)) {
            Zs[[k]][na_rows, ] <- 0
          }
        }
      }
    }
  }
  
  # ---- assemble data list for TMB ----
  data <- list(
    n_i = n_i,
    n_t = length(unique(combined_data$tid)),
    n_v = length(unique(combined_data$vesid)),
    n_f = length(unique(combined_data$flagid)),
    n_g = nrow(key),
    
    b_i = combined_data$cpue,
    e_i = combined_data$encounter,
    t_i = combined_data$tid,
    v_i = combined_data$vesid,
    f_i = combined_data$flagid,
    
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
  
  data$spde <- .prep_anisotropy(mesh = mesh, inla_spde = inla_spde)
  
  list(
    mesh = mesh,
    data = data,
    combined_data = combined_data,
    key = key,
    scales = list(utm_scale = utm$utm_scale, area_scale = key_out$area_scale_val),
    
    # smooth outputs for plotting/prediction
    smooth_basis = sm$basis_out,
    smooth_info = list(
      labels = sm$labels,
      classes = sm$classes,
      sm_dims = sm_dims,
      b_smooth_start = b_smooth_start,
      K_smooth = ncol(Xs),
      n_smooth = length(Zs),
      n_report = n_report
    )
  )
}