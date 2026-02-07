#' Calculate equal-area size (km^2) for each lon/lat grid cell
#'
#' Computes the area (km^2) of grid cells defined by lon/lat cell centers.
#' Cells are constructed as rectangles in lon/lat degrees with width/height
#' given by `cellsize`, then projected to an equal-area CRS to compute areas.
#'
#' @param dxy A data.frame containing lon/lat grid centers. By default, the first
#'   two columns are treated as lon/lat if `lon_name`/`lat_name` are not found.
#'   Must also contain `utm_x_scale` and `utm_y_scale` (used later by the model
#'   key/extrapolation workflow).
#' @param cellsize Numeric length 1 (dx = dy) or length 2 (c(dx, dy)) in degrees,
#'   or `"auto"`. If `"auto"`, the function attempts to infer (dx, dy) from the
#'   spacing of unique lon/lat centers (using the most common spacing and a
#'   robust fallback). **This works best for regular grids. If your grid is
#'   irregular or multi-resolution, supply `cellsize=` explicitly.**
#' @param lon_name,lat_name Column names for lon/lat centers in `dxy`.
#' @param crs_equal_area A CRS string passed to `sf::st_transform()` for area
#'   calculation. Default uses Equal Earth.
#' @param tol Tolerance (degrees) for distinguishing unique spacings and for
#'   neighbor checks.
#' @param check_grid Logical; if TRUE (default), performs a lightweight sanity
#'   check that inferred `cellsize` matches local neighbor spacing and warns if
#'   inconsistent.
#' @param check_n Number of random rows used for the sanity check.
#' @param check_seed Optional seed for reproducible checking.
#' @param quiet Logical; if FALSE, message when a fallback inference method is used.
#'
#' @return A data.frame with columns `order`, `utm_x_scale`, `utm_y_scale`,
#'   and `area_km2` (aligned to the input row order).
#' @export
calculate_area <- function(
    dxy,
    cellsize = "auto",
    lon_name = "X",
    lat_name = "Y",
    crs_equal_area = "+proj=eqearth +datum=WGS84",
    tol = 1e-6,
    check_grid = TRUE,
    check_n = 200,
    check_seed = NULL,
    quiet = TRUE
) {
  dxy <- as.data.frame(dxy)
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for area calculation.", call. = FALSE)
  }
  
  # ---- accept first two columns as lon/lat if lon_name/lat_name not present ----
  if (!(lon_name %in% names(dxy)) || !(lat_name %in% names(dxy))) {
    if (ncol(dxy) < 2) {
      stop("`dxy` must have lon/lat columns (or at least 2 columns).", call. = FALSE)
    }
    names(dxy)[1:2] <- c(lon_name, lat_name)
  }
  
  # ---- required UTM columns (used later by key/A matrices) ----
  if (!("utm_x_scale" %in% names(dxy)) || !("utm_y_scale" %in% names(dxy))) {
    stop("`dxy` must include columns `utm_x_scale` and `utm_y_scale`.", call. = FALSE)
  }
  
  # ---- order column to preserve row order ----
  if (!("order" %in% names(dxy))) dxy$order <- seq_len(nrow(dxy))
  
  lon <- dxy[[lon_name]]
  lat <- dxy[[lat_name]]
  
  if (!is.numeric(lon) || !is.numeric(lat)) {
    stop("lon/lat columns must be numeric.", call. = FALSE)
  }
  if (anyNA(lon) || anyNA(lat)) {
    stop("lon/lat centers must not contain NA.", call. = FALSE)
  }
  
  # ---- cellsize handling ----
  if (is.character(cellsize) && length(cellsize) == 1L) {
    if (!identical(cellsize, "auto")) {
      stop('`cellsize` must be numeric (length 1 or 2) or "auto".', call. = FALSE)
    }
    
    # infer dx/dy: try mode first, then median
    cs_mode <- infer_cellsize_deg(lon = lon, lat = lat, tol = tol, force = "mode")
    dx <- cs_mode[1]; dy <- cs_mode[2]
    
    if (check_grid) {
      ok_mode <- .check_grid_spacing(
        lon = lon, lat = lat, dx = dx, dy = dy,
        tol = tol, n = check_n, seed = check_seed,
        return_ok = TRUE
      )
      
      if (isFALSE(ok_mode)) {
        cs_med <- infer_cellsize_deg(lon = lon, lat = lat, tol = tol, force = "median")
        dx2 <- cs_med[1]; dy2 <- cs_med[2]
        
        ok_med <- .check_grid_spacing(
          lon = lon, lat = lat, dx = dx2, dy = dy2,
          tol = tol, n = check_n, seed = check_seed,
          return_ok = TRUE
        )
        
        if (isTRUE(ok_med)) {
          dx <- dx2; dy <- dy2
          if (!quiet) {
            message(sprintf(
              "Auto cellsize inference: mode-based spacing looked inconsistent; using median-based spacing (dx=%.6g, dy=%.6g).",
              dx, dy
            ))
          }
        } else {
          # keep mode, but warn with actionable guidance
          .check_grid_spacing(
            lon = lon, lat = lat, dx = dx, dy = dy,
            tol = tol, n = check_n, seed = check_seed,
            return_ok = FALSE
          )
        }
      }
    }
    
  } else {
    if (length(cellsize) == 1L) {
      dx <- dy <- as.numeric(cellsize)
    } else if (length(cellsize) == 2L) {
      dx <- as.numeric(cellsize[1])
      dy <- as.numeric(cellsize[2])
    } else {
      stop('`cellsize` must be numeric length 1, 2, or "auto".', call. = FALSE)
    }
    
    if (check_grid) {
      .check_grid_spacing(lon = lon, lat = lat, dx = dx, dy = dy, tol = tol, n = check_n, seed = check_seed)
    }
  }
  
  if (!is.finite(dx) || !is.finite(dy) || dx <= 0 || dy <= 0) {
    stop("`cellsize` must be positive and finite.", call. = FALSE)
  }
  
  hx <- dx / 2
  hy <- dy / 2
  
  # ---- unique grid centers ----
  grid <- unique(dxy[, c(lon_name, lat_name), drop = FALSE])
  grid <- grid[order(grid[[lon_name]], grid[[lat_name]]), , drop = FALSE]
  rownames(grid) <- NULL
  
  # ---- build polygons around centers ----
  make_poly <- function(lon0, lat0) {
    coords <- rbind(
      c(lon0 - hx, lat0 - hy),
      c(lon0 + hx, lat0 - hy),
      c(lon0 + hx, lat0 + hy),
      c(lon0 - hx, lat0 + hy),
      c(lon0 - hx, lat0 - hy)
    )
    sf::st_polygon(list(coords))
  }
  
  polys <- mapply(make_poly, grid[[lon_name]], grid[[lat_name]], SIMPLIFY = FALSE)
  
  grid_sf <- sf::st_sf(
    data.frame(X = grid[[lon_name]], Y = grid[[lat_name]]),
    geometry = sf::st_sfc(polys, crs = 4326)
  )
  
  grid_eq <- sf::st_transform(grid_sf, crs_equal_area)
  area_km2 <- as.numeric(sf::st_area(grid_eq)) / 1e6
  
  # ---- map area back to dxy robustly ----
  grid_key <- paste0(.qkey(grid[[lon_name]], tol), "_", .qkey(grid[[lat_name]], tol))
  dxy_key  <- paste0(.qkey(dxy[[lon_name]],  tol), "_", .qkey(dxy[[lat_name]],  tol))
  
  if (anyDuplicated(grid_key)) {
    warning("Duplicate grid keys detected at tolerance `tol`; grid may be irregular or `tol` too large.", call. = FALSE)
  }
  
  area_map <- stats::setNames(area_km2, grid_key)
  area_out <- unname(area_map[dxy_key])
  
  if (anyNA(area_out)) {
    warning(
      "Some lon/lat centers did not match the inferred grid; areas set to NA. ",
      "This can happen with floating precision. Consider rounding lon/lat or providing `cellsize=` explicitly.",
      call. = FALSE
    )
  }
  
  key <- data.frame(
    order = dxy$order,
    utm_x_scale = dxy$utm_x_scale,
    utm_y_scale = dxy$utm_y_scale,
    area_km2 = area_out
  )
  
  key <- key[order(key$order), , drop = FALSE]
  rownames(key) <- NULL
  key
}


#' Infer grid cell size (dx, dy) in degrees from lon/lat centers
#'
#' Internal helper used by `calculate_area(cellsize = "auto")`.
#' Tries to infer the most likely spacing from unique lon/lat values.
#'
#' @param lon,lat Numeric vectors of lon/lat centers.
#' @param tol Tolerance used to ignore tiny differences.
#' @param force `"auto"`, `"mode"`, or `"median"`. `"auto"` prefers mode then falls back.
#'
#' @return Numeric length-2 vector c(dx, dy).
#' @keywords internal
infer_cellsize_deg <- function(lon, lat, tol = 1e-6, force = c("auto", "mode", "median")) {
  force <- match.arg(force)
  
  infer_1d <- function(x, method) {
    x <- sort(unique(x))
    if (length(x) < 2) return(NA_real_)
    d <- diff(x)
    d <- d[is.finite(d) & d > tol]
    if (!length(d)) return(NA_real_)
    
    # stabilize floating noise
    md <- max(d)
    if (!is.finite(md) || md <= 0) return(NA_real_)
    nd <- max(0L, 8L - floor(log10(md)))
    d <- round(d, nd)
    
    if (method == "median") return(stats::median(d))
    
    # mode, prefer smallest among ties
    tab <- sort(table(d), decreasing = TRUE)
    top_vals <- as.numeric(names(tab)[tab == tab[1]])
    min(top_vals)
  }
  
  if (force == "mode") {
    dx <- infer_1d(lon, "mode");  dy <- infer_1d(lat, "mode")
  } else if (force == "median") {
    dx <- infer_1d(lon, "median"); dy <- infer_1d(lat, "median")
  } else {
    dx <- infer_1d(lon, "mode");  dy <- infer_1d(lat, "mode")
    if (!is.finite(dx) || !is.finite(dy)) {
      dx <- infer_1d(lon, "median"); dy <- infer_1d(lat, "median")
    }
  }
  
  if (!is.finite(dx) || !is.finite(dy)) {
    stop("Cannot infer cell size: not enough unique lon/lat values or diffs are all <= tol.", call. = FALSE)
  }
  c(dx, dy)
}


# Internal: lightweight sanity check that dx/dy match neighbor spacings
# If return_ok=TRUE, returns TRUE/FALSE instead of warning.
.check_grid_spacing <- function(lon, lat, dx, dy, tol = 1e-6, n = 200, seed = NULL, return_ok = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  if (length(unique(lon)) < 2 || length(unique(lat)) < 2) {
    if (return_ok) return(TRUE)
    return(invisible(TRUE))
  }
  
  idx <- sample.int(length(lon), size = min(n, length(lon)))
  
  ok_count <- 0L
  checked <- 0L
  
  for (ii in idx) {
    x <- lon[ii]
    y <- lat[ii]
    
    same_lat <- lon[abs(lat - y) <= tol]
    same_lon <- lat[abs(lon - x) <= tol]
    
    has_x <- FALSE
    has_y <- FALSE
    
    if (length(same_lat) >= 2) {
      dxx <- abs(sort(unique(same_lat)) - x)
      dxx <- dxx[dxx > tol]
      if (length(dxx)) has_x <- abs(min(dxx) - dx) <= max(tol, dx * 1e-3)
    }
    
    if (length(same_lon) >= 2) {
      dyy <- abs(sort(unique(same_lon)) - y)
      dyy <- dyy[dyy > tol]
      if (length(dyy)) has_y <- abs(min(dyy) - dy) <= max(tol, dy * 1e-3)
    }
    
    if (length(same_lat) >= 2 || length(same_lon) >= 2) {
      checked <- checked + 1L
      if (has_x || has_y) ok_count <- ok_count + 1L
    }
  }
  
  ok <- TRUE
  if (checked > 20 && ok_count / checked < 0.7) ok <- FALSE
  
  if (return_ok) return(ok)
  
  if (!ok) {
    warning(
      sprintf(
        paste0(
          "Inferred cellsize (dx=%.6g, dy=%.6g) may be inconsistent with lon/lat centers ",
          "(ok %.0f%%). If your grid is irregular or multi-resolution, supply `cellsize=` explicitly, e.g. ",
          "`cellsize = c(%.6g, %.6g)`."
        ),
        dx, dy, 100 * ok_count / checked, dx, dy
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Internal: key + area scaling + A_gs
#'
#' Builds the extrapolation/key table (unique grid cells), computes cell areas,
#' scales area if requested, and returns the basis matrix A_gs for the key grid.
#'
#' @param combined_data Data after UTM scaling, containing lon/lat and utm_x_scale/utm_y_scale.
#' @param mesh_fm A fmesher mesh object (not an intCPUEmesh wrapper).
#' @param area_scale Numeric or "auto".
#' @param ... Passed to `calculate_area()` (e.g., cellsize, crs_equal_area, check_grid).
#'
#' @return list(key=..., area_scale_val=..., A_gs=...)
#' @keywords internal
.prep_key_area <- function(combined_data, mesh_fm, area_scale = "auto", ...) {
  
  lon_col <- if ("lon_std" %in% names(combined_data)) "lon_std" else "lon"
  
  dxy <- unique(combined_data[, c(lon_col, "lat", "utm_x_scale", "utm_y_scale")])
  names(dxy)[names(dxy) == lon_col] <- "lon"
  rownames(dxy) <- NULL
  dxy$order <- seq_len(nrow(dxy))
  
  key <- calculate_area(dxy, lon_name = "lon", lat_name = "lat", ...)
  
  if (!("area_km2" %in% names(key))) {
    stop("`calculate_area()` must return a column named `area_km2`.", call. = FALSE)
  }
  
  if (identical(area_scale, "auto")) {
    area_scale_val <- .auto_scale_pow10(key$area_km2)
  } else {
    area_scale_val <- as.numeric(area_scale)
    if (!is.finite(area_scale_val) || area_scale_val <= 0) {
      stop("`area_scale` must be a positive number or 'auto'.", call. = FALSE)
    }
  }
  key$area_km2_scaled <- key$area_km2 / area_scale_val
  
  A_gs <- fmesher::fm_basis(mesh_fm, loc = as.matrix(key[, c("utm_x_scale", "utm_y_scale")]))
  
  list(key = key, area_scale_val = area_scale_val, A_gs = A_gs)
}

# ---- internal: quantized numeric key for floating coordinates ----
# Maps numeric values onto an integer lattice with resolution `tol`,
# allowing robust equality matching for floating-point coordinates.
# Two values within ±tol/2 will produce the same key.
#
# This is used to match grid cell centers after geometric operations
# (e.g., sf projection / polygon construction) where tiny floating
# differences would otherwise break exact matching (e.g., 0.3 vs
# 0.30000000000000004). It does NOT modify coordinates used for
# modeling — only the lookup key.
#
# The function assumes `tol` is much smaller than the grid spacing.
.qkey <- function(x, tol) {
  if (!is.numeric(x)) stop("Internal error: non-numeric coordinate.")
  if (!is.finite(tol) || tol <= 0) stop("Internal error: invalid `tol`.")
  as.integer(round(x / tol))
}