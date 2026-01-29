# Calculate area for each lon/lat grid cell

calculate_area <- function(
    dxy,
    cellsize = "auto",
    lon_name = "X",
    lat_name = "Y",
    crs_equal_area = "+proj=eqearth +datum=WGS84",
    infer_method = c("mode", "median"),
    tol = 1e-6,
    check_grid = TRUE,
    check_n = 200,
    check_seed = NULL
) {
  infer_method <- match.arg(infer_method)
  dxy <- as.data.frame(dxy)
  
  # ---- accept first two columns as lon/lat if lon_name/lat_name not present ----
  if (!(lon_name %in% names(dxy)) || !(lat_name %in% names(dxy))) {
    if (ncol(dxy) < 2) stop("`dxy` must have lon/lat columns (or at least 2 columns).", call. = FALSE)
    names(dxy)[1:2] <- c(lon_name, lat_name)
  }
  
  # ---- required UTM columns ----
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
    cs <- infer_cellsize_deg(lon = lon, lat = lat, tol = tol, method = infer_method)
    dx <- cs[1]; dy <- cs[2]
  } else {
    if (length(cellsize) == 1L) {
      dx <- dy <- as.numeric(cellsize)
    } else if (length(cellsize) == 2L) {
      dx <- as.numeric(cellsize[1])
      dy <- as.numeric(cellsize[2])
    } else {
      stop('`cellsize` must be numeric length 1, 2, or "auto".', call. = FALSE)
    }
  }
  
  if (!is.finite(dx) || !is.finite(dy) || dx <= 0 || dy <= 0) {
    stop("`cellsize` must be positive and finite.", call. = FALSE)
  }
  
  # ---- optional sanity check of inferred grid spacing ----
  if (check_grid) {
    .check_grid_spacing(lon = lon, lat = lat, dx = dx, dy = dy, tol = tol, n = check_n, seed = check_seed)
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
  
  # ---- map area back to dxy without merge ----
  grid_key <- interaction(grid_sf$X, grid_sf$Y, drop = TRUE)
  dxy_key  <- interaction(dxy[[lon_name]], dxy[[lat_name]], drop = TRUE)
  
  area_map <- setNames(area_km2, as.character(grid_key))
  area_out <- unname(area_map[as.character(dxy_key)])
  
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


# Infer grid cell size (dx, dy) in degrees from lon/lat centers
infer_cellsize_deg <- function(lon, lat, tol = 1e-6, method = c("mode", "median")) {
  method <- match.arg(method)
  
  infer_1d <- function(x) {
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
    
    if (method == "median") {
      return(stats::median(d))
    }
    
    # mode, but prefer smallest common spacing if multiples exist
    tab <- sort(table(d), decreasing = TRUE)
    top_vals <- as.numeric(names(tab)[tab == tab[1]])
    return(min(top_vals))
  }
  
  dx <- infer_1d(lon)
  dy <- infer_1d(lat)
  
  if (!is.finite(dx) || !is.finite(dy)) {
    stop("Cannot infer cell size: not enough unique lon/lat values or diffs are all <= tol.", call. = FALSE)
  }
  c(dx, dy)
}

# internal: lightweight sanity check that dx/dy match neighbor spacings
.check_grid_spacing <- function(lon, lat, dx, dy, tol = 1e-6, n = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (length(unique(lon)) < 2 || length(unique(lat)) < 2) return(invisible(TRUE))
  
  idx <- sample.int(length(lon), size = min(n, length(lon)))
  
  ok_count <- 0L
  checked <- 0L
  
  for (ii in idx) {
    x <- lon[ii]
    y <- lat[ii]
    
    same_lat <- lon[lat == y]
    same_lon <- lat[lon == x]
    
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
  
  if (checked > 20 && ok_count / checked < 0.7) {
    warning(
      sprintf(
        "Inferred cellsize (dx=%.6g, dy=%.6g) may be inconsistent with lon/lat centers (ok %.0f%%). Supply `cellsize=` explicitly if needed.",
        dx, dy, 100 * ok_count / checked
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}
