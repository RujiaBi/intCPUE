#' Convert lon/lat to UTM coordinates with robust longitude handling
#'
#' Convert geographic coordinates (lon/lat, WGS84) to UTM coordinates, with:
#' - robust handling of lon in either 0..360 or -180..180
#' - automatic UTM zone selection (more stable near the dateline)
#' - optional coordinate scaling (auto or user-specified)
#'
#' @param data_input data.frame with columns lon, lat (and others you need)
#' @param utm_zone integer in 1--60; if NULL, auto-selected
#' @param coord_scale positive numeric or "auto"
#' @param lon_convention "auto", "-180_180", or "0_360"
#'   - "auto": detect and internally standardize lon for projection/zone pick
#'   - "-180_180": treat input as -180..180 (but still normalize safely)
#'   - "0_360": treat input as 0..360 (but still normalize safely)
#' @param quiet logical; if FALSE, may message selected EPSG/zone
#'
#' @return A list with data_utm (adds utm_x/utm_y and scaled versions),
#'   plus utm_scale, utm_epsg, utm_zone.
#' @export
make_utm <- function(
    data_input,
    utm_zone = NULL,
    coord_scale = "auto",
    lon_convention = c("auto", "-180_180", "0_360"),
    quiet = TRUE
) {
  data_input <- as.data.frame(data_input)
  lon_convention <- match.arg(lon_convention)
  
  .check_required_cols(data_input, c("lon", "lat"))
  .check_numeric(data_input, c("lon", "lat"))
  
  if (anyNA(data_input$lon) || anyNA(data_input$lat)) {
    stop("`lon`/`lat` must not contain NA.", call. = FALSE)
  }
  
  out <- .prep_utm_scaled(
    data_input = data_input,
    utm_zone = utm_zone,
    coord_scale = coord_scale,
    lon_convention = lon_convention,
    quiet = quiet
  )
  
  out
}

# ---- internal: UTM transform + shared scaling ----
.prep_utm_scaled <- function(
    data_input,
    utm_zone = NULL,
    coord_scale = "auto",
    lon_convention = c("auto", "-180_180", "0_360"),
    quiet = TRUE
) {
  lon_convention <- match.arg(lon_convention)
  
  # 1) normalize lon based on declared convention (or auto)
  lon_in <- as.numeric(data_input$lon)
  lat_in <- as.numeric(data_input$lat)
  
  if (any(lat_in < -90 | lat_in > 90, na.rm = TRUE)) {
    stop("`lat` must be within [-90, 90].", call. = FALSE)
  }
  
  if (any(lat_in < -80 | lat_in > 84, na.rm = TRUE)) {
    stop("UTM is not recommended for lat outside [-80, 84].", call. = FALSE)
  }
  
  # If user says input is 0_360, normalize accordingly first (optional),
  # then convert to -180..180 for zone selection and projection.
  if (lon_convention == "0_360") {
    lon_tmp <- .normalize_lon(lon_in, convention = "0_360")
    lon_std <- .normalize_lon(lon_tmp, convention = "-180_180")
  } else if (lon_convention == "-180_180") {
    lon_std <- .normalize_lon(lon_in, convention = "-180_180")
  } else {
    # "auto": detect and normalize robustly to -180..180
    lon_std <- .normalize_lon(lon_in, convention = "auto")
  }
  
  # 2) choose UTM zone
  if (is.null(utm_zone)) {
    utm_zone <- .pick_utm_zone(lon_std)
  } else {
    utm_zone <- as.integer(utm_zone)
    if (is.na(utm_zone) || !is.finite(utm_zone) || utm_zone < 1L || utm_zone > 60L) {
      stop("`utm_zone` must be an integer in [1, 60].", call. = FALSE)
    }
  }
  
  # 3) diagnostics
  cross_equator <- any(lat_in < 0) && any(lat_in > 0)
  if (cross_equator && !quiet) {
    message("Data crosses the equator. Using a single UTM hemisphere may introduce small metric distortion near the equator.")
  }
  
  lon_span <- .lon_span_circular(lon_std)
  if (lon_span > 12 && !quiet) {
    message(sprintf(
      "Longitude span is %.1f degrees (> 12). A single UTM zone may distort distances; consider AEQD/LAEA for large domains.",
      lon_span
    ))
  }
  
  # 4) hemisphere -> EPSG
  n_north <- sum(lat_in >= 0)
  n_south <- sum(lat_in < 0)
  epsg_base <- if (n_north >= n_south) 32600L else 32700L
  utm_epsg <- epsg_base + utm_zone
  
  if (!quiet) {
    message(sprintf("Using UTM zone %d (EPSG:%d).", utm_zone, utm_epsg))
  }
  
  # 4) project with sf
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for UTM transforms.", call. = FALSE)
  }
  
  dd_sf <- sf::st_as_sf(
    transform(data_input, lon_std = lon_std, lat = lat_in),
    coords = c("lon_std", "lat"),
    crs = 4326,
    remove = FALSE
  )
  
  dd_utm <- sf::st_transform(dd_sf, crs = paste0("EPSG:", utm_epsg))
  xy_utm <- sf::st_coordinates(dd_utm)
  
  utm_x <- xy_utm[, 1]
  utm_y <- xy_utm[, 2]
  
  # 5) scaling
  if (identical(coord_scale, "auto")) {
    utm_scale <- .auto_scale_pow10_xy(utm_x, utm_y)
  } else {
    utm_scale <- as.numeric(coord_scale)
    if (is.na(utm_scale) || !is.finite(utm_scale) || utm_scale <= 0) {
      stop("`coord_scale` must be a positive number or 'auto'.", call. = FALSE)
    }
  }
  
  data_utm <- data_input
  data_utm$lon_std <- lon_std
  data_utm$utm_x <- utm_x
  data_utm$utm_y <- utm_y
  data_utm$utm_x_scale <- utm_x / utm_scale
  data_utm$utm_y_scale <- utm_y / utm_scale
  
  list(
    data_utm = data_utm,
    utm_scale = utm_scale,
    utm_epsg = utm_epsg,
    utm_zone = utm_zone
  )
}

# ---- internal: normalize longitude to a target convention ----
.normalize_lon <- function(lon, convention = c("auto", "-180_180", "0_360")) {
  convention <- match.arg(convention)
  lon <- as.numeric(lon)
  
  if (anyNA(lon)) stop("`lon` must not contain NA.", call. = FALSE)
  
  if (convention == "auto") {
    # We internally standardize to -180..180 for UTM zone logic
    convention <- "-180_180"
    # If input looks like 0..360, the transformation below will handle it anyway.
  }
  
  if (convention == "0_360") {
    lon2 <- lon %% 360
    lon2[lon2 == 360] <- 0
    return(lon2)
  }
  
  # "-180_180"
  lon2 <- ((lon + 180) %% 360) - 180
  lon2
}

# ---- internal: pick a UTM zone robustly (handles dateline better) ----
.pick_utm_zone <- function(lon_180) {
  lon_180 <- as.numeric(lon_180)
  
  # circular mean for dateline robustness
  rad <- lon_180 * pi / 180
  mean_ang <- atan2(mean(sin(rad), na.rm = TRUE), mean(cos(rad), na.rm = TRUE))
  mean_lon <- mean_ang * 180 / pi
  
  zone <- floor((mean_lon + 180) / 6) + 1
  zone <- max(1L, min(60L, as.integer(zone)))
  zone
}

# ---- internal: circular longitude span (dateline-robust) ----
# Computes the minimal longitudinal span on a circle for longitudes expressed
# in [-180, 180] degrees, avoiding artificial inflation when observations lie
# on both sides of the dateline (e.g. -179° and 179°).
# Unlike `diff(range(lon))`, this measures the shortest arc containing all points.
# The span is computed as 360° minus the largest gap between sorted longitudes
# on the circle.
# This is used only for diagnostics (e.g., warning about using a single UTM zone),
# not for projection itself.
.lon_span_circular <- function(lon_180) {
  lon_180 <- as.numeric(lon_180)
  lon_180 <- lon_180[is.finite(lon_180)]
  if (!length(lon_180)) return(NA_real_)
  
  x <- sort((lon_180 + 360) %% 360)   # 0..360
  gaps <- diff(c(x, x[1] + 360))
  360 - max(gaps)
}


