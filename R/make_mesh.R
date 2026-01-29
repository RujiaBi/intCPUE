# Define a non-convex boundary for INLA/fmesher mesh construction
make_boundary <- function(x, y, concave, convex, resx, resy) {
  border <- INLA::inla.nonconvex.hull(
    cbind(x, y),
    concave = concave,
    convex = convex,
    resolution = c(resx, resy)
  )
  
  border_coords <- border$loc[, c(1, 2)]
  border_polygon <- sp::Polygon(border_coords)
  border_polygons <- sp::Polygons(list(border_polygon), "border")
  border_spatial <- sp::SpatialPolygons(list(border_polygons))
  
  boundary_sf <- sf::st_as_sf(border_spatial)
  boundary_sp <- methods::as(boundary_sf, "Spatial")
  INLA::inla.sp2segment(boundary_sp)
}


# Create a tailored 2D INLA mesh
make_mesh <- function(
    utm_x, utm_y,
    concave = concave, convex = convex,
    resx = resx, resy = resy,
    max_edge = max_edge, bound_outer = bound_outer,
    inner_mesh_length = inner_mesh_length, outer_mesh_length = outer_mesh_length,
    cutoff_set = cutoff_set
) {
  boundary_0 <- make_boundary(
    x = utm_x, y = utm_y,
    concave = concave, convex = convex,
    resx = resx, resy = resy
  )
  
  fmesher::fm_mesh_2d_inla(
    boundary = boundary_0,
    max.edge = c(inner_mesh_length, outer_mesh_length) * max_edge,
    offset = c(max_edge, bound_outer),
    cutoff = max_edge / cutoff_set
  )
}


# Internal: UTM transform + shared scaling
.prep_utm_scaled <- function(data_input, utm_zone = NULL, coord_scale = "auto") {
  zones_each <- floor((data_input$lon + 180) / 6) + 1
  if (is.null(utm_zone)) {
    tab <- table(zones_each)
    utm_zone <- as.integer(names(tab)[which.max(tab)])
  } else {
    utm_zone <- as.integer(utm_zone)
  }
  
  mean_lat <- mean(data_input$lat, na.rm = TRUE)
  epsg_base <- if (mean_lat >= 0) 32600L else 32700L
  utm_epsg  <- epsg_base + utm_zone
  
  dd_sf  <- sf::st_as_sf(data_input, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  dd_utm <- sf::st_transform(dd_sf, crs = paste0("EPSG:", utm_epsg))
  xy_utm <- sf::st_coordinates(dd_utm)
  
  utm_x <- xy_utm[, 1]
  utm_y <- xy_utm[, 2]
  
  if (identical(coord_scale, "auto")) {
    utm_scale <- .auto_scale_pow10_xy(utm_x, utm_y)
  } else {
    utm_scale <- as.numeric(coord_scale)
    if (!is.finite(utm_scale) || utm_scale <= 0) {
      stop("`coord_scale` must be a positive number or 'auto'.", call. = FALSE)
    }
  }
  
  combined_data <- data_input
  combined_data$utm_x <- utm_x
  combined_data$utm_y <- utm_y
  combined_data$utm_x_scale <- utm_x / utm_scale
  combined_data$utm_y_scale <- utm_y / utm_scale
  
  list(
    combined_data = combined_data,
    utm_scale = utm_scale,
    utm_epsg = utm_epsg,
    utm_zone = utm_zone
  )
}


# Internal: key + area scaling + A_gs
.prep_key_area <- function(combined_data, mesh, area_scale = "auto") {
  # Keep both scaled and unscaled UTM coords for area computation
  dxy <- unique(combined_data[, c("lon", "lat", "utm_x_scale", "utm_y_scale")])
  rownames(dxy) <- NULL
  dxy$order <- seq_len(nrow(dxy))
  
  key <- calculate_area(dxy)
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
  
  A_gs <- fmesher::fm_basis(mesh, loc = as.matrix(key[, c("utm_x_scale", "utm_y_scale")]))
  
  list(key = key, area_scale_val = area_scale_val, A_gs = A_gs)
}

# Internal: anisotropy prep
.prep_anisotropy <- function(mesh, inla_spde) {
  Dset <- 1:2
  TV <- mesh$graph$tv
  V0 <- mesh$loc[TV[, 1], Dset]
  V1 <- mesh$loc[TV[, 2], Dset]
  V2 <- mesh$loc[TV[, 3], Dset]
  E0 <- V2 - V1
  E1 <- V0 - V2
  E2 <- V1 - V0
  
  tmp_det <- function(a, b) abs(det(rbind(a, b)))
  Tri_Area <- vapply(seq_len(nrow(E0)), function(i) tmp_det(E0[i, ], E1[i, ]) / 2, numeric(1))
  
  list(
    n_s      = mesh$n,
    n_tri    = nrow(TV),
    Tri_Area = Tri_Area,
    E0       = E0,
    E1       = E1,
    E2       = E2,
    TV       = TV - 1,  # 0-based for C++/TMB
    G0       = inla_spde$param.inla$M0,
    G0_inv   = methods::as(diag(1 / diag(inla_spde$param.inla$M0)), "TsparseMatrix")
  )
}