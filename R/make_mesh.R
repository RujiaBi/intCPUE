# ------------------------------------------------------------------------------
# Mesh construction workflow inspired by sdmTMB:
# https://github.com/sdmTMB/sdmTMB
#
# Anisotropy helper adapted from the TMB examples repository:
# https://github.com/kaskr/adcomp/tree/master/tmb_examples
# ------------------------------------------------------------------------------

#' Construct an SPDE mesh for intCPUE
#'
#' mesh constructor with convenience options:
#' - "cutoff": simplest mesh given `cutoff`
#' - "kmeans": place vertices using kmeans centers (`n_knots`)
#' - "tailored": build from all points but pass advanced args via `...`
#' - "custom": supply your own `mesh`
#'
#' @param data_utm A data.frame with utm coords.
#' @param xy_cols Character vector length 2 giving x/y column names in `data`.
#' @param type Mesh construction method.
#' @param cutoff Minimum allowed triangle edge length. Required for type="cutoff".
#' @param n_knots Number of knots if type="kmeans".
#' @param seed Seed for kmeans if type="kmeans".
#' @param mesh Optional custom mesh (required for type="custom").
#' @param fmesher_func Which fmesher function to use. Default is
#'   [fmesher::fm_rcdt_2d_inla()] (supports very simple cutoff meshes).
#' @param convex,concave If supplied, passed to [fmesher::fm_nonconvex_hull()]
#'   and used as a boundary.
#' @param ... Passed to `fmesher_func`. Examples: `max.edge`, `offset`,
#'   `cutoff`, `extend`, `refine`.
#'
#' @return An object of class `intCPUEmesh` with elements `mesh`, `spde`, `A`,
#'   `loc_xy`, `xy_cols`, and `loc_centers` (for kmeans).
#' @export
make_mesh <- function(
    data_utm,
    xy_cols,
    type = c("kmeans", "cutoff", "tailored", "custom"),
    cutoff,
    n_knots,
    seed = 42,
    mesh = NULL,
    fmesher_func = fmesher::fm_rcdt_2d_inla,
    convex = NULL,
    concave = convex,
    ...
) {
  type <- match.arg(type)
  
  if (!requireNamespace("fmesher", quietly = TRUE)) {
    stop("Package 'fmesher' is required.", call. = FALSE)
  }
  if (!is.data.frame(data_utm) || nrow(data_utm) == 0L) {
    stop("`data` must be a non-empty data.frame.", call. = FALSE)
  }
  if (missing(xy_cols) || length(xy_cols) != 2L) {
    stop("`xy_cols` must be a character vector of length 2.", call. = FALSE)
  }
  
  .check_required_cols(data_utm, xy_cols)
  .check_numeric(data_utm, xy_cols)
  
  all_x_non_na <- sum(is.na(data_utm[[xy_cols[[1]]]])) == 0L
  all_y_non_na <- sum(is.na(data_utm[[xy_cols[[2]]]])) == 0L
  if (!all_x_non_na || !all_y_non_na) {
    stop("Some coordinates in `xy_cols` were NA. Remove/fix these rows.", call. = FALSE)
  }
  
  loc_xy <- as.matrix(data_utm[, xy_cols, drop = FALSE])
  loc_centers <- NULL
  
  # boundary only if convex/concave provided
  nch <- NULL
  if (!is.null(convex) || !is.null(concave)) {
    nch <- fmesher::fm_nonconvex_hull(loc_xy, convex = convex, concave = concave)
  }
  
  # Decide type automatically: if cutoff provided and n_knots missing
  if (!missing(cutoff) && missing(n_knots) && type == "kmeans") {
    type <- "cutoff"
  }
  
  # argument checks
  if (type == "cutoff" && is.null(mesh) && missing(cutoff)) {
    stop("You need to specify `cutoff` when type = 'cutoff'.", call. = FALSE)
  }
  if (type == "kmeans" && is.null(mesh) && missing(n_knots)) {
    stop("You need to specify `n_knots` when type = 'kmeans'.", call. = FALSE)
  }
  if (type == "custom" && is.null(mesh)) {
    stop("You need to supply `mesh` when type = 'custom'.", call. = FALSE)
  }
  
  # Build mesh
  if (is.null(mesh)) {
    dots <- list(...)
    
    if (type == "kmeans") {
      if (n_knots >= nrow(loc_xy)) {
        warning("Reducing `n_knots` to be one less than the number of data points.", call. = FALSE)
        n_knots <- nrow(loc_xy) - 1
      }
      set.seed(seed)
      km <- stats::kmeans(x = loc_xy, centers = n_knots)
      loc_centers <- km$centers
      
      base_args <- list(
        loc_centers,
        refine = list(),
        extend = list()
      )
      if (!is.null(nch)) base_args$boundary <- nch
      
      mesh <- do.call(fmesher_func, c(base_args, dots))
      
    } else if (type == "cutoff") {
      base_args <- list(
        loc_xy,
        refine = list(),
        extend = list(),
        cutoff = cutoff
      )
      if (!is.null(nch)) base_args$boundary <- nch
      
      mesh <- do.call(fmesher_func, c(base_args, dots))
      
    } else if (type == "tailored") {
      # "tailored": same as cutoff style but expect users to pass max.edge/offset/etc via ...
      # We allow cutoff either from argument or via ...
      base_args <- list(
        loc_xy,
        refine = list(),
        extend = list()
      )
      if (!is.null(nch)) base_args$boundary <- nch
      if (!missing(cutoff)) base_args$cutoff <- cutoff
      
      mesh <- do.call(fmesher_func, c(base_args, dots))
    }
  }
  
  # FEM / SPDE + basis matrix
  spde <- fmesher::fm_fem(mesh)
  A <- fmesher::fm_basis(mesh, loc = loc_xy)
  
  out <- structure(list(
    loc_xy = loc_xy,
    xy_cols = xy_cols,
    mesh = mesh,
    spde = spde,
    A = A,
    loc_centers = loc_centers
  ), class = "intCPUEmesh")
  
  if (!is.null(mesh$n) && mesh$n > 1000L) {
    message("This mesh has > 1000 vertices; consider simplifying for speed.")
  }
  
  out
}

#' Plot an intCPUE mesh
#'
#' @param x An intCPUEmesh object
#' @param ... passed to plot()
#' @export
#' @method plot intCPUEmesh
plot.intCPUEmesh <- function(x, ...) {
  graphics::plot(x$mesh, main = NA, asp = 1, ...)
  graphics::points(x$loc_xy, pch = 21, cex = 0.3, col = "#00000080")
  if (!is.null(x$loc_centers) && all(is.finite(x$loc_centers))) {
    graphics::points(x$loc_centers, pch = 20, col = "red")
  }
  invisible(x)
}

# Internal: anisotropy (from TMB examples repository)
.prep_anisotropy <- function(mesh, spde) {
    Dset <- 1:2
    # Triangle info
    TV <- mesh$graph$tv # Triangle to vertex indexing
    V0 <- mesh$loc[TV[, 1], Dset] # V = vertices for each triangle
    V1 <- mesh$loc[TV[, 2], Dset]
    V2 <- mesh$loc[TV[, 3], Dset]
    E0 <- V2 - V1 # E = edge for each triangle
    E1 <- V0 - V2
    E2 <- V1 - V0
    # Calculate Areas
    TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
    Tri_Area <- rep(NA, nrow(E0))
    for (i in seq_len(length(Tri_Area))) Tri_Area[i] <- TmpFn(E0[i, ], E1[i, ]) / 2
    
    list(
      n_s = mesh$n, 
      n_tri = nrow(TV),
      Tri_Area = Tri_Area,
      E0 = E0,
      E1 = E1,
      E2 = E2,
      TV = TV - 1,
      G0 = spde$c0,
      G0_inv = methods::as(Matrix::diag(1 / Matrix::diag(spde$c0)), "TsparseMatrix")
    )
}