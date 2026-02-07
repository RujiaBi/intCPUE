#' @keywords internal
.as_intCPUEmesh <- function(mesh, loc_xy = NULL, xy_cols = NULL, recompute_A = c("auto", "always", "never")) {
  recompute_A <- match.arg(recompute_A)
  
  if (!requireNamespace("fmesher", quietly = TRUE)) {
    stop("Package 'fmesher' is required.", call. = FALSE)
  }
  
  # Case 1: already intCPUEmesh
  if (inherits(mesh, "intCPUEmesh")) {
    m <- mesh$mesh
    if (is.null(m)) stop("`mesh` is an intCPUEmesh but `mesh$mesh` is NULL.", call. = FALSE)
    
    # Ensure spde exists
    spde <- mesh$spde
    if (is.null(spde)) spde <- fmesher::fm_fem(m)
    
    # Decide if A can be reused
    A <- mesh$A
    need_A <- FALSE
    
    if (is.null(loc_xy)) {
      # no loc requested -> keep existing A (if present)
      if (is.null(A)) need_A <- TRUE
    } else {
      # loc requested -> A must correspond to loc_xy
      if (is.null(A)) {
        need_A <- TRUE
      } else if (recompute_A == "always") {
        need_A <- TRUE
      } else if (recompute_A == "auto") {
        # Heuristic: if dimensions mismatch, recompute
        if (!identical(nrow(A), nrow(loc_xy))) need_A <- TRUE
      }
    }
    
    if (need_A) {
      if (is.null(loc_xy)) stop("Need to compute A but `loc_xy` was not provided.", call. = FALSE)
      A <- fmesher::fm_basis(m, loc = loc_xy)
    }
    
    out <- mesh
    out$mesh <- m
    out$spde <- spde
    out$A <- A
    if (!is.null(loc_xy)) out$loc_xy <- loc_xy
    if (!is.null(xy_cols)) out$xy_cols <- xy_cols
    return(out)
  }
  
  # Case 2: bare fmesher mesh
  m <- mesh
  if (is.null(m)) stop("`mesh` is NULL.", call. = FALSE)
  
  if (is.null(loc_xy)) stop("When `mesh` is not an intCPUEmesh, you must provide `loc_xy`.", call. = FALSE)
  
  spde <- fmesher::fm_fem(m)
  A <- fmesher::fm_basis(m, loc = loc_xy)
  
  structure(list(
    loc_xy = loc_xy,
    xy_cols = if (is.null(xy_cols)) c(NA_character_, NA_character_) else xy_cols,
    mesh = m,
    spde = spde,
    A = A,
    loc_centers = NULL
  ), class = "intCPUEmesh")
}
