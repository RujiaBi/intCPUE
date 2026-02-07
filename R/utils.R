# Internal utilities
.check_required_cols <- function(data, cols) {
  miss <- setdiff(cols, names(data))
  if (length(miss)) {
    stop("Missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.check_numeric <- function(data, cols) {
  bad <- cols[!vapply(data[cols], is.numeric, logical(1))]
  if (length(bad)) {
    stop("These columns must be numeric: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# Scale utm_x/y
.auto_scale_pow10_xy <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  m <- max(abs(x), abs(y), na.rm = TRUE)
  if (!is.finite(m) || m <= 0) return(1)
  
  # choose scale = 10^k so that m/scale is roughly 1..100
  k <- floor(log10(m)) - 1
  scale <- 10^k
  
  # guard
  if (!is.finite(scale) || scale <= 0) scale <- 1
  scale
}

# Scale area of grid
.auto_scale_pow10 <- function(z) {
  z <- z[is.finite(z)]
  m <- max(abs(z))
  if (!is.finite(m) || m <= 0) return(1)
  k <- floor(log10(m))
  10^k
}