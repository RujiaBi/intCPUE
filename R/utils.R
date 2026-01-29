# Internal utilities

.check_required_cols <- function(df, req) {
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) {
    stop("`data` is missing required columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.check_numeric <- function(df, cols) {
  for (nm in cols) {
    if (!nm %in% names(df)) next
    if (!is.numeric(df[[nm]])) {
      stop("`", nm, "` must be numeric.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

.auto_scale_pow10_xy <- function(x, y) {
  z <- c(x, y)
  z <- z[is.finite(z)]
  m <- max(abs(z))
  if (!is.finite(m) || m <= 0) return(1)
  k <- floor(log10(m))
  10^k
}

.auto_scale_pow10 <- function(z) {
  z <- z[is.finite(z)]
  m <- max(abs(z))
  if (!is.finite(m) || m <= 0) return(1)
  k <- floor(log10(m))
  10^k
}