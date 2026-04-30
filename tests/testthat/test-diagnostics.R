make_small_intcpue_input <- function() {
  data.frame(
    cpue = c(1.2, 0.0, 0.8, 1.1, 0.0, 0.7, 0.9, 0.0, 1.3, 0.6, 0.0, 1.0),
    encounter = c(1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L),
    lon = c(150.0, 150.1, 150.2, 150.3, 150.4, 150.5, 150.6, 150.7, 150.8, 150.9, 151.0, 151.1),
    lat = c(40.0, 40.0, 40.1, 40.1, 40.2, 40.2, 40.3, 40.3, 40.4, 40.4, 40.5, 40.5),
    vesid = rep(0:1, 6),
    tid = rep(0:2, each = 4),
    flagid = rep(c(0L, 1L), 6),
    depth = c(100, 110, 120, NA, 140, 150, 160, 170, 180, NA, 200, 210)
  )
}

test_that("smooth design matrices preserve row count with missing covariates", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("mgcv")

  data_input <- make_small_intcpue_input()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm

  mesh <- make_mesh(
    data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  prep <- make_data(
    formula_catchability = ~ s(depth),
    data_utm = data_utm,
    mesh = mesh
  )

  expect_equal(prep$data$has_smooths_catch, 1L)
  expect_equal(nrow(prep$data$Xs_catch), nrow(data_utm))
  expect_true(length(prep$data$Zs_catch) >= 1L)
  expect_true(all(vapply(prep$data$Zs_catch, nrow, integer(1)) == nrow(data_utm)))
  expect_equal(prep$data$has_smooths_pop, 0L)
})

test_that("population smooths are parsed separately and projected to extrapolation grid", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("mgcv")

  data_input <- make_small_intcpue_input()
  data_input$temp <- seq(10, 21, length.out = nrow(data_input))
  data_input$chl <- seq(0.5, 1.6, length.out = nrow(data_input))
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm

  mesh <- make_mesh(
    data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  prep <- make_data(
    formula_catchability = ~ s(depth),
    formula_population = ~ s(temp) + s(chl),
    data_utm = data_utm,
    mesh = mesh
  )

  expect_equal(prep$data$has_smooths_catch, 1L)
  expect_equal(prep$data$has_smooths_pop, 1L)
  expect_equal(nrow(prep$data$Xs_catch), nrow(data_utm))
  expect_equal(nrow(prep$data$Xs_pop_i), nrow(data_utm))
  expect_equal(nrow(prep$data$Xs_pop_g), prep$data$n_g * prep$data$n_t)
  expect_equal(length(prep$data$Zs_pop_i), length(prep$data$Zs_pop_g))
  expect_true(length(prep$smooth_basis$population) >= 1L)
  expect_equal(nrow(prep$projection_data), prep$data$n_g * prep$data$n_t)
})

test_that("population smooths require projection data when covariates vary within cell", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("mgcv")

  data_input <- make_small_intcpue_input()
  data_input$temp <- seq(10, 21, length.out = nrow(data_input))
  data_input$chl <- seq(0.5, 1.6, length.out = nrow(data_input))
  data_input$lon[2] <- data_input$lon[1]
  data_input$lat[2] <- data_input$lat[1]
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm

  mesh <- make_mesh(
    data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  expect_error(
    make_data(
      formula_population = ~ s(temp),
      data_utm = data_utm,
      mesh = mesh
    ),
    "varies within an extrapolation grid cell"
  )

  projection_data <- unique(data_utm[, c("utm_x_scale", "utm_y_scale"), drop = FALSE])
  projection_data$temp <- seq(50, 50 + nrow(projection_data) - 1)

  prep <- make_data(
    formula_population = ~ s(temp),
    data_utm = data_utm,
    mesh = mesh,
    projection_data = projection_data
  )

  expect_equal(prep$data$has_smooths_pop, 1L)
  expect_equal(nrow(prep$data$Xs_pop_g), nrow(prep$key) * prep$data$n_t)
})

test_that("population smooths accept time-varying projection data with tid", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("mgcv")

  data_input <- make_small_intcpue_input()
  data_input$temp <- seq(10, 21, length.out = nrow(data_input))
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm

  mesh <- make_mesh(
    data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  prep0 <- make_data(
    data_utm = data_utm,
    mesh = mesh
  )
  key <- prep0$key
  tid_values <- sort(unique(data_utm$tid))
  projection_data <- key[rep(seq_len(nrow(key)), times = length(tid_values)), c("utm_x_scale", "utm_y_scale"), drop = FALSE]
  projection_data$tid <- rep(tid_values, each = nrow(key))
  projection_data$temp <- projection_data$tid + seq_len(nrow(projection_data)) / 100

  prep <- make_data(
    formula_population = ~ s(temp),
    data_utm = data_utm,
    mesh = mesh,
    projection_data = projection_data
  )

  expect_equal(nrow(prep$projection_data), prep$data$n_g * prep$data$n_t)
  expect_equal(prep$projection_data$tid, rep(seq.int(0L, prep$data$n_t - 1L), each = prep$data$n_g))
  expect_equal(nrow(prep$data$Xs_pop_g), prep$data$n_g * prep$data$n_t)
})

test_that("projection_data NA values trigger a warning but still parse", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("mgcv")

  data_input <- make_small_intcpue_input()
  data_input$temp <- seq(10, 21, length.out = nrow(data_input))
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm

  mesh <- make_mesh(
    data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  prep0 <- make_data(
    data_utm = data_utm,
    mesh = mesh
  )

  projection_data <- prep0$key[, c("utm_x_scale", "utm_y_scale"), drop = FALSE]
  projection_data$temp <- seq(50, 50 + nrow(projection_data) - 1)
  projection_data$temp[1] <- NA_real_

  expect_warning(
    prep <- make_data(
      formula_population = ~ s(temp),
      data_utm = data_utm,
      mesh = mesh,
      projection_data = projection_data
    ),
    "contains NA in population covariates"
  )

  expect_equal(nrow(prep$data$Xs_pop_g), prep$data$n_g * prep$data$n_t)
  expect_true(anyNA(prep$projection_data$temp))
})

test_that("automatic population projection data preserves non-numeric covariate types", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("mgcv")

  data_input <- make_small_intcpue_input()
  data_input$temp <- seq(10, 21, length.out = nrow(data_input))
  data_input$season <- factor(rep(c("spring", "summer"), length.out = nrow(data_input)))
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm

  mesh <- make_mesh(
    data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  prep <- make_data(
    formula_population = ~ s(temp, by = season),
    data_utm = data_utm,
    mesh = mesh
  )

  expect_true(is.factor(prep$projection_data$season))
  expect_identical(levels(prep$projection_data$season), levels(data_utm$season))
  expect_equal(nrow(prep$data$Xs_pop_g), prep$data$n_g * prep$data$n_t)
})

test_that("diagnostic helpers return expected objects for a fitted model", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("TMB")
  skip_if_not_installed("ggplot2")

  data_input <- make_small_intcpue_input()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm

  mesh <- make_mesh(
    data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  fit <- suppressWarnings(
    intCPUE(
      formula_catchability = ~ 1,
      data_utm = data_utm,
      mesh = mesh,
      q_diffs_time = "off",
      q_diffs_spatial = "off",
      control = list(eval.max = 200, iter.max = 200),
      silent = TRUE
    )
  )

  conv <- check_convergence(fit)
  expect_true(all(c("convergence", "message", "max_grad") %in% names(conv)))

  aic_df <- calc_marginal_aic(fit)
  expect_true(all(c("logLik_marginal", "k", "AIC_marginal") %in% names(aic_df)))

  pred <- get_predicted(fit, data = data_utm)
  expect_true(all(c(
    "encounter_prob", "positive_mean", "fitted",
    "residual_raw", "residual_log"
  ) %in% names(pred)))
  expect_equal(sum(is.na(pred$residual_log)), sum(pred$observed <= 0))

  plots <- plot_residuals(pred)
  expect_s3_class(plots$observed_predicted, "ggplot")
  expect_s3_class(plots$spatial_residual, "ggplot")

  idx <- get_index(fit)
  p_idx <- plot_index(idx, time_values = c(2001, 2002, 2003))
  expect_s3_class(p_idx, "ggplot")

  p_aniso <- plot_anisotropy(fit)
  expect_s3_class(p_aniso, "ggplot")

  fitted_alias <- get_fitted(fit)
  expect_equal(fitted_alias$fitted, pred$fitted)
})
