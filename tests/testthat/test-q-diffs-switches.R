test_that("q-diff switches are honored in TMB when turned off", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("TMB")

  data_input <- data.frame(
    cpue = c(1.2, 1e-6, 0.8, 1.1, 1e-6, 0.7, 0.9, 1e-6, 1.3, 0.6, 1e-6, 1.0),
    encounter = c(1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L),
    lon = c(150.0, 150.1, 150.2, 150.3, 150.4, 150.5, 150.6, 150.7, 150.8, 150.9, 151.0, 151.1),
    lat = c(40.0, 40.0, 40.1, 40.1, 40.2, 40.2, 40.3, 40.3, 40.4, 40.4, 40.5, 40.5),
    vesid = rep(0:1, 6),
    tid = rep(0:2, each = 4),
    flagid = rep(c(0L, 1L), 6)
  )

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

  rep <- fit$obj$report()

  expect_equal(as.integer(rep$use_q_diffs_time), 0L)
  expect_equal(as.integer(rep$q_diffs_spatial_enabled), 0L)
  expect_equal(as.integer(rep$use_pop_spatiotemporal_rw), 1L)
  expect_equal(as.integer(rep$use_pop_spatiotemporal_ar1), 0L)
  expect_s3_class(fit, "intCPUE")
})

test_that("flags with fewer than two time points are excluded from flag_t constraints", {
  has_tf <- matrix(
    c(
      TRUE,  FALSE,
      FALSE, TRUE,
      FALSE, TRUE,
      FALSE, FALSE
    ),
    nrow = 4,
    byrow = TRUE
  )

  constraint <- intCPUE:::.build_flag_t_constraint(has_tf)

  expect_equal(constraint$estimable_flag, c(FALSE, TRUE))
  expect_equal(constraint$n_free, 1L)
  expect_equal(
    constraint$flag_t_index,
    matrix(
      c(
        -1L, -1L,
        -1L, -1L,
        -1L,  0L,
        -1L, -1L
      ),
      nrow = 4,
      byrow = TRUE
    )
  )
})

test_that("observation-SD mapping is honored", {
  parameters <- intCPUE:::.make_parameters_intCPUE(
    n_t = 3L,
    n_v = 1L,
    n_f = 2L,
    n_s = 4L,
    K_smooth_catch = 0L,
    n_smooth_catch = 0L,
    sum_k_catch = 0L,
    K_smooth_pop = 0L,
    n_smooth_pop = 0L,
    sum_k_pop = 0L
  )

  map <- intCPUE:::.make_map_intCPUE(
    parameters = parameters,
    n_f = 2L,
    n_v = 1L,
    q_diffs_time = "off",
    obs_sd = "flag"
  )

  expect_true(is.na(map$ln_sd))
  expect_null(map$ln_sd_flag)

  shared_map <- intCPUE:::.make_map_intCPUE(
    parameters = parameters,
    n_f = 2L,
    n_v = 1L,
    q_diffs_time = "off",
    obs_sd = "shared"
  )

  expect_true(all(is.na(shared_map$ln_sd_flag)))
  expect_null(shared_map[["ln_sd", exact = TRUE]])
})

test_that("RW maps out unused AR1 correlation parameters", {
  parameters <- intCPUE:::.make_parameters_intCPUE(
    n_t = 3L,
    n_v = 1L,
    n_f = 2L,
    n_s = 4L,
    K_smooth_catch = 0L,
    n_smooth_catch = 0L,
    sum_k_catch = 0L,
    K_smooth_pop = 0L,
    n_smooth_pop = 0L,
    sum_k_pop = 0L
  )

  rw_map <- intCPUE:::.make_map_intCPUE(
    parameters = parameters,
    n_f = 2L,
    n_v = 1L,
    q_diffs_time = "off",
    obs_sd = "shared",
    pop_spatiotemporal_type = "rw"
  )

  expect_true(is.na(rw_map$transf_rho_1))
  expect_true(is.na(rw_map$transf_rho_2))

  ar1_map <- intCPUE:::.make_map_intCPUE(
    parameters = parameters,
    n_f = 2L,
    n_v = 1L,
    q_diffs_time = "off",
    obs_sd = "shared",
    pop_spatiotemporal_type = "ar1"
  )

  expect_null(ar1_map[["transf_rho_1", exact = TRUE]])
  expect_null(ar1_map[["transf_rho_2", exact = TRUE]])
})

test_that("observation-SD and AR1 spatiotemporal settings are reported by TMB", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("TMB")

  data_input <- data.frame(
    cpue = c(1.2, 1e-6, 0.8, 1.1, 1e-6, 0.7, 0.9, 1e-6, 1.3, 0.6, 1e-6, 1.0),
    encounter = c(1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L),
    lon = c(150.0, 150.1, 150.2, 150.3, 150.4, 150.5, 150.6, 150.7, 150.8, 150.9, 151.0, 151.1),
    lat = c(40.0, 40.0, 40.1, 40.1, 40.2, 40.2, 40.3, 40.3, 40.4, 40.4, 40.5, 40.5),
    vesid = rep(0:1, 6),
    tid = rep(0:2, each = 4),
    flagid = rep(c(0L, 1L), 6)
  )

  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  data_utm <- utm$data_utm

  mesh <- make_mesh(
    data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  fit_flag_sd <- suppressWarnings(
    intCPUE(
      formula_catchability = ~ 1,
      data_utm = data_utm,
      mesh = mesh,
      obs_sd = "flag",
      q_diffs_time = "off",
      q_diffs_spatial = "off",
      control = list(eval.max = 200, iter.max = 200),
      silent = TRUE
    )
  )

  rep_flag_sd <- fit_flag_sd$obj$report()
  expect_equal(as.integer(rep_flag_sd$use_pop_spatiotemporal_rw), 1L)
  expect_equal(as.integer(rep_flag_sd$use_pop_spatiotemporal_ar1), 0L)
  expect_equal(as.integer(rep_flag_sd$use_flag_sd), 1L)
  expect_length(rep_flag_sd$sd_flag, fit_flag_sd$data_tmb$n_f)

  fit_ar1 <- suppressWarnings(
    intCPUE(
      formula_catchability = ~ 1,
      data_utm = data_utm,
      mesh = mesh,
      pop_spatiotemporal_type = "ar1",
      obs_sd = "shared",
      q_diffs_time = "off",
      q_diffs_spatial = "off",
      control = list(eval.max = 200, iter.max = 200),
      silent = TRUE
    )
  )

  rep_ar1 <- fit_ar1$obj$report()
  expect_equal(as.integer(rep_ar1$use_pop_spatiotemporal_rw), 0L)
  expect_equal(as.integer(rep_ar1$use_pop_spatiotemporal_ar1), 1L)
  expect_equal(as.integer(rep_ar1$use_flag_sd), 0L)
})
