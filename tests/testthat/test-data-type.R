test_that("positive-only data ignores delta_link", {
  expect_identical(
    intCPUE:::.validate_delta_link_intCPUE("Poisson", data_type = "positive"),
    "traditional"
  )
  expect_identical(
    intCPUE:::.validate_delta_link_intCPUE("traditional", data_type = "positive"),
    "traditional"
  )
})

test_that("mixed data requires traditional delta link", {
  expect_error(
    intCPUE:::.validate_delta_link_intCPUE("Poisson", data_type = "mixed"),
    "requires `delta_link='traditional'`",
    fixed = TRUE
  )
  expect_identical(
    intCPUE:::.validate_delta_link_intCPUE("traditional", data_type = "mixed"),
    "traditional"
  )
})
