#' Tests for prepare_niche_data
#'
#' Validates that the data preparation function produces correctly structured
#' Stan data lists for the presence-only model.
test_that("prepare_niche_data works for presence_only model", {
  set.seed(42)
  occ <- data.frame(
    bio01 = rnorm(30, 20, 3),
    bio12 = rnorm(30, 800, 100)
  )

  stan_data <- prepare_niche_data(occ)

  expect_equal(stan_data$N_occ, 30)
  expect_equal(stan_data$P, 2)
  expect_equal(dim(stan_data$occ), c(30, 2))
  expect_equal(length(stan_data$mu_prior), 2)
  expect_equal(length(stan_data$sigma_prior), 2)
  expect_equal(stan_data$eta, 1)
  expect_equal(stan_data$sigma_rate, 0.1)

  # Check that priors are data-adaptive
  expect_equal(stan_data$mu_prior, colMeans(occ))
  expect_equal(stan_data$sigma_prior,
               pmax(apply(occ, 2, stats::sd) * 10, 1e-6))
})

test_that("prepare_niche_data validates minimum occurrence points", {
  occ1 <- data.frame(bio01 = rnorm(1, 20, 3), bio12 = rnorm(1, 800, 100))

  expect_error(prepare_niche_data(occ1),
               "Need at least 2 occurrence points")
})

test_that("prepare_niche_data validates input dimensions", {
  occ1d <- data.frame(bio01 = rnorm(30, 20, 3))

  # 1D should work
  stan_data_1d <- prepare_niche_data(occ1d)
  expect_equal(stan_data_1d$P, 1)

  # Should work with matrix input
  mat <- matrix(rnorm(60), ncol = 2)
  stan_data_mat <- prepare_niche_data(mat)
  expect_equal(stan_data_mat$P, 2)
})

test_that("prepare_niche_data rejects invalid eta", {
  occ <- data.frame(bio01 = rnorm(30, 20, 3), bio12 = rnorm(30, 800, 100))

  expect_error(prepare_niche_data(occ, eta = "invalid"),
               "eta must be a positive finite numeric scalar")
  expect_error(prepare_niche_data(occ, eta = 0),
               "eta must be a positive finite numeric scalar")
  expect_error(prepare_niche_data(occ, eta = c(1, 2)),
               "eta must be a positive finite numeric scalar")
})

test_that("prepare_niche_data rejects invalid sigma_rate", {
  occ <- data.frame(bio01 = rnorm(30, 20, 3), bio12 = rnorm(30, 800, 100))

  expect_error(prepare_niche_data(occ, sigma_rate = -1),
               "sigma_rate must be a positive finite numeric scalar")
  expect_error(prepare_niche_data(occ, sigma_rate = 0),
               "sigma_rate must be a positive finite numeric scalar")
  expect_error(prepare_niche_data(occ, sigma_rate = NA_real_),
               "sigma_rate must be a positive finite numeric scalar")
})

test_that("prepare_niche_data accepts custom eta and sigma_rate", {
  occ <- data.frame(bio01 = rnorm(30, 20, 3), bio12 = rnorm(30, 800, 100))

  stan_data <- prepare_niche_data(occ, eta = 2, sigma_rate = 0.5)

  expect_equal(stan_data$eta, 2)
  expect_equal(stan_data$sigma_rate, 0.5)
})

test_that("prepare_niche_data rejects missing and non-finite values", {
  expect_error(
    prepare_niche_data(data.frame(bio01 = c(20, 21, NA_real_),
                                  bio12 = c(800, 810, 820))),
    "env_occ must contain only finite, non-missing values"
  )
  expect_error(
    prepare_niche_data(data.frame(bio01 = c(20, 21),
                                  bio12 = c(800, Inf))),
    "env_occ must contain only finite, non-missing values"
  )
  expect_error(
    prepare_niche_data(data.frame(bio01 = c("20", "21"),
                                  bio12 = c(800, 810))),
    "env_occ must contain numeric environmental values"
  )
})
