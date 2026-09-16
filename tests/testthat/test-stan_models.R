#' Tests for Stan model structure
#'
#' Validates that the Stan model file has the required blocks and declares
#' the expected data, parameters, and generated quantities.
test_that("niche_presence_only Stan model has required blocks", {
  model_path <- system.file("stan", "niche_presence_only.stan",
                            package = "arfun")
  content <- paste(readLines(model_path), collapse = "\n")

  expect_true(grepl("data \\{", content),
              info = "missing data block")
  expect_true(grepl("parameters \\{", content),
              info = "missing parameters block")
  expect_true(grepl("model \\{", content),
              info = "missing model block")
  expect_true(grepl("generated quantities \\{", content),
              info = "missing generated quantities block")
})

test_that("niche_presence_only Stan model declares required data fields", {
  model_path <- system.file("stan", "niche_presence_only.stan",
                            package = "arfun")
  content <- paste(readLines(model_path), collapse = "\n")

  expect_true(grepl("int.*N_occ", content),
              info = "missing N_occ (occurrence count)")
  expect_true(grepl("int.*P", content),
              info = "missing P (dimensions)")
  expect_true(grepl("matrix.*N_occ.*P.*occ", content),
              info = "missing occ (occurrence data matrix)")
  expect_true(grepl("mu_prior", content),
              info = "missing mu_prior (prior mean for centroid)")
  expect_true(grepl("sigma_prior", content),
              info = "missing sigma_prior (prior SD for centroid)")
  expect_true(grepl("sigma_rate", content),
              info = "missing sigma_rate (exponential prior rate)")
  expect_true(grepl("eta", content),
              info = "missing eta (LKJ shape parameter)")
})

test_that("niche_presence_only Stan model declares required parameters", {
  model_path <- system.file("stan", "niche_presence_only.stan",
                            package = "arfun")
  content <- paste(readLines(model_path), collapse = "\n")

  expect_true(grepl("vector\\[P\\] mu", content),
              info = "missing mu parameter (centroid)")
  expect_true(grepl("vector<lower=0>\\[P\\] sigma", content),
              info = "missing sigma parameter (standard deviations)")
  expect_true(grepl("cholesky_factor_corr\\[P\\] L_corr", content),
              info = "missing L_corr (Cholesky of correlation)")
})

test_that("niche_presence_only Stan model uses LKJ prior", {
  model_path <- system.file("stan", "niche_presence_only.stan",
                            package = "arfun")
  content <- paste(readLines(model_path), collapse = "\n")

  expect_true(grepl("lkj_corr_cholesky", content),
              info = "missing LKJ prior on correlation")
})

test_that("niche_presence_only Stan model uses multivariate normal likelihood", {
  model_path <- system.file("stan", "niche_presence_only.stan",
                            package = "arfun")
  content <- paste(readLines(model_path), collapse = "\n")

  expect_true(grepl("multi_normal_cholesky_lpdf", content),
              info = "missing multivariate normal likelihood")
})

test_that("niche_presence_only Stan model generates Sigma and R_corr", {
  model_path <- system.file("stan", "niche_presence_only.stan",
                            package = "arfun")
  content <- paste(readLines(model_path), collapse = "\n")

  expect_true(grepl("Sigma", content),
              info = "missing Sigma in generated quantities")
  expect_true(grepl("R_corr", content),
              info = "missing R_corr in generated quantities")
})

test_that("niche_presence_only Stan model does not require background data", {
  model_path <- system.file("stan", "niche_presence_only.stan",
                            package = "arfun")
  content <- paste(readLines(model_path), collapse = "\n")

  # Should NOT have N_m or env_m for presence_only
  expect_false(grepl("int.*N_m", content),
               info = "presence_only should not have N_m")
  expect_false(grepl("env_m", content),
               info = "presence_only should not have env_m")
})

test_that("stan model file exists", {
  model_path <- system.file("stan", "niche_presence_only.stan",
                            package = "arfun")
  expect_true(file.exists(model_path),
              info = "niche_presence_only.stan not found")
})
