ou_vcv_reference <- function(tanc, alpha, sigma_sq) {
  t_root <- diag(tanc)
  dist <- outer(t_root, t_root, "+") - 2 * tanc
  varanc <- -0.5 * expm1(-2 * alpha * tanc) / alpha
  sigma_sq * exp(-alpha * dist) * varanc
}

test_that("OU covariance has the BM limit for small alpha", {
  skip_if_not_installed("ape")

  set.seed(20260826)
  tree <- ape::rtree(6)
  tanc <- ape::vcv(tree)
  sigma_sq <- 0.7
  alpha <- 1e-6

  ou_cov <- ou_vcv_reference(tanc, alpha, sigma_sq)
  bm_limit <- sigma_sq * tanc

  expect_true(all(is.finite(ou_cov)))
  expect_equal(ou_cov, t(ou_cov), tolerance = 1e-12)
  expect_equal(ou_cov, bm_limit, tolerance = 1e-5)
})

test_that("OU covariance is positive definite across supported parameter scales", {
  skip_if_not_installed("ape")

  set.seed(20260827)
  tree <- ape::rtree(8)
  tanc <- ape::vcv(tree)
  parameter_grid <- expand.grid(
    alpha = c(1e-6, 1e-3, 0.1, 1, 10),
    sigma_sq = c(1e-8, 1e-4, 1, 10)
  )

  for (i in seq_len(nrow(parameter_grid))) {
    ou_cov <- ou_vcv_reference(
      tanc,
      parameter_grid$alpha[[i]],
      parameter_grid$sigma_sq[[i]]
    )
    eigenvalues <- eigen((ou_cov + t(ou_cov)) / 2,
                         symmetric = TRUE,
                         only.values = TRUE)$values

    expect_true(all(is.finite(ou_cov)))
    expect_true(min(eigenvalues) > 0)
  }
})

test_that("null model keeps species parameters independent of phylogeny", {
  path <- model_path("null")
  content <- paste(readLines(path), collapse = "\n")

  expect_false(grepl("T_anc|matrix\\[S,\\s*S\\]\\s+C\\b", content))
  expect_match(content, "matrix\\[S, P\\] mu;")
  expect_match(content, "matrix\\[S, P\\] log_sigma;")
  expect_match(content, "vector\\[S\\] z_rho;")
  expect_match(content, "array\\[S\\] matrix\\[P, P\\] Sigma;")
  expect_match(content, "vector\\[S\\] log_lik;")
})
