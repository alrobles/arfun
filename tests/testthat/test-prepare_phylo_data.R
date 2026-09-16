test_that("prepare_phylo_data rejects unsupported dimensions", {
  skip_if_not_installed("ape")

  phy <- ape::rtree(3, tip.label = c("a", "b", "c"))
  niche <- data.frame(
    mu_1 = rep(0, 3), mu_2 = rep(0, 3), sigma_1 = rep(1, 3),
    sigma_2 = rep(1, 3), rho = rep(0, 3)
  )
  occ <- replicate(3, matrix(rnorm(6), ncol = 2), simplify = FALSE)
  env_m <- replicate(3, matrix(rnorm(6), ncol = 2), simplify = FALSE)

  expect_error(
    prepare_phylo_data(
      phy = phy,
      species = c("a", "b", "c"),
      niche_params = niche,
      occ_data = lapply(occ, function(x) cbind(x, extra = 0)),
      env_m_data = env_m
    ),
    "exactly 2 dimensions"
  )
})

test_that("prepare_phylo_data rejects malformed species-level inputs", {
  skip_if_not_installed("ape")

  phy <- ape::rtree(3, tip.label = c("a", "b", "c"))
  niche <- data.frame(
    mu_1 = rep(0, 3), mu_2 = rep(0, 3), sigma_1 = rep(1, 3),
    sigma_2 = rep(1, 3), rho = rep(0, 3)
  )
  occ <- replicate(3, matrix(rnorm(6), ncol = 2), simplify = FALSE)
  env_m <- replicate(3, matrix(rnorm(6), ncol = 2), simplify = FALSE)

  expect_error(
    prepare_phylo_data(
      phy = phy,
      species = c("a", "b", "c"),
      niche_params = niche,
      occ_data = occ,
      env_m_data = env_m[-1]
    ),
    "one entry per species"
  )
  expect_error(
    prepare_phylo_data(
      phy = phy,
      species = c("a", "b", "c"),
      niche_params = niche,
      occ_data = occ,
      env_m_data = env_m,
      mu_anc_prior = c(0, NA_real_)
    ),
    "finite"
  )
})

test_that("prepare_phylo_data aligns inputs and prunes extra tips", {
  skip_if_not_installed("ape")

  phy <- ape::rtree(4, tip.label = c("a", "b", "c", "extra"))
  species <- c("c", "a", "b")
  niche <- data.frame(
    mu_1 = c(30, 10, 20), mu_2 = c(3, 1, 2),
    sigma_1 = rep(1, 3), sigma_2 = rep(1, 3), rho = rep(0, 3)
  )
  occ <- list(
    matrix(c(30, 300), ncol = 2),
    matrix(c(10, 100), ncol = 2),
    matrix(c(20, 200), ncol = 2)
  )
  env_m <- list(
    matrix(c(31, 301), ncol = 2),
    matrix(c(11, 101), ncol = 2),
    matrix(c(21, 201), ncol = 2)
  )

  result <- prepare_phylo_data(
    phy = phy,
    species = species,
    niche_params = niche,
    occ_data = occ,
    env_m_data = env_m
  )

  expect_equal(result$S, 3L)
  expect_equal(dim(result$C), c(3L, 3L))
  expect_equal(rownames(result$C), c("a", "b", "c"))
  expect_equal(colnames(result$C), c("a", "b", "c"))
  expect_equal(result$occ[, 1], c(10, 20, 30))
  expect_equal(result$env_m[, 1], c(11, 21, 31))
  expect_equal(result$N_occ, c(1L, 1L, 1L))
  expect_equal(result$N_m, c(1L, 1L, 1L))
})
