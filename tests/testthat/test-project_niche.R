test_that("project_ancestral_niche output_dir=NULL returns in-memory rasters", {
  skip_if_not_installed("terra")
  skip_if_not_installed("xnicher")
  env <- c(terra::rast(nrows = 10, ncols = 10, vals = stats::runif(100)),
           terra::rast(nrows = 10, ncols = 10, vals = stats::runif(100)))
  names(env) <- c("e1", "e2")
  anc_mu <- matrix(c(0, 0, 1, 1), 2, 2, byrow = TRUE)
  anc_Sigma <- array(0, c(2, 2, 2))
  anc_Sigma[1, , ] <- diag(2); anc_Sigma[2, , ] <- diag(2)
  tr <- ape::rtree(2, tip.label = c("sA", "sB"))
  tr$node.label <- c("n1")
  tmp <- tempfile()
  r <- project_ancestral_niche(phy = tr, anc_mu = anc_mu,
                               anc_Sigma = anc_Sigma, env = env,
                               node_labels = c("sA", "n1"),
                               output_dir = NULL)
  expect_length(r, 2)
  expect_s4_class(r[[1]], "SpatRaster")
  # in-memory: no file source attached
  expect_false(any(grepl(tmp, terra::sources(r[[1]]))))
  # bad node label caught
  expect_error(
    project_ancestral_niche(phy = tr, anc_mu = anc_mu,
                            anc_Sigma = anc_Sigma, env = env,
                            node_labels = c("sA", "WRONG")),
    "not found in 'phy'")
})
