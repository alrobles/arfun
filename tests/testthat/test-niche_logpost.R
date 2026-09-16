#' Tests for the math-scale joint log-posterior and multistart machinery
#'
#' Covers R/C++ parity, the masked-species likelihood semantics, the
#' shared-background contract, boundary diagnostics, and deterministic
#' multistart behaviour.
#'
#' @noRd
library(testthat)
library(arfun)

skip_if_not_installed("ape")

make_tiny_clade <- function(S = 4, n_occ = 60, n_bg = 800, seed = 7) {
  set.seed(seed)
  tr <- ape::rtree(S, tip.label = paste0("sp", seq_len(S)))
  mu_t <- matrix(c(0, 0, 0.5, 0.2, -0.4, 0.3, 0.1, -0.5)[seq_len(2 * S)],
                 S, 2, byrow = TRUE)
  occ <- lapply(seq_len(S), function(s)
    matrix(rnorm(n_occ * 2, rep(mu_t[s, ], each = n_occ), 0.3), n_occ, 2))
  names(occ) <- tr$tip.label
  bg <- matrix(runif(n_bg * 2, -2, 2), n_bg, 2)
  bgl <- lapply(seq_len(S), function(s) bg)
  names(bgl) <- tr$tip.label
  list(occ = occ, bgl = bgl, tree = tr)
}

test_that("prepare_phylo_niche_data builds consistent blocks", {
  cl <- make_tiny_clade()
  d <- prepare_phylo_niche_data(cl$occ, cl$bgl, cl$tree)
  expect_equal(d$S, 4)
  expect_equal(d$P, 2)
  expect_equal(d$N_occ, rep(60L, 4))
  expect_equal(d$N_occ_total, 240L)
  expect_equal(dim(d$C), c(4, 4))
  expect_true(all(!d$masked))
  # starts tile the concatenated rows with no gaps
  expect_equal(d$occ_start, c(1L, 61L, 121L, 181L))
})

test_that("shared_background_sample returns identical rows for all species", {
  bg <- matrix(runif(2000), 1000, 2)
  a <- shared_background_sample(bg, n = 500, seed = 3)
  b <- shared_background_sample(bg, n = 500, seed = 3)
  expect_identical(a, b)
  expect_equal(nrow(a), 500)
  # same index vector for every species -> commensurate normalizer
  bgl <- lapply(1:3, function(s) a)
  expect_true(all(vapply(bgl, function(x) identical(x, bgl[[1]]), TRUE)))
})

test_that("mask_species drops likelihood rows but keeps the phylo prior", {
  cl <- make_tiny_clade()
  d0 <- prepare_phylo_niche_data(cl$occ, cl$bgl, cl$tree)
  dm <- prepare_phylo_niche_data(cl$occ, cl$bgl, cl$tree,
                                 mask_species = "sp2")
  i2 <- match("sp2", dm$species)
  expect_equal(dm$N_occ[i2], 0L)
  expect_equal(dm$N_m[i2], 0L)
  expect_true(dm$masked[i2])
  expect_equal(sum(dm$masked), 1L)
  expect_equal(dim(dm$C), c(4, 4))  # masked species stays in the tree
  th <- niche_start_ranges(d0)$center
  lp0 <- niche_logpost(th, d0, "noncentered_bounded")
  lpm <- niche_logpost(th, dm, "noncentered_bounded")
  expect_gt(lpm, lp0)  # removing a likelihood term raises the lp
  # masking a non-tip errors
  expect_error(prepare_phylo_niche_data(cl$occ, cl$bgl, cl$tree,
                                        mask_species = "nope"),
               "not in tree")
})

test_that("R and C++ objectives agree to machine precision", {
  cl <- make_tiny_clade()
  d <- prepare_phylo_niche_data(cl$occ, cl$bgl, cl$tree)
  rng <- niche_start_ranges(d, "noncentered_bounded")
  xp <- niche_logpost_xptr(d, "noncentered_bounded")
  for (param in c("noncentered_bounded", "noncentered",
                  "centered", "centered_lograte")) {
    rng <- niche_start_ranges(d, param)
    xp <- niche_logpost_xptr(d, param)
    probe <- rbind(rng$center,
                   rng$center + 0.1,
                   (rng$lo + rng$hi) / 2)
    for (i in seq_len(nrow(probe))) {
      lr <- niche_logpost(probe[i, ], d, param)
      lc <- niche_logpost_eval_xptr(xp, probe[i, ])
      if (is.finite(lr) && is.finite(lc)) {
        expect_lt(abs(lr - lc), 1e-8)
      } else {
        expect_identical(is.finite(lr), is.finite(lc))
      }
    }
  }
})

test_that("theta layout matches Stan column-major serialization", {
  nm <- niche_theta_layout(3, 2, "noncentered_bounded")
  # order: anc blocks | mu_raw[S,P] col-major | log_sigma_raw | rho block
  expect_equal(nm[1:4], c("mu_anc[1]", "mu_anc[2]",
                          "log_sigma_anc[1]", "log_sigma_anc[2]"))
  expect_true("mu_raw[3,1]" %in% nm)
  expect_true("log_sigma_raw[3,2]" %in% nm)
  expect_equal(tail(nm, 3), c("z_rho_raw[1]", "z_rho_raw[2]", "z_rho_raw[3]"))
  # length: 4P + 2*S*P + 2 + S = 8 + 12 + 2 + 3
  expect_length(nm, 4 * 2 + 2 * 3 * 2 + 2 + 3)
})

test_that("sigma floor rejects the degenerate point-mass mode", {
  cl <- make_tiny_clade()
  d <- prepare_phylo_niche_data(cl$occ, cl$bgl, cl$tree)
  th <- niche_start_ranges(d, "noncentered_bounded")$center
  expect_true(is.finite(niche_logpost(th, d, "noncentered_bounded")))
  # drive all log_sigma_raw strongly negative -> sigma below floor
  idx <- grep("log_sigma_raw", niche_theta_layout(4, 2, "noncentered_bounded"))
  th2 <- th; th2[idx] <- -30
  expect_equal(niche_logpost(th2, d, "noncentered_bounded",
                             sigma_floor = 0.02), -Inf)
})

test_that("niche_multistart recovers interior optima deterministically", {
  skip_if_not_installed("ucminfcpp")
  cl <- make_tiny_clade()
  d <- prepare_phylo_niche_data(cl$occ, cl$bgl, cl$tree)
  f1 <- niche_multistart(d, "noncentered_bounded", n_starts = 8,
                         seed = 11, verbose = FALSE)
  f2 <- niche_multistart(d, "noncentered_bounded", n_starts = 8,
                         seed = 11, verbose = FALSE)
  expect_s3_class(f1, "niche_multistart")
  expect_identical(f1$starts, f2$starts)
  expect_equal(f1$map_lp, f2$map_lp, tolerance = 1e-6)
  expect_true(all(f1$boundary$flag == "ok"))
  # recovered tolerances near the simulated sigma = 0.3
  expect_true(all(f1$boundary$sigma1 > 0.15 & f1$boundary$sigma1 < 0.6))
})

test_that("masked species are predicted but contribute no likelihood", {
  skip_if_not_installed("ucminfcpp")
  cl <- make_tiny_clade()
  dm <- prepare_phylo_niche_data(cl$occ, cl$bgl, cl$tree,
                                 mask_species = "sp2")
  fm <- niche_multistart(dm, "noncentered_bounded", n_starts = 8,
                         seed = 11, verbose = FALSE)
  expect_equal(fm$boundary$flag[fm$boundary$species == "sp2"], "masked")
  # the masked species still has finite predicted niche parameters
  row <- fm$boundary[fm$boundary$species == "sp2", ]
  expect_true(is.finite(row$sigma1) && row$sigma1 > 0)
})
