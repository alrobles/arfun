# Regression tests for the warm-start defects documented in the
# paper-astra audit (arfun-paper/paper-astra-audit.R):
#
# (iii) atanh(rho / rho_cap) left its domain when the per-species MLE
#       correlation exceeded rho_cap (|rho| could reach 0.999 while
#       rho_cap = 0.98), injecting NaN into the start vector.
# (iv)  raw residuals of masked/unfitted species were zeroed AFTER the
#       dense solve through L_C, which shifted the implied traits of
#       every *fitted* species below them in the Cholesky ordering.

skip_if_not_installed("ape")
skip_if_not_installed("xnicher")

# Small clade with shared paths. The masked species must sit at position
# >= 2 in the C ordering: masking the first species is a no-op under the
# defect (its solved coordinate is already 0), while masking a middle
# species leaves a nonzero z that used to be zeroed post hoc, shifting
# every species below it.
make_audit_clade <- function(n_occ = 40, seed = 7) {
  set.seed(seed)
  tr <- ape::read.tree(text = "(((sp1:1,sp2:1):1,sp3:2):1,sp4:3);")
  occ <- lapply(seq_len(4), function(s)
    matrix(rnorm(n_occ * 2, rep(c(0.3 * s, -0.2 * s), each = n_occ), 0.3),
           n_occ, 2))
  names(occ) <- tr$tip.label
  bg <- matrix(runif(2000 * 2, -2, 2), 2000, 2)
  bgl <- lapply(seq_len(4), function(s) bg)
  names(bgl) <- tr$tip.label
  list(occ = occ, bgl = bgl, tree = tr)
}

# Reconstruct the implied trait vectors from a theta in
# niche_theta_layout("bnd") order: trait = anc + sqrt(rate) * (L_C %*% raw).
implied_traits <- function(theta, data) {
  P <- data$P
  grab <- function(pat) theta[grepl(pat, names(theta))]
  anc  <- theta[paste0("mu_anc[", seq_len(P), "]")]
  rate <- exp(theta[paste0("lrate_mu[", seq_len(P), "]")])
  raw  <- matrix(grab("^mu_raw"), data$S, P)
  L_C  <- t(chol(data$C))
  mu_impl <- sweep(sweep(L_C %*% raw, 2, sqrt(rate), `*`), 2, anc, `+`)

  ls_anc  <- theta[paste0("log_sigma_anc[", seq_len(P), "]")]
  ls_rate <- exp(theta[paste0("lrate_ls[", seq_len(P), "]")])
  ls_raw  <- matrix(grab("^log_sigma_raw"), data$S, P)
  ls_impl <- sweep(sweep(L_C %*% ls_raw, 2, sqrt(ls_rate), `*`),
                   2, ls_anc, `+`)

  zr_anc  <- theta["z_rho_anc"]
  zr_rate <- exp(theta["lrate_rho"])
  zr_raw  <- grab("^z_rho_raw")
  zr_impl <- zr_anc + sqrt(zr_rate) * drop(L_C %*% zr_raw)

  list(mu = mu_impl, ls = ls_impl, zr = zr_impl)
}

test_that("warm start stays finite when the MLE correlation exceeds rho_cap", {
  cl <- make_audit_clade()
  d <- prepare_phylo_niche_data(cl$occ, cl$bgl, cl$tree)

  # deterministic per-species fit whose implied correlation is ~0.9993
  # (v = 8 -> rho = tanh(v / 2) > rho_cap = 0.98)
  local_mocked_bindings(
    optimize_niche = function(...) {
      list(best = list(theta = c(0.1, -0.2, log(0.4), log(0.5), 8),
                       convergence = 1L))
    },
    .package = "xnicher"
  )

  ws <- niche_warmstart(d, num_starts = 2)
  expect_true(all(is.finite(ws$theta)))
  expect_true(all(is.finite(ws$per_species$rho)))
  # implied correlation respects the cap
  expect_true(all(abs(ws$per_species$rho) <= 0.98 + 1e-9))
})

test_that("masking a species does not shift implied fitted traits", {
  cl <- make_audit_clade()
  d <- prepare_phylo_niche_data(cl$occ, cl$bgl, cl$tree,
                                mask_species = "sp2")

  i <- 0L
  local_mocked_bindings(
    optimize_niche = function(...) {
      i <<- i + 1L
      list(best = list(theta = c(0.4 * i, -0.3 * i,
                                 log(0.4 + 0.1 * i), log(0.5 + 0.1 * i),
                                 0.5),
                       convergence = 1L))
    },
    .package = "xnicher"
  )

  ws <- niche_warmstart(d, num_starts = 2)

  # The implied traits must equal the per-species targets EXACTLY for
  # every species -- including the masked one. Under the defect, zeroing
  # the solved coordinates post hoc displaced the fitted species.
  impl <- implied_traits(ws$theta, d)
  ps   <- ws$per_species
  expect_equal(unname(impl$mu[, 1]), ps$mu1, tolerance = 1e-10)
  expect_equal(unname(impl$mu[, 2]), ps$mu2, tolerance = 1e-10)
  expect_equal(unname(exp(impl$ls[, 1])), ps$sigma1, tolerance = 1e-10)
  expect_equal(unname(exp(impl$ls[, 2])), ps$sigma2, tolerance = 1e-10)
  expect_equal(unname(0.98 * tanh(impl$zr)), ps$rho, tolerance = 1e-10)
})

test_that("masked species start at the phylogenetic conditional mean", {
  cl <- make_audit_clade()
  d <- prepare_phylo_niche_data(cl$occ, cl$bgl, cl$tree,
                                mask_species = "sp2")

  i <- 0L
  local_mocked_bindings(
    optimize_niche = function(...) {
      i <<- i + 1L
      list(best = list(theta = c(0.4 * i, -0.3 * i,
                                 log(0.5), log(0.6), 0.5),
                       convergence = 1L))
    },
    .package = "xnicher"
  )

  ws <- niche_warmstart(d, num_starts = 2)
  impl <- implied_traits(ws$theta, d)
  ps   <- ws$per_species
  uf <- which(!ps$fitted); f <- which(ps$fitted)

  # E[y_uf | y_f] under the warm-start covariance (rate cancels):
  # anc + C_uf C_ff^{-1} (y_f - anc)
  mu_anc <- ws$theta[paste0("mu_anc[", 1:2, "]")]
  W <- d$C[uf, f, drop = FALSE] %*% solve(d$C[f, f, drop = FALSE])
  expected <- outer(rep(1, length(uf)), mu_anc) +
    W %*% sweep(impl$mu[f, , drop = FALSE], 2, mu_anc)
  expect_equal(unname(impl$mu[uf, , drop = FALSE]), unname(expected),
               tolerance = 1e-10)
})
