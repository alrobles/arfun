#' Per-species warm-start for the joint phylo-niche objective
#'
#' Fits each species independently with
#' \code{\link[xnicher]{optimize_niche}} and assembles a starting vector
#' in \code{\link{niche_theta_layout}()} order for the joint model:
#' ancestral states are set to the cross-species means, rates to the
#' cross-species variances, and (for non-centred parameterizations) the
#' raw residuals are solved through the Cholesky factor of \eqn{C} so the
#' implied species values equal the per-species MLEs.
#'
#' This gives \code{\link{niche_multistart}()} a biologically informed
#' start -- each species begins at its own independent niche optimum --
#' alongside the default Sobol/LHS design. Per-species fits that fail to
#' converge, or that land on the degenerate \eqn{\sigma \to 0} mode
#' (far below the empirical spread), fall back to empirical moments;
#' masked or otherwise unfitted species start at the phylogenetic
#' conditional mean given the fitted species (the ancestor mean when no
#' species fits).
#'
#' @param data A list from \code{\link{prepare_phylo_niche_data}()}.
#' @param likelihood \code{xnicher} likelihood family used for the
#'   per-species fits. \code{"weighted"} (default) is the M-restricted
#'   J&S family; \code{"presence_only"} fits the unnormalized model.
#' @param prior_log_sigma_lambda Penalty strength on
#'   \eqn{\log\sigma} passed to \code{\link[xnicher]{optimize_niche}};
#'   regularises the per-species fits away from the degenerate
#'   \eqn{\sigma \to 0} mode. Increase if warm-start species collapse.
#' @param num_starts Number of multi-start points per species (passed to
#'   \code{\link[xnicher]{optimize_niche}}).
#' @param min_cells Minimum number of distinct environmental cells to
#'   attempt an independent fit; species below the threshold (and masked
#'   species) use empirical moments or the ancestor mean.
#' @param verbose Print per-species fit results.
#' @inheritParams niche_logpost
#'
#' @return A list with \code{theta} (named vector in
#'   \code{\link{niche_theta_layout}()} order, ready for
#'   \code{\link{niche_multistart}(extra_starts = ...)}) and
#'   \code{per_species} (data frame of the independent fits).
#' @export
niche_warmstart <- function(data,
                            parameterization = "noncentered_bounded",
                            rho_cap = 0.98,
                            likelihood = "weighted",
                            prior_log_sigma_lambda = 10,
                            num_starts = 50,
                            min_cells = 5,
                            verbose = FALSE) {
  if (!requireNamespace("xnicher", quietly = TRUE)) {
    stop("Package 'xnicher' is required for warm starts. ",
         "Install it from CRAN: install.packages('xnicher')",
         call. = FALSE)
  }
  model <- .canon_param(parameterization)
  S <- data$S; P <- data$P
  emp <- niche_empirical(data)

  mu_hat <- emp$mu
  ls_hat <- emp$ls
  zr_hat <- emp$zr
  fitted <- rep(FALSE, S)

  for (s in seq_len(S)) {
    if (data$N_occ[s] == 0L || emp$n_cells[s] < min_cells) next
    occ_s <- data$occ[data$occ_start[s]:
                      (data$occ_start[s] + data$N_occ[s] - 1L), ,
                      drop = FALSE]
    bg_s  <- data$env_m[data$m_start[s]:
                        (data$m_start[s] + data$N_m[s] - 1L), ,
                        drop = FALSE]
    # try the M-restricted family first; fall back to the unnormalized
    # presence-only fit if it lands on the degenerate sigma -> 0 mode
    th <- NULL; conv <- NA_integer_
    for (lik in unique(c(likelihood, "presence_only"))) {
      fit <- tryCatch(suppressWarnings(
        xnicher::optimize_niche(env_occ = as.data.frame(occ_s),
                                env_m = as.data.frame(bg_s),
                                likelihood = lik,
                                num_starts = num_starts,
                                prior_log_sigma_lambda =
                                  prior_log_sigma_lambda,
                                verbose = FALSE)),
        error = function(e) NULL)
      cand <- if (!is.null(fit)) fit$best$theta else NULL
      if (is.null(cand) || length(cand) < 2 * P || !all(is.finite(cand)))
        next
      ls_cand <- cand[(P + 1):(2 * P)]
      # degenerate fits: sigma implausibly below the empirical spread
      if (any(exp(ls_cand) < 0.25 * exp(emp$ls[s, ]))) next
      th <- cand
      conv <- fit$best$convergence
      break
    }
    if (is.null(th)) next
    mu_hat[s, ] <- th[seq_len(P)]
    ls_hat[s, ] <- th[(P + 1):(2 * P)]
    if (P == 2L && length(th) > 2 * P) {
      L <- xnicher::cvine_cholesky(th[(2 * P + 1):length(th)], d = P)
      # clamp to the model's admissible correlation range BEFORE dividing
      # by rho_cap: |rho| > rho_cap would leave atanh() outside its
      # domain and inject NaN into the start vector
      rho <- tcrossprod(L)[1, 2]
      rho <- sign(rho) * min(abs(rho), rho_cap * (1 - 1e-4))
      zr_hat[s] <- atanh(rho / rho_cap)
    }
    fitted[s] <- TRUE
    if (verbose)
      message("  ", data$species[s], " (conv ", conv, "): mu=(",
              paste(round(mu_hat[s, ], 3), collapse = ","),
              ") sigma=(", paste(round(exp(ls_hat[s, ]), 3),
                                 collapse = ","), ")")
  }

  use <- which(emp$n_cells >= min_cells & is.finite(mu_hat[, 1]))
  if (!length(use)) use <- seq_len(S)
  mu_anc <- colMeans(mu_hat[use, , drop = FALSE])
  ls_anc <- colMeans(ls_hat[use, , drop = FALSE])
  zr_anc <- mean(zr_hat[use])
  rate_mu <- pmax(apply(mu_hat[use, , drop = FALSE], 2, stats::var), 0.01)
  rate_ls <- pmax(apply(ls_hat[use, , drop = FALSE], 2, stats::var), 0.01)
  rate_rho <- max(stats::var(zr_hat[use]), 0.01)

  # unfitted/masked species start at the phylogenetic conditional mean
  # given the fitted species, E[y_uf | y_f] = anc + C_uf C_ff^{-1}(y_f -
  # anc): the BM/OU rate cancels in the conditioning weights. Falls back
  # to the ancestor mean when no species fitted or C_ff is singular.
  for (s in which(!fitted)) {
    mu_hat[s, ] <- mu_anc
    ls_hat[s, ] <- ls_anc
    zr_hat[s] <- zr_anc
  }
  fit_idx <- which(fitted)
  uf <- which(!fitted)
  if (length(fit_idx) && length(uf)) {
    W <- tryCatch(
      data$C[uf, fit_idx, drop = FALSE] %*%
        solve(data$C[fit_idx, fit_idx, drop = FALSE]),
      error = function(e) NULL)
    if (!is.null(W) && all(is.finite(W))) {
      mu_hat[uf, ] <- outer(rep(1, length(uf)), mu_anc) +
        W %*% sweep(mu_hat[fit_idx, , drop = FALSE], 2, mu_anc)
      ls_hat[uf, ] <- outer(rep(1, length(uf)), ls_anc) +
        W %*% sweep(ls_hat[fit_idx, , drop = FALSE], 2, ls_anc)
      zr_hat[uf] <- zr_anc + drop(W %*% (zr_hat[fit_idx] - zr_anc))
    }
  }

  L_C <- t(chol(data$C))
  if (.is_noncentered(model)) {
    # Solve the raw residuals ONCE against the complete target vector:
    # L_C is dense, so zeroing solved coordinates post hoc would also
    # shift the implied traits of every species below the masked one.
    mu_raw <- matrix(0, S, P); ls_raw <- matrix(0, S, P)
    for (k in seq_len(P)) {
      mu_raw[, k] <- solve(L_C, (mu_hat[, k] - mu_anc[k]) /
                               sqrt(rate_mu[k]))
      ls_raw[, k] <- solve(L_C, (ls_hat[, k] - ls_anc[k]) /
                               sqrt(rate_ls[k]))
    }
    zr_raw <- solve(L_C, (zr_hat - zr_anc) / sqrt(rate_rho))
    theta <- c(mu_anc, ls_anc, log(rate_mu), log(rate_ls),
               as.vector(mu_raw), as.vector(ls_raw),
               zr_anc, log(rate_rho), zr_raw)
  } else {
    theta <- c(mu_anc, ls_anc, log(rate_mu), log(rate_ls),
               as.vector(mu_hat), as.vector(ls_hat),
               zr_anc, log(rate_rho), zr_hat)
  }
  names(theta) <- niche_theta_layout(S, P, model)

  per_species <- data.frame(
    species = data$species, fitted = fitted,
    n_cells = emp$n_cells,
    mu1 = mu_hat[, 1], sigma1 = exp(ls_hat[, 1]),
    sigma2 = if (P > 1) exp(ls_hat[, 2]) else NA_real_,
    rho = if (P == 2) rho_cap * tanh(zr_hat) else NA_real_)
  list(theta = theta, per_species = per_species)
}
