#' Extract niche parameters from a fitted presence-only model
#'
#' Extracts posterior summaries of the niche centroid (\code{mu}), standard
#' deviations (\code{sigma}), covariance matrix (\code{Sigma}) and correlation
#' matrix (\code{R_corr}) from a fitted Bayesian presence-only niche model.
#'
#' @param fit A \code{CmdStanMCMC} object returned by \code{fit_niche}.
#' @param prob Numeric, credible interval probability (default 0.95).
#'
#' @return A list with components:
#'   \describe{
#'     \item{mu}{Data frame with posterior summary of the centroid
#'       (mean, sd, and credible interval for each dimension).}
#'     \item{sigma}{Data frame with posterior summary of the standard
#'       deviations (mean, sd, and credible interval for each dimension).}
#'     \item{Sigma}{Posterior mean of the covariance matrix (P x P matrix).}
#'     \item{R_corr}{Posterior mean of the correlation matrix (P x P matrix).}
#'   }
#'
#' @examples
#' \dontrun{
#' fit <- fit_niche(occ, chains = 2)
#' pars <- extract_pars(fit)
#' pars$mu       # centroid posterior summary
#' pars$sigma    # tolerance posterior summary
#' pars$Sigma    # posterior mean covariance matrix
#' pars$R_corr   # posterior mean correlation matrix
#'
#' # Geographic interpretation: project the 95% Mahalanobis ellipsoid
#' # onto the environmental raster to get spatial habitat suitability.
#' }
#'
#' @export
extract_pars <- function(fit, prob = 0.95) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop(
      "Package 'posterior' is required for parameter extraction.\n",
      "Install it with: install.packages('posterior')",
      call. = FALSE
    )
  }

  draws <- fit$draws(format = "df")
  alpha <- (1 - prob) / 2

  # Extract mu parameters
  mu_cols <- grep("^mu\\[", names(draws), value = TRUE)
  mu_draws <- as.matrix(draws[, mu_cols])
  mu_summary <- data.frame(
    dimension = seq_len(ncol(mu_draws)),
    mean = colMeans(mu_draws),
    sd = apply(mu_draws, 2, stats::sd),
    lower = apply(mu_draws, 2, stats::quantile, probs = alpha),
    upper = apply(mu_draws, 2, stats::quantile, probs = 1 - alpha),
    row.names = NULL
  )
  names(mu_summary)[4:5] <- c(
    paste0("q", alpha * 100),
    paste0("q", (1 - alpha) * 100)
  )

  # Extract sigma parameters
  sigma_cols <- grep("^sigma\\[", names(draws), value = TRUE)
  sigma_draws <- as.matrix(draws[, sigma_cols])
  sigma_summary <- data.frame(
    dimension = seq_len(ncol(sigma_draws)),
    mean = colMeans(sigma_draws),
    sd = apply(sigma_draws, 2, stats::sd),
    lower = apply(sigma_draws, 2, stats::quantile, probs = alpha),
    upper = apply(sigma_draws, 2, stats::quantile, probs = 1 - alpha),
    row.names = NULL
  )
  names(sigma_summary)[4:5] <- c(
    paste0("q", alpha * 100),
    paste0("q", (1 - alpha) * 100)
  )

  # Extract covariance matrix (posterior mean)
  sigma_mat_cols <- grep("^Sigma\\[", names(draws), value = TRUE)
  P <- length(mu_cols)
  Sigma_mean <- matrix(0, P, P)
  for (col_name in sigma_mat_cols) {
    idx <- as.integer(
      regmatches(col_name, gregexpr("[0-9]+", col_name))[[1]]
    )
    Sigma_mean[idx[1], idx[2]] <- mean(draws[[col_name]])
  }

  # Extract correlation matrix (posterior mean)
  corr_cols <- grep("^R_corr\\[", names(draws), value = TRUE)
  R_corr_mean <- matrix(0, P, P)
  for (col_name in corr_cols) {
    idx <- as.integer(
      regmatches(col_name, gregexpr("[0-9]+", col_name))[[1]]
    )
    R_corr_mean[idx[1], idx[2]] <- mean(draws[[col_name]])
  }

  list(
    mu = mu_summary,
    sigma = sigma_summary,
    Sigma = Sigma_mean,
    R_corr = R_corr_mean
  )
}
