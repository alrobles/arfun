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
  qnames <- c(paste0("q", alpha * 100), paste0("q", (1 - alpha) * 100))

  # Posterior summary per dimension for a vector parameter
  vec_summary <- function(prefix) {
    cols <- grep(paste0("^", prefix, "\\["), names(draws), value = TRUE)
    d <- as.matrix(draws[, cols])
    out <- data.frame(
      dimension = seq_len(ncol(d)),
      mean = colMeans(d),
      sd = apply(d, 2, stats::sd),
      lower = apply(d, 2, stats::quantile, probs = alpha),
      upper = apply(d, 2, stats::quantile, probs = 1 - alpha),
      row.names = NULL
    )
    names(out)[4:5] <- qnames
    out
  }

  # Posterior-mean matrix for a matrix parameter (Stan "[i,j]" columns)
  mat_mean <- function(prefix, P) {
    cols <- grep(paste0("^", prefix, "\\["), names(draws), value = TRUE)
    m <- matrix(0, P, P)
    for (nm in cols) {
      idx <- as.integer(regmatches(nm, gregexpr("[0-9]+", nm))[[1]])
      m[idx[1], idx[2]] <- mean(draws[[nm]])
    }
    m
  }

  mu_summary <- vec_summary("mu")
  sigma_summary <- vec_summary("sigma")
  P <- nrow(mu_summary)
  Sigma_mean <- mat_mean("Sigma", P)
  R_corr_mean <- mat_mean("R_corr", P)

  list(
    mu = mu_summary,
    sigma = sigma_summary,
    Sigma = Sigma_mean,
    R_corr = R_corr_mean
  )
}
