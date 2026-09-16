#' Prepare niche data for Stan
#'
#' Converts occurrence environmental data into the list format required by
#' the Stan presence-only model. The function validates inputs and builds
#' a data-adaptive prior on the niche centroid.
#'
#' @param env_occ Data frame or matrix of environmental values at presence
#'   points (rows = points, columns = environmental variables).
#' @param eta LKJ shape parameter for correlation matrix prior (default 1,
#'   uniform over correlation matrices).
#' @param sigma_rate Rate parameter for exponential prior on standard
#'   deviations (default 0.1).
#'
#' @return A named list suitable for passing to \code{fit_niche} or directly
#'   to CmdStan. Contains:
#'   \describe{
#'     \item{N_occ}{Number of occurrence points}
#'     \item{P}{Number of environmental dimensions}
#'     \item{occ}{Occurrence data matrix}
#'     \item{mu_prior}{Prior mean for centroid (data-adaptive, equal to
#'       column means of env_occ)}
#'     \item{sigma_prior}{Prior SD for centroid (data-adaptive, equal to
#'       10 times column SDs of env_occ)}
#'     \item{sigma_rate}{Rate for exponential prior on sigma (default 0.1)}
#'     \item{eta}{LKJ shape parameter (default 1, uniform over
#'       correlation matrices)}
#'   }
#'
#' @examples
#' \dontrun{
#' # Simulate some 2D environmental data
#' set.seed(42)
#' occ <- data.frame(
#'   bio01 = rnorm(50, 20, 3),   # annual mean temperature
#'   bio12 = rnorm(50, 800, 100) # annual precipitation
#' )
#'
#' # Prepare for Stan
#' stan_data <- prepare_niche_data(occ)
#' str(stan_data)
#' }
#'
#' @export
prepare_niche_data <- function(env_occ,
                               eta = 1,
                               sigma_rate = 0.1) {
  if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta) ||
      eta <= 0) {
    stop("eta must be a positive finite numeric scalar")
  }
  if (!is.numeric(sigma_rate) || length(sigma_rate) != 1L ||
      !is.finite(sigma_rate) || sigma_rate <= 0) {
    stop("sigma_rate must be a positive finite numeric scalar")
  }

  env_occ <- as.matrix(env_occ)

  if (!is.numeric(env_occ)) {
    stop("env_occ must contain numeric environmental values")
  }
  if (anyNA(env_occ) || any(!is.finite(env_occ))) {
    stop("env_occ must contain only finite, non-missing values")
  }

  N_occ <- nrow(env_occ)
  P <- ncol(env_occ)

  if (N_occ < 2) {
    stop("Need at least 2 occurrence points")
  }
  if (P < 1) {
    stop("Need at least 1 environmental dimension")
  }

  # Data-adaptive priors: center on data mean, scale by data range
  mu_prior <- colMeans(env_occ)
  sigma_prior <- apply(env_occ, 2, stats::sd) * 10
  sigma_prior <- pmax(sigma_prior, 1e-6)
  # Preserve dimension names for interpretability
  names(mu_prior) <- colnames(env_occ)
  names(sigma_prior) <- colnames(env_occ)

  # Build Stan data list
  list(
    N_occ = N_occ,
    P = P,
    occ = env_occ,
    mu_prior = mu_prior,
    sigma_prior = sigma_prior,
    sigma_rate = sigma_rate,
    eta = eta
  )
}
