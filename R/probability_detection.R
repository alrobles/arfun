#' Probability of detection under the ellipsoidal niche model
#'
#' Per-point probability of detecting an occurrence, ported from the
#' \pkg{xsdm} detection kernel (\code{xsdm::log_prob_detect}) to the
#' symmetric ellipsoidal niche used by \pkg{arfun}/\pkg{nicherbayes}:
#' the suitability surface is the Mahalanobis kernel of the niche ellipse,
#' passed through a softplus detection link.
#'
#' @param env Numeric matrix or data.frame of environmental values,
#'   \eqn{n \times p} (rows = points, columns = variables).
#' @param mu Numeric vector of length \eqn{p}: niche optimum (ellipse centre).
#' @param sigma Either a \eqn{p \times p} covariance matrix \code{Sigma}, or a
#'   numeric vector of length \eqn{p} of per-axis standard deviations
#'   (then combined with \code{rho}).
#' @param rho Scalar correlation in \eqn{(-1, 1)}; used only when \code{sigma}
#'   is a vector. Default 0.
#' @param ctil Scalar detection-link offset: larger values push the
#'   detection curve toward zero everywhere (harder to detect).
#'   Default 0.
#' @param pd Maximum detection probability in \eqn{(0, 1]} at the niche
#'   optimum. Default 1.
#'
#' @return Numeric vector of length \code{nrow(env)}: detection
#'   probabilities.
#'
#' @details The detection probability is
#'   \deqn{p(x) = \mathrm{pd} \, / \, (1 + \exp(\mathrm{ctil} + h(x)))}
#'   with \eqn{h(x) = \tfrac12 (x-\mu)^\top \Sigma^{-1} (x-\mu)}, the
#'   half-squared Mahalanobis distance. This matches the xsdm kernel
#'   \code{log p = log pd - log1pexp(ctil + h)} with symmetric widths.
#'   With \code{pd = 1} and \code{ctil -> -Inf} the probability tends to
#'   the pure J&S suitability \eqn{\exp(-0.5\,q^2)} up to normalization.
#'
#' @seealso [virtual_species_occurrences()] to sample occurrences from
#'   this surface.
#' @export
#' @examples
#' env <- expand.grid(x = seq(-3, 3, 0.5), y = seq(-3, 3, 0.5))
#' probability_detection(as.matrix(env), mu = c(0, 0),
#'                       sigma = c(1, 1.5), rho = 0.5, ctil = 2, pd = 0.9)
probability_detection <- function(env, mu, sigma, rho = 0,
                                  ctil = 0, pd = 1) {
  env <- as.matrix(env)
  p <- ncol(env)
  stopifnot(length(mu) == p, pd > 0, pd <= 1, is.finite(ctil))
  Sigma <- if (is.matrix(sigma)) {
    sigma
  } else {
    stopifnot(length(sigma) == p, abs(rho) < 1)
    D <- diag(sigma)
    R <- diag(p); if (p >= 2) { R[1, 2] <- R[2, 1] <- rho }
    D %*% R %*% D
  }
  L <- t(chol(Sigma))
  q <- forwardsolve(L, t(env) - mu)        # whitened coords
  h <- 0.5 * colSums(q^2)                  # half Mahalanobis^2
  pd / (1 + exp(ctil + h))                 # softplus link (xsdm kernel)
}

#' Sample occurrence points for a virtual species
#'
#' Draws occurrence coordinates from the niche intensity over an
#' environmental domain — the generative process assumed by the
#' Jetz-and-Styles-style likelihood used in \pkg{nicherbayes}: intensity
#' \eqn{\propto \exp(-0.5\,q^2(x))} truncated to the domain. An optional
#' detection layer thins the draws through
#' [probability_detection()].
#'
#' @param env Numeric matrix \eqn{n \times p}: the environmental domain
#'   (e.g., a grid over the accessible area M).
#' @param mu Niche centre, length \eqn{p}.
#' @param sigma Covariance matrix \eqn{p \times p} or per-axis sds
#'   (with \code{rho}).
#' @param rho Correlation, used when \code{sigma} is a vector.
#' @param n_occ Number of occurrences to draw.
#' @param ctil,pd Optional detection-layer parameters. When \code{ctil} is
#'   \code{NULL} (default) occurrences are drawn with probability
#'   proportional to the niche intensity only (pure J&S process).
#'
#' @return A matrix with \code{n_occ} rows and \eqn{p} columns of sampled
#'   environmental coordinates.
#'
#' @details
#' Without the detection layer this is exactly the process simulated by
#' the recovery harness in \code{arfun-paper-scripts}
#' (\code{11_sim_recovery.R}): occurrence intensity proportional to the
#' niche kernel over the shared domain. With \code{ctil}/\code{pd} set,
#' each candidate point is kept with probability
#' \code{probability_detection(x)} — a detection-thinned variant useful
#' for realistic virtual-species benchmarks.
#'
#' @export
#' @examples
#' env <- as.matrix(expand.grid(x = seq(-5, 5, 0.2), y = seq(-5, 5, 0.2)))
#' occ <- virtual_species_occurrences(env, mu = c(0.5, -0.3),
#'                                  sigma = c(1, 1.4), rho = -0.4,
#'                                  n_occ = 200)
#' plot(env, pch = ".", col = "grey"); points(occ, col = "red")
virtual_species_occurrences <- function(env, mu, sigma, rho = 0,
                                        n_occ, ctil = NULL, pd = 1) {
  env <- as.matrix(env); p <- ncol(env)
  Sigma <- if (is.matrix(sigma)) sigma else {
    D <- diag(sigma); R <- diag(p); if (p >= 2) { R[1, 2] <- R[2, 1] <- rho }
    D %*% R %*% D
  }
  L <- t(chol(Sigma))
  h <- 0.5 * colSums(forwardsolve(L, t(env) - mu)^2)
  w <- exp(-h); w <- w / sum(w)            # J&S intensity over domain
  out <- matrix(NA_real_, 0, p)
  while (nrow(out) < n_occ) {
    cand <- env[sample.int(nrow(env), n_occ * 2, replace = TRUE,
                           prob = w), , drop = FALSE]
    if (!is.null(ctil)) {
      keep <- stats::runif(nrow(cand)) <
        probability_detection(cand, mu, sigma, rho, ctil = ctil, pd = pd)
      cand <- cand[keep, , drop = FALSE]
    }
    out <- rbind(out, cand)
  }
  colnames(out) <- colnames(env)
  out[seq_len(n_occ), , drop = FALSE]
}
