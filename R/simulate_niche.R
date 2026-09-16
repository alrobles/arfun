#' Simulate niche evolution along a phylogeny
#'
#' Generates virtual-clade niche parameters (\eqn{\mu}, \eqn{\log\sigma},
#' \eqn{z_\rho = \mathrm{atanh}(\rho)}) under Brownian motion,
#' Ornstein--Uhlenbeck, or an independent-species null process. The output
#' matches the parameterization estimated by the joint models and consumed by
#' \code{\link{virtual_species_occurrences}()}, enabling reverse-engineering
#' benchmarks: simulate a niche-evolution process inside a real accessible
#' area, then check that \code{\link{niche_multistart}()} /
#' \code{\link{fit_evolution}()} recover it.
#'
#' @param tree A \code{phylo} object.
#' @param process One of \code{"bm"} (Brownian motion), \code{"ou"}
#'   (Ornstein--Uhlenbeck toward the ancestral optimum), or \code{"null"}
#'   (independent species: every tip draws its niche independently of the
#'   phylogeny).
#' @param mu_anc Length-2 ancestral centroid.
#' @param log_sigma_anc Length-2 ancestral log tolerances.
#' @param z_rho_anc Ancestral Fisher-z correlation.
#' @param rate_mu,rate_log_sigma,rate_rho Evolutionary rates (BM variance per
#'   unit root-to-tip time; under \code{"ou"} the stationary variance is
#'   \code{rate/(2*alpha)}).
#' @param alpha OU constraint strength per trait block (recycled if length 1):
#'   \code{c(mu, log_sigma, z_rho)}. Ignored under \code{"bm"}/\code{"null"}.
#' @param seed Optional RNG seed.
#'
#' @return A list with \code{mu} (S x P), \code{sigma} (S x P, on the natural
#'   scale), \code{rho} (S), \code{tree}, and \code{process}. Row order
#'   matches \code{tree$tip.label}.
#'
#' @details Tips evolve by recursive branch traversal from the root: under BM
#'   each branch adds \code{N(0, rate * branch_length)}; under OU the child
#'   value is \code{theta + (parent - theta) * exp(-alpha * t) + N(0,
#'   rate/(2 alpha) * (1 - exp(-2 alpha t)))}; under \code{"null"} each tip
#'   draws independently of the root and of each other, \code{N(theta,
#'   rate * t_tip)} where \code{t_tip} is its terminal branch length.
#'
#' @seealso \code{\link{virtual_species_occurrences}},
#'   \code{\link{niche_multistart}}
#' @export
simulate_niche_evolution <- function(tree,
                                     process = c("bm", "ou", "null"),
                                     mu_anc = c(0, 0),
                                     log_sigma_anc = log(c(0.4, 0.4)),
                                     z_rho_anc = 0,
                                     rate_mu = c(0.3, 0.3),
                                     rate_log_sigma = c(0.3, 0.3),
                                     rate_rho = 0.3,
                                     alpha = c(1, 1, 1),
                                     seed = NULL) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required", call. = FALSE)
  }
  process <- match.arg(process)
  if (!is.null(seed)) set.seed(seed)
  P <- length(mu_anc)
  S <- length(tree$tip.label)
  alpha <- rep_len(alpha, 3)

  edge <- tree$edge
  el <- tree$edge.length
  n_node <- tree$Nnode
  tip_id <- seq_len(S)

  # trait values at every node (rows), per block
  anc <- list(mu = mu_anc, ls = log_sigma_anc, zr = z_rho_anc)
  rates <- list(mu = rate_mu, ls = rate_log_sigma, zr = rate_rho)
  dims <- c(mu = P, ls = P, zr = 1L)
  node_val <- lapply(dims, function(d)
    matrix(NA_real_, S + n_node, d))
  for (b in names(dims)) node_val[[b]][S + 1, ] <- anc[[b]]

  evolve <- function(parent, t, theta, rate, a) {
    if (process == "ou") {
      m <- theta + (parent - theta) * exp(-a * t)
      v <- rate / (2 * a) * (1 - exp(-2 * a * t))
      m + stats::rnorm(length(parent), 0, sqrt(pmax(v, 0)))
    } else if (process == "bm") {
      parent + stats::rnorm(length(parent), 0, sqrt(rate * t))
    } else {
      parent + stats::rnorm(length(parent), 0, sqrt(rate * t))
    }
  }

  # traverse edges in an order where parents precede children
  ord <- order(edge[, 1])
  for (e in ord) {
    parent <- edge[e, 1]; child <- edge[e, 2]; t <- el[e]
    for (b in names(dims)) {
      theta_b <- anc[[b]]
      rate_b <- rates[[b]]
      a_b <- switch(b, mu = alpha[1], ls = alpha[2], zr = alpha[3])
      if (process == "null") {
        # each tip independent of the root state
        if (child <= S) {
          node_val[[b]][child, ] <- theta_b +
            stats::rnorm(dims[b], 0, sqrt(rate_b * max(t, 1e-8)))
        } else {
          node_val[[b]][child, ] <- theta_b
        }
      } else {
        node_val[[b]][child, ] <-
          evolve(node_val[[b]][parent, ], t, theta_b, rate_b, a_b)
      }
    }
  }

  mu <- node_val$mu[tip_id, , drop = FALSE]
  ls <- node_val$ls[tip_id, , drop = FALSE]
  zr <- node_val$zr[tip_id, 1]
  rownames(mu) <- rownames(ls) <- tree$tip.label
  names(zr) <- tree$tip.label

  list(mu = mu, sigma = exp(ls), rho = tanh(zr), z_rho = zr,
       tree = tree, process = process)
}
