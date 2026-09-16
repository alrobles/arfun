#' Joint phylo-niche log-posterior on the unconstrained scale
#'
#' Exact R implementation of the joint niche-evolution log-posterior used by
#' the \code{bm_evolving} / \code{ou_evolving} Stan models, evaluated on the
#' unconstrained ("math") parameter scale. This is the objective optimized by
#' \code{\link{niche_multistart}()} and the reference for the compiled C++
#' backend (\code{\link{niche_logpost_xptr}()}); both agree to ~1e-12.
#'
#' @param theta Numeric vector in Stan declaration order (see
#'   \code{\link{niche_theta_layout}()} for names and length).
#' @param data A prepared list from \code{\link{prepare_phylo_niche_data}()}.
#' @param parameterization One of \code{"noncentered"} (raw deviations, no
#'   bound on rho), \code{"noncentered_bounded"} (raw deviations,
#'   \code{rho = rho_cap * tanh(z_rho)}; matches the Stan models),
#'   \code{"centered"}, or \code{"centered_lograte"} (centered traits with
#'   declared log-rate parameters and lognormal rate priors).
#' @param sigma_floor Hard lower bound on niche tolerances. The J&S
#'   presence-only likelihood has an unbounded spike as
#'   \eqn{\sigma \to 0} when the background sample does not cover the
#'   occurrence cluster (the "bg-hole" mode); the floor keeps the optimizer
#'   out of that degenerate region. Set to \code{0} to disable (not
#'   recommended on real data).
#' @param rho_cap Bound on \eqn{|\rho|} for the bounded parameterizations.
#' @param debug Print the prior/Jacobian/likelihood decomposition.
#'
#' @return Scalar log-posterior density (up to a constant), or \code{-Inf}
#'   outside the support.
#'
#' @details Parameter layout (Stan declaration order; rate slots hold
#'   \code{log(rate)}):
#' \preformatted{
#'   mu_anc[P] | log_sigma_anc[P] | lrate_mu[P] | lrate_ls[P] |
#'   trait_matrix[S*P col-major] | shape_matrix[S*P] |
#'   z_rho_anc | lrate_rho | z_rho_vec[S]
#' }
#'   For non-centered parameterizations the matrices are \code{mu_raw} and
#'   \code{log_sigma_raw} (standard-normal deviations mapped through
#'   \code{chol(C)}); for centered ones they are \code{mu} and
#'   \code{log_sigma} directly, with \code{MVN(anc, rate*C + jitter*I)}
#'   priors.
#'
#'   Species listed in \code{data$masked} contribute exactly zero to the
#'   likelihood (they carry no occurrence/background rows) but remain in the
#'   phylogenetic prior, so their niche is predicted from their relatives.
#'
#' @seealso \code{\link{niche_multistart}}, \code{\link{niche_logpost_xptr}},
#'   \code{\link{prepare_phylo_niche_data}}
#' @export
niche_logpost <- function(theta, data,
                          parameterization = c("noncentered_bounded",
                                               "noncentered", "centered",
                                               "centered_lograte",
                                               "unb", "bnd", "ctr", "ctrlr"),
                          sigma_floor = 0.02, rho_cap = 0.98,
                          debug = FALSE) {
  model <- match.arg(parameterization)
  model <- switch(model,
                  noncentered = "unb", noncentered_bounded = "bnd",
                  centered = "ctr", centered_lograte = "ctrlr",
                  model)
  pr <- .niche_unpack(theta, data$S, data$P, model)

  lp <- sum(stats::dnorm(pr$mu_anc, data$mu_anc_prior,
                         data$mu_anc_sigma, log = TRUE)) +
        sum(stats::dnorm(pr$log_sigma_anc, 0, 3, log = TRUE)) +
        stats::dnorm(pr$z_rho_anc, data$z_rho_anc_prior,
                     data$z_rho_anc_sigma, log = TRUE)
  if (debug) cat("  lp anc-priors:", lp, "\n")

  if (model == "ctrlr") {
    lp <- lp +
      sum(stats::dnorm(pr$lrate_mu, log(data$rate_mu_scale), 1, log = TRUE)) +
      sum(stats::dnorm(pr$lrate_ls, log(data$rate_log_sigma_scale), 1,
                       log = TRUE)) +
      stats::dnorm(pr$lrate_rho, log(data$rate_rho_scale), 1, log = TRUE)
  } else {
    lp <- lp +
      sum(stats::dnorm(pr$rate_mu, 0, data$rate_mu_scale, log = TRUE)) +
      sum(stats::dnorm(pr$rate_ls, 0, data$rate_log_sigma_scale,
                       log = TRUE)) +
      stats::dnorm(pr$rate_rho, 0, data$rate_rho_scale, log = TRUE) +
      sum(pr$lrate_mu) + sum(pr$lrate_ls) + pr$lrate_rho
    if (debug) cat("  + rate priors+jac:", lp, "\n")
  }

  noncentered <- model %in% c("unb", "bnd")
  Jmat <- diag(1e-9, data$S)
  if (noncentered) {
    lp <- lp + sum(stats::dnorm(pr$mu_raw, log = TRUE)) +
                sum(stats::dnorm(pr$log_sigma_raw, log = TRUE)) +
                sum(stats::dnorm(pr$z_rho_raw, log = TRUE))
    L_C <- t(chol(data$C))
    S <- data$S; P <- data$P
    mu <- matrix(NA_real_, S, P)
    log_sigma <- matrix(NA_real_, S, P)
    for (k in seq_len(P)) {
      mu[, k] <- pr$mu_anc[k] + sqrt(pr$rate_mu[k]) *
                 (L_C %*% pr$mu_raw[, k])
      log_sigma[, k] <- pr$log_sigma_anc[k] + sqrt(pr$rate_ls[k]) *
                        (L_C %*% pr$log_sigma_raw[, k])
    }
    z_rho <- drop(pr$z_rho_anc + sqrt(pr$rate_rho) *
                  (L_C %*% pr$z_rho_raw))
  } else {
    mu <- pr$mu; log_sigma <- pr$log_sigma; z_rho <- pr$z_rho
    S <- data$S; P <- data$P
    for (k in seq_len(P)) {
      lp <- lp + .dmvn_lp(mu[, k], rep(pr$mu_anc[k], S),
                          pr$rate_mu[k] * data$C + Jmat) +
                 .dmvn_lp(log_sigma[, k], rep(pr$log_sigma_anc[k], S),
                          pr$rate_ls[k] * data$C + Jmat)
    }
    lp <- lp + .dmvn_lp(z_rho, rep(pr$z_rho_anc, S),
                        pr$rate_rho * data$C + Jmat)
    if (debug) cat("  + centered MVN priors:", lp, "\n")
  }
  if (!is.finite(lp)) return(-Inf)

  cap <- if (model == "unb") 1 else rho_cap
  rho_v <- cap * tanh(z_rho)
  if (any(!is.finite(rho_v))) return(-Inf)
  if (any(exp(log_sigma) < sigma_floor)) return(-Inf)

  for (s in seq_len(S)) {
    ll <- .species_ll_2d(data, s, mu[s, ], log_sigma[s, ], rho_v[s])
    if (!is.finite(ll)) return(-Inf)
    lp <- lp + ll
  }
  if (debug) cat("  total:", lp, "\n")
  lp
}

#' Compiled C++ objective for the phylo-niche log-posterior
#'
#' Builds an external pointer to the compiled version of
#' \code{\link{niche_logpost}()} (identical math, ~1e-12 agreement) for use
#' with \code{ucminfcpp::ucminf_xptr}. This is the fast path used by
#' \code{\link{niche_multistart}()}; it can also be passed to any optimizer
#' that accepts an \code{ObjFun} external pointer.
#'
#' @param data A list from \code{\link{prepare_phylo_niche_data}()}.
#' @param grad Finite-difference scheme for the optimizer gradient:
#'   \code{"central"} (default, more accurate) or \code{"forward"} (cheaper).
#' @inheritParams niche_logpost
#'
#' @return An external pointer (\code{ucminf::ObjFun}) evaluating
#'   \eqn{-\mathrm{lp}} plus a finite-difference gradient, suitable for
#'   \code{ucminfcpp::ucminf_xptr}. Evaluate the log-posterior at a point
#'   with \code{\link{niche_logpost_eval_xptr}()}.
#' @export
niche_logpost_xptr <- function(data,
                               parameterization = "noncentered_bounded",
                               sigma_floor = 0.02, rho_cap = 0.98,
                               grad = c("central", "forward")) {
  model <- match.arg(parameterization,
                     c("noncentered_bounded", "noncentered", "centered",
                       "centered_lograte", "unb", "bnd", "ctr", "ctrlr"))
  model <- switch(model, noncentered = "unb",
                  noncentered_bounded = "bnd", centered = "ctr",
                  centered_lograte = "ctrlr", model)
  grad <- match.arg(grad)
  if (!requireNamespace("ucminfcpp", quietly = TRUE)) {
    stop("Package 'ucminfcpp' is required for the compiled objective.\n",
         "Install it with: remotes::install_github('alrobles/ucminfcpp')",
         call. = FALSE)
  }
  make_niche_logpost_xptr(
    occ = data$occ, occ_start = data$occ_start, N_occ = data$N_occ,
    env_m = data$env_m, m_start = data$m_start, N_m = data$N_m,
    C = data$C,
    mu_anc_prior = data$mu_anc_prior, mu_anc_sigma = data$mu_anc_sigma,
    z_rho_anc_prior = data$z_rho_anc_prior,
    z_rho_anc_sigma = data$z_rho_anc_sigma,
    rate_mu_scale = data$rate_mu_scale,
    rate_ls_scale = data$rate_log_sigma_scale,
    rate_rho_scale = data$rate_rho_scale,
    model = model, rho_cap = rho_cap, sigma_floor = sigma_floor,
    grad = grad)
}

#' Evaluate the compiled phylo-niche log-posterior
#'
#' @param xptr External pointer from \code{\link{niche_logpost_xptr}()}.
#' @param theta Numeric vector in \code{\link{niche_theta_layout}()} order.
#'
#' @return Scalar log-posterior (\code{-Inf} outside support).
#' @export
niche_logpost_eval_xptr <- function(xptr, theta) {
  xptr_niche_logpost_eval(xptr, as.numeric(theta))
}

# ---- internals -------------------------------------------------------------

.lse <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# J&S species log-likelihood (P = 2). Masked species (N_occ == 0) return 0.
.species_ll_2d <- function(data, s, mu_s, log_sigma_s, rho) {
  no <- data$N_occ[s]
  if (no == 0L) return(0)
  sg <- exp(log_sigma_s)
  L <- matrix(c(sg[1], rho * sg[2], 0,
                sg[2] * sqrt(max(1 - rho^2, 1e-12))), 2, 2)
  os <- data$occ_start[s]; oe <- os + no - 1L
  y <- forwardsolve(L, t(data$occ[os:oe, , drop = FALSE]) - mu_s)
  bs <- data$m_start[s]; be <- bs + data$N_m[s] - 1L
  ym <- forwardsolve(L, t(data$env_m[bs:be, , drop = FALSE]) - mu_s)
  -0.5 * sum(y^2) - no * .lse(-0.5 * colSums(ym^2))
}

.dmvn_lp <- function(x, mean_v, Sigma) {
  L <- tryCatch(t(chol(Sigma)), error = function(e) NULL)
  if (is.null(L)) return(-Inf)
  z <- forwardsolve(L, x - mean_v)
  -0.5 * sum(z^2) - sum(log(diag(L))) - length(x) / 2 * log(2 * pi)
}

.niche_unpack <- function(theta, S, P, model) {
  i <- 1
  take <- function(n) {
    v <- theta[i:(i + n - 1)]
    i <<- i + n
    v
  }
  pr <- list()
  pr$mu_anc <- take(P)
  pr$log_sigma_anc <- take(P)
  pr$lrate_mu <- take(P)
  pr$lrate_ls <- take(P)
  pr$rate_mu <- exp(pr$lrate_mu)
  pr$rate_ls <- exp(pr$lrate_ls)
  if (model %in% c("unb", "bnd")) {
    pr$mu_raw <- matrix(take(S * P), S, P)
    pr$log_sigma_raw <- matrix(take(S * P), S, P)
  } else {
    pr$mu <- matrix(take(S * P), S, P)
    pr$log_sigma <- matrix(take(S * P), S, P)
  }
  pr$z_rho_anc <- take(1)
  pr$lrate_rho <- take(1)
  pr$rate_rho <- exp(pr$lrate_rho)
  if (model %in% c("unb", "bnd")) pr$z_rho_raw <- take(S) else pr$z_rho <- take(S)
  if (i - 1 != length(theta)) {
    stop("`theta` has length ", length(theta), " but the '", model,
         "' layout needs ", i - 1, " (S = ", S, ", P = ", P, ")")
  }
  pr
}

#' Parameter layout of the phylo-niche objective
#'
#' Names and order of the \code{theta} vector consumed by
#' \code{\link{niche_logpost}()} / \code{\link{niche_logpost_xptr}()},
#' matching Stan matrix serialization (column-major).
#'
#' @param S Number of species (tree tips).
#' @param P Number of environmental dimensions (currently 2).
#' @inheritParams niche_logpost
#'
#' @return Character vector of parameter names in \code{theta} order.
#' @export
niche_theta_layout <- function(S, P,
                               parameterization = "noncentered_bounded") {
  model <- match.arg(parameterization,
                     c("noncentered_bounded", "noncentered", "centered",
                       "centered_lograte", "unb", "bnd", "ctr", "ctrlr"))
  model <- switch(model, noncentered = "unb",
                  noncentered_bounded = "bnd", centered = "ctr",
                  centered_lograte = "ctrlr", model)
  mat_names <- function(prefix) {
    as.vector(outer(seq_len(S), seq_len(P),
                    function(s, k) paste0(prefix, "[", s, ",", k, "]")))
  }
  trait <- if (model %in% c("unb", "bnd")) "mu_raw" else "mu"
  shape <- if (model %in% c("unb", "bnd")) "log_sigma_raw" else "log_sigma"
  zvec <- if (model %in% c("unb", "bnd")) "z_rho_raw" else "z_rho"
  c(paste0("mu_anc[", seq_len(P), "]"),
    paste0("log_sigma_anc[", seq_len(P), "]"),
    paste0("lrate_mu[", seq_len(P), "]"),
    paste0("lrate_ls[", seq_len(P), "]"),
    mat_names(trait), mat_names(shape),
    "z_rho_anc", "lrate_rho",
    paste0(zvec, "[", seq_len(S), "]"))
}

#' Empirical per-species niche moments
#'
#' Mean, standard deviation, correlation and Fisher-z of the occurrence
#' cluster of each species in a prepared data list -- used for start designs
#' and boundary diagnostics.
#'
#' @param data A list from \code{\link{prepare_phylo_niche_data}()}.
#'
#' @return A list with matrices \code{mu} (S x P), \code{sd} (S x P),
#'   \code{ls} (log sd), and vectors \code{cor}, \code{zr} (Fisher-z),
#'   \code{n_occ}, \code{n_cells}. Masked species get \code{NA} rows.
#' @export
niche_empirical <- function(data) {
  S <- data$S; P <- data$P
  mu <- matrix(NA_real_, S, P); sdev <- mu
  cor_v <- rep(NA_real_, S)
  n_cells <- integer(S)
  for (s in seq_len(S)) {
    no <- data$N_occ[s]
    if (no == 0L) next
    r <- data$occ[data$occ_start[s]:(data$occ_start[s] + no - 1L), ,
                  drop = FALSE]
    mu[s, ] <- colMeans(r)
    sdev[s, ] <- pmax(apply(r, 2, stats::sd), 0.05)
    cor_v[s] <- if (P == 2) stats::cor(r[, 1], r[, 2]) else NA_real_
    n_cells[s] <- nrow(unique(round(r, 3)))
  }
  cor_v <- pmin(pmax(ifelse(is.finite(cor_v), cor_v, 0), -0.9), 0.9)
  list(mu = mu, sd = sdev, ls = log(sdev), cor = cor_v, zr = atanh(cor_v),
       n_occ = data$N_occ, n_cells = n_cells, species = data$species)
}

#' Data-driven start ranges for the phylo-niche objective
#'
#' Builds lower/upper bounds for each unconstrained coordinate from the
#' empirical occurrence moments -- the design space sampled by
#' \code{\link{niche_multistart}()}.
#'
#' @param data A list from \code{\link{prepare_phylo_niche_data}()}.
#' @inheritParams niche_logpost
#'
#' @return A list with \code{lo}, \code{hi} (named numeric vectors in
#'   \code{\link{niche_theta_layout}()} order) and \code{center} (the
#'   data-driven centre of the ranges).
#' @export
niche_start_ranges <- function(data,
                               parameterization = "noncentered_bounded") {
  model <- match.arg(parameterization,
                     c("noncentered_bounded", "noncentered", "centered",
                       "centered_lograte", "unb", "bnd", "ctr", "ctrlr"))
  model <- switch(model, noncentered = "unb",
                  noncentered_bounded = "bnd", centered = "ctr",
                  centered_lograte = "ctrlr", model)
  S <- data$S; P <- data$P
  emp <- niche_empirical(data)
  occ_in <- do.call(rbind, lapply(seq_len(S), function(s) {
    if (data$N_occ[s] == 0L) return(NULL)
    data$occ[data$occ_start[s]:(data$occ_start[s] + data$N_occ[s] - 1L), ,
             drop = FALSE]
  }))
  pool_sd <- pmax(apply(occ_in, 2, stats::sd), 1e-3)
  lo <- c(); hi <- c(); nm <- c()
  add <- function(l, h, n) {
    lo <<- c(lo, l); hi <<- c(hi, h); nm <<- c(nm, n)
  }
  add(apply(occ_in, 2, min), apply(occ_in, 2, max),
      paste0("mu_anc[", seq_len(P), "]"))
  add(log(pool_sd / 4), log(pool_sd * 2),
      paste0("log_sigma_anc[", seq_len(P), "]"))
  add(rep(log(0.02), P), rep(log(2), P),
      paste0("lrate_mu[", seq_len(P), "]"))
  add(rep(log(0.02), P), rep(log(1.5), P),
      paste0("lrate_ls[", seq_len(P), "]"))
  if (model %in% c("unb", "bnd")) {
    add(rep(-3, S * P), rep(3, S * P),
        as.vector(outer(seq_len(S), seq_len(P),
                        function(s, k) paste0("mu_raw[", s, ",", k, "]"))))
    add(rep(-3, S * P), rep(3, S * P),
        as.vector(outer(seq_len(S), seq_len(P),
                        function(s, k) paste0("log_sigma_raw[", s, ",", k, "]"))))
  } else {
    mu_lo <- emp$mu; mu_hi <- emp$mu
    sd_lo <- emp$ls; sd_hi <- emp$ls
    for (s in seq_len(S)) {
      if (data$N_occ[s] == 0L) {
        mu_lo[s, ] <- apply(occ_in, 2, min)
        mu_hi[s, ] <- apply(occ_in, 2, max)
        sd_lo[s, ] <- log(pool_sd / 4)
        sd_hi[s, ] <- log(pool_sd * 2)
      }
    }
    add(as.vector(mu_lo - 2 * emp$sd),
        as.vector(mu_hi + 2 * emp$sd),
        as.vector(outer(seq_len(S), seq_len(P),
                        function(s, k) paste0("mu[", s, ",", k, "]"))))
    add(as.vector(sd_lo - log(3)), as.vector(sd_hi + log(3)),
        as.vector(outer(seq_len(S), seq_len(P),
                        function(s, k) paste0("log_sigma[", s, ",", k, "]"))))
  }
  pc <- atanh(pmin(pmax(stats::cor(occ_in)[1, 2], -0.95), 0.95))
  add(pc - 1, pc + 1, "z_rho_anc")
  add(log(0.02), log(1.5), "lrate_rho")
  if (model %in% c("unb", "bnd")) {
    add(rep(-3, S), rep(3, S),
        paste0("z_rho_raw[", seq_len(S), "]"))
  } else {
    zr_lo <- emp$zr; zr_hi <- emp$zr
    for (s in seq_len(S)) {
      if (data$N_occ[s] == 0L) {
        zr_lo[s] <- pc - 1; zr_hi[s] <- pc + 1
      }
    }
    add(as.vector(zr_lo - 2), as.vector(zr_hi + 2),
        paste0("z_rho[", seq_len(S), "]"))
  }
  # masked species' empirical moments are NA -- widen them to the pool range
  fix_na <- function(v, l, h) ifelse(is.na(v), l, ifelse(is.na(h), v, v))
  lo[is.na(lo)] <- -3; hi[is.na(hi)] <- 3
  names(lo) <- names(hi) <- nm
  bad <- !is.finite(lo) | !is.finite(hi) | lo >= hi
  lo[bad] <- -3; hi[bad] <- 3
  list(lo = lo, hi = hi, center = (lo + hi) / 2)
}
