#' Multi-start MAP optimization of the joint phylo-niche posterior
#'
#' Optimizes \code{\link{niche_logpost}()} from a space-filling design over
#' \code{\link{niche_start_ranges}()}. By default it uses the compiled C++
#' objective through \code{ucminfcpp::ucminf_xptr} (finite-difference
#' gradients, no R callbacks); set \code{use_xptr = FALSE} to optimize the R
#' implementation instead. Start points that are not finite under the
#' objective are skipped. The result includes per-species boundary flags so
#' collapsed (\eqn{\sigma \to \sigma_{\mathrm{floor}}}), sliver-like
#' (\eqn{|\rho| \to \rho_{\mathrm{cap}}}) and inflated
#' (\eqn{\sigma \to \infty}) niches are reported explicitly rather than
#' silently contaminating downstream fits.
#'
#' @param data A list from \code{\link{prepare_phylo_niche_data}()}.
#' @param n_starts Integer number of start points (Sobol/LHS design plus one
#'   data-driven centre). Increase for large clades.
#' @param use_xptr Logical; use the compiled C++ objective via
#'   \code{ucminfcpp::ucminf_xptr} (default \code{TRUE}). Falls back to the R
#'   objective and \code{ucminf} if the compiled backend or package is
#'   unavailable.
#' @param seed Integer seed for the start design (LHS fallback only; Sobol is
#'   deterministic).
#' @param control List of control arguments passed to the optimizer
#'   (defaults: \code{grtol = 1e-5, xtol = 1e-10, stepmax = 5,
#'   maxeval = 3000}).
#' @param verbose Print per-start progress.
#' @inheritParams niche_logpost
#'
#' @return A \code{niche_multistart} object: a list with \code{map} (named
#'   best \code{theta}), \code{map_lp}, \code{solutions} (per-start
#'   \code{par}, \code{lp}, \code{convergence}, \code{maxgradient},
#'   \code{start_lp}), \code{starts}, \code{boundary} (per-species boundary
#'   report at the MAP; see \code{\link{niche_boundary_report}()}),
#'   \code{parameterization}, and the \code{sigma_floor}/\code{rho_cap} used.
#' @export
niche_multistart <- function(data,
                             parameterization = c("noncentered_bounded",
                                                  "noncentered", "centered",
                                                  "centered_lograte",
                                                  "unb", "bnd", "ctr",
                                                  "ctrlr"),
                             n_starts = 32,
                             sigma_floor = 0.02, rho_cap = 0.98,
                             use_xptr = TRUE, seed = 1,
                             control = list(), verbose = TRUE) {
  parameterization <- match.arg(parameterization)
  model <- switch(parameterization,
                  noncentered = "unb", noncentered_bounded = "bnd",
                  centered = "ctr", centered_lograte = "ctrlr",
                  parameterization)
  rng <- niche_start_ranges(data, model)
  n_par <- length(rng$lo)
  theta_names <- names(rng$lo)

  starts <- .niche_start_design(rng$lo, rng$hi, n_starts - 1, seed)
  starts <- rbind(starts, rng$center)

  lp_r <- function(th) {
    niche_logpost(th, data, parameterization = model,
                  sigma_floor = sigma_floor, rho_cap = rho_cap)
  }

  xptr <- NULL
  if (use_xptr) {
    xptr <- tryCatch(
      .make_xptr_objective(data, model, sigma_floor, rho_cap),
      error = function(e) {
        warning("C++ objective unavailable (", conditionMessage(e),
                "); falling back to the R objective", call. = FALSE)
        NULL
      })
    if (!is.null(xptr)) {
      probe <- rbind(rng$center,
                     starts[seq_len(min(3, nrow(starts) - 1)), ,
                            drop = FALSE])
      dmax <- 0
      for (ii in seq_len(nrow(probe))) {
        lr <- lp_r(probe[ii, ])
        lc <- xptr_niche_logpost_eval(xptr, probe[ii, ])
        dd <- if (is.finite(lr) && is.finite(lc)) abs(lr - lc)
              else if (identical(lr, lc)) 0 else Inf
        dmax <- max(dmax, dd)
      }
      if (!is.finite(dmax) || dmax > 0.5) {
        warning("C++ objective deviates from R by ", dmax,
                " at probe points; using the R objective", call. = FALSE)
        xptr <- NULL
      }
    }
  }
  use_xptr <- !is.null(xptr)

  start_lp <- apply(starts, 1, function(th)
    tryCatch(lp_r(th), error = function(e) -Inf))
  keep <- is.finite(start_lp)
  if (verbose) {
    message("niche_multistart: ", sum(keep), "/", nrow(starts),
            " finite starts | objective: ",
            if (use_xptr) "C++ (ucminf_xptr)" else "R (ucminf)")
  }
  if (!any(keep)) {
    stop("No finite start points; check data ranges and 'sigma_floor'")
  }

  ctrl <- utils::modifyList(list(grtol = 1e-5, xtol = 1e-10, stepmax = 5,
                          maxeval = 3000), control)
  obj <- function(th) {
    v <- -lp_r(th)
    if (!is.finite(v)) 1e12 else v
  }
  sols <- vector("list", nrow(starts))
  ord <- order(start_lp, decreasing = TRUE)
  for (i in ord) {
    if (!is.finite(start_lp[i])) next
    r <- tryCatch(
      if (use_xptr) {
        ucminfcpp::ucminf_xptr(par = unname(starts[i, ]), xptr = xptr,
                               control = ctrl)
      } else {
        ucminf_compat(par = starts[i, ], fn = obj, control = ctrl)
      },
      error = function(e) NULL)
    if (is.null(r)) next
    sols[[i]] <- list(par = r$par, lp = -r$value, conv = r$convergence,
                      maxgrad = r$info[["maxgradient"]],
                      start_i = i, start_lp = start_lp[i])
    if (verbose) {
      message(sprintf("start %2d -> lp %.2f conv %s maxgrad %.2e",
                      i, -r$value, r$convergence,
                      ifelse(is.null(r$info[["maxgradient"]]), NA_real_,
                             r$info[["maxgradient"]])))
    }
  }
  sols <- sols[!vapply(sols, is.null, logical(1))]
  if (length(sols) == 0) stop("All optimizations failed")
  sols <- sols[order(vapply(sols, `[[`, numeric(1), "lp"),
                     decreasing = TRUE)]
  map <- sols[[1]]$par
  names(map) <- theta_names
  names(start_lp) <- NULL

  bnd <- niche_boundary_report(data, map, model, sigma_floor, rho_cap)
  if (verbose && any(bnd$flag != "ok")) {
    message("boundary flags: ",
            paste(bnd$species[bnd$flag != "ok"],
                  bnd$flag[bnd$flag != "ok"], sep = ":", collapse = ", "))
  }

  out <- list(map = map, map_lp = sols[[1]]$lp, solutions = sols,
              starts = starts, start_lp = start_lp, boundary = bnd,
              parameterization = model, sigma_floor = sigma_floor,
              rho_cap = rho_cap, species = data$species,
              masked = data$masked)
  class(out) <- "niche_multistart"
  out
}

# Sobol design if the 'sobol' package is available, else scrambled LHS.
.niche_start_design <- function(lo, hi, n, seed) {
  n_par <- length(lo)
  if (n <= 0) return(matrix(nrow = 0, ncol = n_par))
  if (requireNamespace("sobol", quietly = TRUE)) {
    return(as.matrix(sobol::sobol_design(lower = lo, upper = hi, nseq = n)))
  }
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit({
      if (is.null(old)) rm(".Random.seed", envir = .GlobalEnv)
      else assign(".Random.seed", old, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  u <- matrix(stats::runif(n * n_par), n, n_par)
  for (j in seq_len(n_par)) u[, j] <- u[order(stats::runif(n)), j]
  sweep(u, 2, hi - lo, `*`) + rep(lo, each = n)
}

# Build the compiled objective via the Rcpp-exported factory.
.make_xptr_objective <- function(data, model, sigma_floor, rho_cap) {
  if (!requireNamespace("ucminfcpp", quietly = TRUE)) {
    stop("package 'ucminfcpp' is not installed")
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
    grad = "central")
}

# Thin wrapper so niche_multistart works with either the Fortran 'ucminf'
# package or a stats::optim fallback when neither is installed.
ucminf_compat <- function(par, fn, control = list()) {
  if (requireNamespace("ucminf", quietly = TRUE)) {
    ctrl <- control
    ctrl$grad <- ctrl$grad %||% "forward"
    return(ucminf::ucminf(par = par, fn = fn, control = ctrl))
  }
  o <- stats::optim(par, fn, method = "BFGS",
                    control = list(fnscale = 1, maxit = control$maxeval %||% 3000))
  list(par = o$par, value = o$value, convergence = o$convergence,
       info = list(maxgradient = NA_real_))
}

#' Per-species boundary diagnostics at a fitted point
#'
#' Flags niche parameters that sit on an identifiability boundary of the J&S
#' presence-only likelihood: \code{at_floor} (\eqn{\sigma} pinned at
#' \code{sigma_floor}, the degenerate point-mass mode), \code{rho_at_cap}
#' (sliver-like ellipse, \eqn{|\rho| \to \rho_{\mathrm{cap}}}), and
#' \code{sigma_huge} (nearly-uniform niche over the accessible area). These
#' are \emph{diagnostics}, not errors: a species at the floor with few
#' environmental cells is microendemic -- a real identifiability limit --
#' whereas one with many cells and good background coverage was likely a
#' sampling artefact of the background draw.
#'
#' @param data A list from \code{\link{prepare_phylo_niche_data}()}.
#' @param theta A \code{theta} vector (see \code{\link{niche_theta_layout}()}).
#' @inheritParams niche_logpost
#'
#' @return A \code{data.frame} with per-species \code{sigma1}, \code{sigma2},
#'   \code{rho}, empirical SD floor (\code{emp_sd_min}), \code{n_occ},
#'   \code{n_cells}, \code{masked}, and a \code{flag} column
#'   (\code{"ok"}, \code{"masked"}, \code{"at_floor"}, \code{"rho_at_cap"},
#'   \code{"sigma_huge"}, \code{"few_cells"}).
#' @export
niche_boundary_report <- function(data, theta,
                                  parameterization = "noncentered_bounded",
                                  sigma_floor = 0.02, rho_cap = 0.98) {
  model <- switch(parameterization,
                  noncentered = "unb", noncentered_bounded = "bnd",
                  centered = "ctr", centered_lograte = "ctrlr",
                  parameterization)
  pr <- .niche_unpack(theta, data$S, data$P, model)
  S <- data$S; P <- data$P

  if (model %in% c("unb", "bnd")) {
    L_C <- t(chol(data$C))
    log_sigma <- matrix(NA_real_, S, P)
    z_rho <- numeric(S)
    for (k in seq_len(P)) {
      log_sigma[, k] <- pr$log_sigma_anc[k] + sqrt(pr$rate_ls[k]) *
                        (L_C %*% pr$log_sigma_raw[, k])
    }
    z_rho <- drop(pr$z_rho_anc + sqrt(pr$rate_rho) *
                  (L_C %*% pr$z_rho_raw))
  } else {
    log_sigma <- pr$log_sigma
    z_rho <- pr$z_rho
  }
  cap <- if (model == "unb") 1 else rho_cap
  sigma <- exp(log_sigma)
  rho <- cap * tanh(z_rho)
  emp <- niche_empirical(data)
  emp_sd_min <- vapply(seq_len(S), function(s) {
    v <- emp$sd[s, ]
    v <- v[is.finite(v)]
    if (length(v) == 0) NA_real_ else min(v)
  }, numeric(1))

  flag <- rep("ok", S)
  flag[data$masked] <- "masked"
  huge <- apply(sigma, 1, max) > 5 * pmax(apply(data$occ, 2, stats::sd), 1e-3)[1]
  atfloor <- apply(sigma, 1, min) <= sigma_floor * 1.05
  rhocap <- abs(rho) >= cap * 0.98
  few <- !is.na(emp$n_cells) & emp$n_cells < 15 &
         apply(sigma, 1, min) < pmax(emp_sd_min, sigma_floor) * 1.5
  flag[!data$masked & atfloor] <- "at_floor"
  flag[!data$masked & !atfloor & rhocap] <- "rho_at_cap"
  flag[!data$masked & !atfloor & !rhocap & huge] <- "sigma_huge"
  flag[!data$masked & flag == "ok" & few] <- "few_cells"

  data.frame(species = data$species,
             sigma1 = sigma[, 1], sigma2 = sigma[, min(2, P)],
             rho = rho, emp_sd_min = emp_sd_min,
             n_occ = data$N_occ, n_cells = emp$n_cells,
             masked = data$masked, flag = flag,
             stringsAsFactors = FALSE)
}

#' @export
print.niche_multistart <- function(x, ...) {
  cat("niche_multistart (", x$parameterization, ")\n", sep = "")
  cat("  species:", length(x$species),
      "| masked:", sum(x$masked), "\n")
  cat("  map lp:", signif(x$map_lp, 6), "| solutions:",
      length(x$solutions), "\n")
  flags <- x$boundary$flag
  if (any(flags != "ok")) {
    cat("  boundary:", paste(x$boundary$species[flags != "ok"],
                              flags[flags != "ok"], sep = ":",
                              collapse = ", "), "\n")
  }
  invisible(x)
}
