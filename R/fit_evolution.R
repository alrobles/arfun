#' Prepare phylogenetic data for evolutionary models
#'
#' Converts a phylogeny, species niche estimates, and environmental data into
#' the list format required by the Stan evolutionary models.
#'
#' @param phy An object of class \code{phylo} from the \pkg{ape} package,
#'   representing the phylogenetic tree of the species.
#' @param species Character vector of species names, in the same order as the
#'   rows of \code{niche_params}.
#' @param niche_params A matrix or data.frame with columns:
#'   \code{mu_1}, \code{mu_2} (centroid dimensions),
#'   \code{sigma_1}, \code{sigma_2} (tolerances),
#'   \code{rho} (correlation), one row per species.
#' @param occ_data A list of matrices, one per species, with environmental
#'   values at occurrence points.
#' @param env_m_data A list of matrices, one per species, with environmental
#'   values at background points.
#' @param mu_anc_prior Numeric vector of length 2, prior mean for ancestral
#'   centroid.
#' @param mu_anc_sigma Numeric vector of length 2, prior SD for ancestral
#'   centroid.
#' @param rate_mu_scale Numeric, scale for half-normal prior on BM rate of
#'   centroid.
#' @param rate_log_sigma_scale Numeric, scale for half-normal prior on BM rate
#'   of log tolerances.
#' @param z_rho_anc_prior Numeric, prior mean for ancestral z_rho (atanh(rho)).
#' @param z_rho_anc_sigma Numeric, prior SD for ancestral z_rho.
#' @param rate_rho_scale Numeric, scale for half-normal prior on BM rate of
#'   z_rho.
#' @param alpha_mu_scale Numeric vector of length 2, scale for prior on OU
#'   alpha (centroid).
#' @param sigma_sq_mu_scale Numeric vector of length 2, scale for prior on OU
#'   sigma^2 (centroid).
#' @param alpha_log_sigma_scale Numeric vector of length 2, scale for prior on
#'   OU alpha (log tolerances).
#' @param sigma_sq_log_sigma_scale Numeric vector of length 2, scale for prior
#'   on OU sigma^2 (log tolerances).
#' @param alpha_rho_scale Numeric, scale for prior on OU alpha (rho).
#' @param sigma_sq_rho_scale Numeric, scale for prior on OU sigma^2 (rho).
#' @param grain_size Integer, grainsize for \code{reduce_sum} parallelisation.
#'
#' @return A named list suitable for passing to \code{fit_evolution()} as the
#'   \code{stan_data} argument.
#'
#' @export
prepare_phylo_data <- function(phy,
                                species,
                                niche_params,
                                occ_data,
                                env_m_data,
                                mu_anc_prior = c(0, 0),
                                mu_anc_sigma = c(5, 5),
                                rate_mu_scale = 1,
                                rate_log_sigma_scale = 0.5,
                                z_rho_anc_prior = 0,
                                z_rho_anc_sigma = 1,
                                rate_rho_scale = 0.5,
                                alpha_mu_scale = c(1, 1),
                                sigma_sq_mu_scale = c(1, 1),
                                alpha_log_sigma_scale = c(1, 1),
                                sigma_sq_log_sigma_scale = c(1, 1),
                                alpha_rho_scale = 1,
                                sigma_sq_rho_scale = 1,
                                grain_size = 1) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required but not installed.",
         call. = FALSE)
  }

  if (!inherits(phy, "phylo")) {
    stop("phy must be an object of class 'phylo'", call. = FALSE)
  }
  if (!is.character(species) || length(species) < 2L ||
      anyNA(species) || anyDuplicated(species)) {
    stop("species must be at least two unique, non-missing names",
         call. = FALSE)
  }
  if (!is.matrix(niche_params) && !is.data.frame(niche_params)) {
    stop("niche_params must be a matrix or data.frame", call. = FALSE)
  }
  if (!is.numeric(mu_anc_prior) || length(mu_anc_prior) != 2L ||
      anyNA(mu_anc_prior) ||
      any(!is.finite(mu_anc_prior)) ||
      !is.numeric(mu_anc_sigma) || length(mu_anc_sigma) != 2L ||
      anyNA(mu_anc_sigma) ||
      any(!is.finite(mu_anc_sigma)) || any(mu_anc_sigma <= 0)) {
    stop("ancestral centroid priors must be finite vectors of length 2",
         call. = FALSE)
  }
  positive_scalar_inputs <- list(
    rate_mu_scale, rate_log_sigma_scale, z_rho_anc_sigma,
    rate_rho_scale, alpha_rho_scale, sigma_sq_rho_scale
  )
  positive_vector_inputs <- list(
    alpha_mu_scale, sigma_sq_mu_scale,
    alpha_log_sigma_scale, sigma_sq_log_sigma_scale
  )
  if (any(!vapply(positive_scalar_inputs, is.numeric, logical(1))) ||
      any(!vapply(positive_vector_inputs, is.numeric, logical(1))) ||
      any(!vapply(positive_scalar_inputs, function(x) {
        length(x) == 1L && is.finite(x) && x > 0
      }, logical(1))) ||
      any(!vapply(positive_vector_inputs, function(x) {
        length(x) == 2L && all(is.finite(x)) && all(x > 0)
      }, logical(1))) ||
      !is.numeric(z_rho_anc_prior) || length(z_rho_anc_prior) != 1L ||
      !is.finite(z_rho_anc_prior) ||
      length(grain_size) != 1L || !is.numeric(grain_size) ||
      !is.finite(grain_size) || grain_size < 1) {
    stop("evolutionary prior scales and grain_size must be positive finite values",
         call. = FALSE)
  }
  if (nrow(niche_params) != length(species)) {
    stop("niche_params must have one row per species", call. = FALSE)
  }
  required_niche_cols <- c("mu_1", "mu_2", "sigma_1", "sigma_2", "rho")
  if (!all(required_niche_cols %in% colnames(niche_params))) {
    stop(
      "niche_params must contain columns: ",
      paste(required_niche_cols, collapse = ", "),
      call. = FALSE
    )
  }
  niche_values <- as.matrix(niche_params[, required_niche_cols, drop = FALSE])
  if (!is.numeric(niche_values) || anyNA(niche_values) ||
      any(!is.finite(niche_values))) {
    stop("niche_params must contain only finite, non-missing values",
         call. = FALSE)
  }
  if (!is.list(occ_data) || length(occ_data) != length(species) ||
      !is.list(env_m_data) || length(env_m_data) != length(species)) {
    stop("occ_data and env_m_data must be lists with one entry per species",
         call. = FALSE)
  }
  valid_table <- function(x) is.matrix(x) || is.data.frame(x)
  if (any(!vapply(occ_data, valid_table, logical(1))) ||
      any(!vapply(env_m_data, valid_table, logical(1)))) {
    stop("occ_data and env_m_data entries must be matrices or data.frames",
         call. = FALSE)
  }
  occ_dims <- vapply(occ_data, ncol, integer(1))
  env_m_dims <- vapply(env_m_data, ncol, integer(1))
  if (any(occ_dims != 2L) || any(env_m_dims != 2L)) {
    stop("evolutionary Stan models currently require exactly 2 dimensions",
         call. = FALSE)
  }
  if (any(vapply(occ_data, nrow, integer(1)) < 1L) ||
      any(vapply(env_m_data, nrow, integer(1)) < 1L)) {
    stop("occ_data and env_m_data entries must contain at least one row",
         call. = FALSE)
  }
  all_tables <- c(occ_data, env_m_data)
  if (any(vapply(all_tables, function(x) {
    x <- as.matrix(x)
    !is.numeric(x) || anyNA(x) || any(!is.finite(x))
  }, logical(1)))) {
    stop("occ_data and env_m_data must contain only finite numeric values",
         call. = FALSE)
  }

  # Match species to tree tips
  if (!all(species %in% phy$tip.label)) {
    missing <- setdiff(species, phy$tip.label)
    stop("Species not found in phylogeny: ", paste(missing, collapse = ", "))
  }

  # Order all species-level inputs by tree tip order
  tip_idx <- match(species, phy$tip.label)
  ordered_idx <- order(tip_idx)
  species_ordered <- species[ordered_idx]
  occ_data <- occ_data[ordered_idx]
  env_m_data <- env_m_data[ordered_idx]
  niche_params <- niche_params[ordered_idx, , drop = FALSE]

  # Compute phylogenetic variance-covariance matrix
  phy <- ape::multi2di(phy)  # ensure dichotomous
  phy <- ape::keep.tip(phy, species_ordered)
  C <- ape::vcv(phy)

  # Concatenate occurrence data
  N_occ_total <- sum(sapply(occ_data, nrow))
  P <- 2L
  occ_list <- lapply(occ_data, as.matrix)
  occ_concat <- do.call(rbind, occ_list)
  occ_start <- cumsum(c(0, sapply(occ_data, nrow)))[-(length(occ_data) + 1)] + 1
  N_occ <- sapply(occ_data, nrow)

  # Concatenate background data
  N_m_total <- sum(sapply(env_m_data, nrow))
  env_m_list <- lapply(env_m_data, as.matrix)
  env_m_concat <- do.call(rbind, env_m_list)
  m_start <- cumsum(c(0, sapply(env_m_data, nrow)))[-(length(env_m_data) + 1)] + 1
  N_m <- sapply(env_m_data, nrow)

  # Extract niche parameters (already ordered)
  mu <- niche_params[, c("mu_1", "mu_2"), drop = FALSE]
  sigma <- niche_params[, c("sigma_1", "sigma_2"), drop = FALSE]
  rho <- niche_params[, "rho"]

  list(
    S = length(species_ordered),
    P = P,
    N_occ_total = N_occ_total,
    occ = occ_concat,
    occ_start = occ_start,
    N_occ = N_occ,
    N_m_total = N_m_total,
    env_m = env_m_concat,
    m_start = m_start,
    N_m = N_m,
    C = C,
    T_anc = C,
    mu_prior = colMeans(as.matrix(niche_params[, c("mu_1", "mu_2"),
                                                drop = FALSE])),
    mu_sigma = rep(5, P),
    z_rho_prior = 0,
    z_rho_sigma = 1,
    mu_anc_prior = mu_anc_prior,
    mu_anc_sigma = mu_anc_sigma,
    rate_mu_scale = rate_mu_scale,
    rate_log_sigma_scale = rate_log_sigma_scale,
    z_rho_anc_prior = z_rho_anc_prior,
    z_rho_anc_sigma = z_rho_anc_sigma,
    rate_rho_scale = rate_rho_scale,
    alpha_mu_scale = alpha_mu_scale,
    sigma_sq_mu_scale = sigma_sq_mu_scale,
    alpha_log_sigma_scale = alpha_log_sigma_scale,
    sigma_sq_log_sigma_scale = sigma_sq_log_sigma_scale,
    alpha_rho_scale = alpha_rho_scale,
    sigma_sq_rho_scale = sigma_sq_rho_scale,
    grainsize = as.integer(grain_size)
  )
}

#' Fit an evolutionary model for niche parameters
#'
#' Unified interface for fitting phylogenetic evolutionary models (BM, OU,
#' null) to the fundamental niche parameters estimated by \code{fit_niche()}.
#'
#' @param stan_data A named list of Stan data, as returned by
#'   \code{prepare_phylo_data()} or constructed manually. For the null model,
#'   the phylogenetic matrix \code{C} should be omitted.
#' @param model Character string specifying the model. One of:
#'   \code{"bm_constant"} (BM with constant Sigma),
#'   \code{"bm_evolving"} (BM with evolving Sigma),
#'   \code{"ou_constant"} (OU with constant Sigma),
#'   \code{"ou_evolving"} (OU with evolving Sigma),
#'   \code{"null"} (independent species, no phylogeny).
#' @param chains Integer, number of MCMC chains (default 4).
#' @param iter_warmup Integer, warmup iterations per chain (default 1000).
#' @param iter_sampling Integer, sampling iterations per chain (default 1000).
#' @param seed Integer, random seed for reproducibility.
#' @param grainsize Integer, grainsize for \code{reduce_sum}.
#' @param threads_per_chain Integer, threads per chain for within-chain
#'   parallelisation (default 1). Values > 1 require the model compiled with
#'   \code{stan_threads = TRUE}.
#' @param mod Optional pre-compiled \code{CmdStanModel}. If \code{NULL}, the
#'   model is compiled internally via \code{compile_model()}.
#' @param ... Additional arguments passed to \code{mod$sample()}.
#'
#' @return A \code{CmdStanMCMC} object. Use \code{$summary()} to obtain
#'   posterior summaries. Generated quantities include \code{Sigma},
#'   \code{R_corr}, and \code{log_lik} per species.
#'
#' @details
#' The model names map to the following Stan files:
#' \tabular{ll}{
#'   \code{bm_constant}   \tab \code{inst/stan/niche_bm_constant_sigma.stan} \cr
#'   \code{bm_evolving}  \tab \code{inst/stan/niche_bm_evolving_sigma.stan} \cr
#'   \code{ou_constant}  \tab \code{inst/stan/niche_ou_constant_sigma.stan} \cr
#'   \code{ou_evolving}  \tab \code{inst/stan/niche_ou_evolving_sigma.stan} \cr
#'   \code{null}         \tab \code{inst/stan/niche_null_model.stan}
#' }
#'
#' For HPC workflows, pre-compile the model with \code{compile_model()} and
#' pass the returned \code{CmdStanModel} object to the \code{mod} argument to
#' avoid redundant compilation across parallel jobs.
#'
#' @examples
#' \dontrun{
#'   # Assume stan_data is prepared via prepare_phylo_data()
#'   mod <- compile_model("bm_constant")
#'   fit <- fit_evolution(stan_data, model = "bm_constant", mod = mod, chains = 2)
#'   fit$summary()
#' }
#'
#' @export
fit_evolution <- function(stan_data,
                          model = c("bm_constant", "bm_evolving",
                                    "ou_constant", "ou_evolving", "null"),
                          chains = 4,
                          iter_warmup = 1000,
                          iter_sampling = 1000,
                          seed = NULL,
                          grainsize = 1,
                          threads_per_chain = 1,
                          mod = NULL,
                          ...) {
  model <- match.arg(model)

  if (is.null(stan_data$grainsize)) {
    stan_data$grainsize <- as.integer(grainsize)
  }

  if (is.null(mod)) {
    mod <- compile_model(model)
  }

  # Pass the setting explicitly when the model was compiled with threads.
  sample_args <- list(
    data = stan_data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    threads_per_chain = as.integer(threads_per_chain)
  )

  do.call(mod$sample, c(sample_args, list(...)))
}

#' Extract evolutionary parameters from a fitted model
#'
#' Extracts posterior summaries of phylogenetic evolutionary parameters
#' (ancestral states, evolutionary rates, OU parameters) from a fitted
#' \code{fit_evolution()} object.
#'
#' @param fit A \code{CmdStanMCMC} object returned by \code{fit_evolution()}.
#' @param variables Character vector of variable names to extract. If
#'   \code{NULL}, extracts all parameters.
#'
#' @return A data.frame with posterior summaries (mean, sd, 2.5\%, 50\%, 97.5\%,
#'   effective sample size, R-hat) for each requested variable.
#'
#' @export
extract_evolution_pars <- function(fit,
                                    variables = NULL) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required but not installed.",
         call. = FALSE)
  }

  # Get all parameter names from the fit
  all_params <- colnames(posterior::rhat(fit))

  if (is.null(variables)) {
    variables <- all_params
  }

  # Filter to requested variables
  available <- intersect(variables, all_params)
  if (length(available) == 0) {
    stop("None of the requested variables found in the fitted model.",
         call. = FALSE)
  }

  # Extract and summarise
  draws <- posterior::as_draws_rvars(fit)
  sum <- posterior::summarise_draws(draws[available],
                                    "mean", "sd", "median",
                                    "quantile", name = "q5",
                                    lwr = 0.025, upr = 0.975,
                                    "ess_bulk", "rhat")

  as.data.frame(sum)
}
