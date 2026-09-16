#' Fit a Bayesian presence-only niche model
#'
#' Compiles and fits a Bayesian fundamental niche model using Stan via
#' \pkg{cmdstanr}. Estimates the fundamental niche as a multivariate normal
#' ellipsoid from occurrence data only.
#'
#' @param env_occ Data frame or matrix of environmental values at presence
#'   points (rows = points, columns = environmental variables).
#' @param eta Numeric, shape parameter for the LKJ prior on the correlation
#'   matrix. Default is 1 (uniform over correlation matrices). Larger values
#'   concentrate mass near the identity matrix.
#' @param sigma_rate Numeric, rate parameter for the exponential prior on
#'   standard deviations. Default is 0.1 (weakly informative).
#' @param model Optional; either \code{NULL} (default, compile on the fly),
#'   a character string with the model name (\code{"presence_only"}), or a
#'   pre-compiled \code{CmdStanModel} object. Passing a pre-compiled model
#'   avoids redundant compilation when running multiple fits in parallel
#'   (e.g., on HPC).
#' @param grainsize Integer, grainsize for \code{reduce_sum} parallelization.
#'   If \code{NULL} or 0, Stan chooses automatically. Used only when the model
#'   is compiled with \code{stan_threads = TRUE}.
#' @param chains Integer, number of MCMC chains (default 4).
#' @param iter_warmup Integer, number of warmup iterations per chain
#'   (default 1000).
#' @param iter_sampling Integer, number of sampling iterations per chain
#'   (default 1000).
#' @param seed Integer, random seed for reproducibility.
#' @param ... Additional arguments passed to \code{cmdstanr::CmdStanModel$sample()}.
#'
#' @return A \code{CmdStanMCMC} object with the fitted model.
#'
#' @details
#' Requires the \pkg{cmdstanr} package and a working CmdStan installation.
#' Install cmdstanr with:
#' \preformatted{
#' install.packages("cmdstanr", repos = c(
#'   "https://stan-dev.r-universe.dev", getOption("repos")
#' ))
#' cmdstanr::install_cmdstan()
#' }
#'
#' For parallel HPC workflows, pre-compile the model once and pass the
#' \code{CmdStanModel} object to avoid concurrent compilation:
#' \preformatted{
#' # Single compilation (e.g., in a setup job) with thread support
#' mod <- compile_model(stan_threads = TRUE)
#'
#' # Run with threading (e.g., 4 threads per chain)
#' fit_niche(occ, model = mod, grainsize = 10, chains = 2,
#'           threads = 4)  # passed via ... to mod$sample()
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate 2D environmental data
#' set.seed(42)
#' occ <- data.frame(
#'   bio01 = rnorm(50, 20, 3),   # annual mean temperature (C)
#'   bio12 = rnorm(50, 800, 100) # annual precipitation (mm)
#' )
#'
#' # Fit the presence-only model (compiles on the fly)
#' fit <- fit_niche(occ, chains = 2)
#'
#' # Summarize posterior
#' fit$summary()
#'
#' # Extract parameters
#' pars <- extract_pars(fit)
#' pars$mu     # centroid (optimal environment)
#' pars$sigma  # tolerances (breadth per dimension)
#' pars$R_corr # correlation (trade-offs)
#'
#' # 95% most suitable region (Mahalanobis ellipsoid) in niche space:
#' # points x such that (x - mu)' Sigma^-1 (x - mu) <= qchisq(0.95, df=2)
#'
#' # ---- Parallel HPC workflow ----
#' # Compile once with thread support
#' mod <- compile_model(stan_threads = TRUE)
#'
#' # Run parallel fit (grainsize = 10, 4 threads per chain)
#' fit_par <- fit_niche(occ, model = mod, grainsize = 10,
#'                      chains = 2, threads = 4)
#' }
#'
#' @export
fit_niche <- function(env_occ,
                      eta = 1,
                      sigma_rate = 0.1,
                      model = NULL,
                      grainsize = NULL,
                      chains = 4,
                      iter_warmup = 1000,
                      iter_sampling = 1000,
                      seed = NULL,
                      ...) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required but not installed.\n",
      "Install it with:\n",
      '  install.packages("cmdstanr", ',
      'repos = c("https://stan-dev.r-universe.dev", getOption("repos")))\n',
      "  cmdstanr::install_cmdstan()",
      call. = FALSE
    )
  }

  # Prepare data for Stan
  stan_data <- prepare_niche_data(env_occ, eta = eta, sigma_rate = sigma_rate)

  # Add grainsize to Stan data if provided
  if (!is.null(grainsize) && grainsize > 0) {
    stan_data$grainsize <- as.integer(grainsize)
  } else {
    stan_data$grainsize <- 0L  # 0 = auto/no threading
  }

  # Get or use pre-compiled model
  if (!is.null(model)) {
    if (inherits(model, "CmdStanModel")) {
      # User passed a pre-compiled model object
      mod <- model
    } else if (is.character(model) && length(model) == 1) {
      # User passed a model name — compile it
      mod <- cmdstanr::cmdstan_model(model_path(model))
    } else {
      stop("model must be NULL, a character string, or a CmdStanModel object")
    }
  } else {
    # Compile on the fly (default)
    mod <- cmdstanr::cmdstan_model(model_path("presence_only"))
  }

  # Run MCMC sampling
  mod$sample(
    data = stan_data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    ...
  )
}

#' Compile a Stan model for the fundamental niche
#'
#' Compiles one of the bundled Stan models and returns a \code{CmdStanModel}
#' object. Use this to pre-compile a model before running multiple fits in
#' parallel (e.g., on HPC).
#'
#' @param model Character string specifying the model to compile. One of:
#'   \code{"presence_only"} (default), \code{"bm_constant"}, \code{"bm_evolving"},
#'   \code{"ou_constant"}, \code{"ou_evolving"}, or \code{"null"}.
#' @param stan_threads Logical, whether to compile with threading support
#'   for \code{reduce_sum} parallelisation (default \code{FALSE}).
#' @param ... Additional arguments passed to \code{cmdstanr::cmdstan_model()}.
#'
#' @return A \code{CmdStanModel} object.
#'
#' @examples
#' \dontrun{
#' mod <- compile_model("presence_only")
#' mod$exe_file()
#'
#' mod_threads <- compile_model("bm_constant", stan_threads = TRUE)
#' }
#'
#' @export
compile_model <- function(model = "presence_only", stan_threads = FALSE, ...) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required but not installed.\n",
      "Install it with:\n",
      '  install.packages("cmdstanr", ',
      'repos = c("https://stan-dev.r-universe.dev", getOption("repos")))\n',
      "  cmdstanr::install_cmdstan()",
      call. = FALSE
    )
  }

  model_path <- model_path(model)
  cmdstanr::cmdstan_model(
    model_path,
    cpp_options = list(stan_threads = stan_threads),
    ...
  )
}

#' Get path to a Stan model file
#'
#' Returns the full path to a Stan model file included in the package.
#'
#' @param model Character string specifying the model. Currently only
#'   \code{"presence_only"} is supported in arfun.
#'
#' @return A character string with the full path to the Stan model file.
#'
#' @examples
#' model_path("presence_only")
#'
#' @export
model_path <- function(model = c("presence_only", "bm_constant", "bm_evolving",
                                  "ou_constant", "ou_evolving", "null")) {
  model <- match.arg(model)

  # Map model names to actual Stan file names
  fname <- switch(
    model,
    presence_only  = "niche_presence_only.stan",
    bm_constant    = "niche_bm_constant_sigma.stan",
    bm_evolving    = "niche_bm_evolving_sigma.stan",
    ou_constant    = "niche_ou_constant_sigma.stan",
    ou_evolving    = "niche_ou_evolving_sigma.stan",
    null           = "niche_null_model.stan"
  )

  path <- system.file("stan", fname, package = "arfun")
  if (path == "") {
    stop("Stan model file not found: ", fname,
         ". Is the package installed correctly?")
  }
  path
}
