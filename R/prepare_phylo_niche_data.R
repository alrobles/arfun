#' Prepare joint phylogenetic niche data for evolutionary models
#'
#' Builds the concatenated data list used by the joint niche-evolution Stan
#' models (e.g. \code{bm_evolving}, \code{ou_evolving}) and by
#' \code{\link{niche_logpost}()} / \code{\link{niche_multistart}()}. Unlike
#' \code{\link{prepare_niche_data}()} (single-species) and
#' \code{prepare_phylo_data()} (two-stage niche_params interface), this
#' function takes per-species occurrence and background environmental
#' matrices directly and returns the object consumed by the joint models.
#'
#' @param env_occ_list Named list of numeric matrices (one per species):
#'   environmental values at occurrence points. Names must match the tip
#'   labels of \code{tree}.
#' @param env_m_list Named list of numeric matrices with environmental values
#'   at background (accessible-area, \eqn{M}) points for each species. For
#'   commensurate niche estimates across species, pass the \emph{same} shared
#'   sample to every species (see Details).
#' @param tree A \code{phylo} object covering all species in
#'   \code{env_occ_list}.
#' @param mask_species Optional character vector of tip labels to keep in the
#'   phylogenetic prior but \emph{exclude from the likelihood}. Masked
#'   species contribute no occurrence/background rows, so the evolutionary
#'   model predicts their niche from their relatives -- the recommended
#'   treatment for species whose realised occurrences occupy too few
#'   environmental cells to identify a two-dimensional ellipse
#'   (microendemics; see Details).
#' @param eta Reserved for LKJ priors in higher-dimensional extensions.
#' @param rate_mu_scale Scale of the half-normal prior on the BM rate of the
#'   centroid.
#' @param rate_log_sigma_scale Scale of the half-normal prior on the BM rate
#'   of the log tolerances.
#' @param rate_rho_scale Scale of the half-normal prior on the BM rate of
#'   \code{z_rho}. If \code{NULL} (default), set to \code{1 / max(diag(C))}.
#' @param z_rho_anc_prior Prior mean of the ancestral Fisher-z correlation.
#' @param z_rho_anc_sigma Prior SD of the ancestral Fisher-z correlation.
#'
#' @return A named list with fields \code{S, P, N_occ_total, occ, occ_start,
#'   N_occ, N_m_total, env_m, m_start, N_m, C, mu_anc_prior, mu_anc_sigma,
#'   rate_mu_scale, rate_log_sigma_scale, rate_rho_scale, z_rho_anc_prior,
#'   z_rho_anc_sigma, eta, species, masked}. \code{masked} is a logical vector
#'   (length \code{S}, tree-tip order) flagging likelihood-excluded species.
#'
#' @details
#' \strong{Shared background.} The presence-only likelihood normalises each
#' species' ellipse against its background sample. If every species receives
#' an \emph{independent} subsample of a shared accessible area, Monte Carlo
#' holes near small-range clusters create a spurious \eqn{\sigma \to 0}
#' likelihood mode (the "bg-hole" artefact). Always pass the identical
#' background matrix (same rows) for every species, or use
#' \code{\link{shared_background_sample}()} to draw it once.
#'
#' \strong{Masking microendemics.} Presence-only data only bound the
#' fundamental niche from above; a species whose occurrences fall in
#' \eqn{\lesssim 10} unique environmental cells carries essentially no
#' information about niche shape, and its likelihood is dominated by the
#' \eqn{\sigma \to 0} boundary. \code{mask_species} removes such species
#' from the likelihood while keeping them in the phylogenetic covariance, so
#' their niche is predicted from their relatives -- turning an
#' identifiability limit into a testable prediction target.
#'
#' @seealso \code{\link{niche_logpost}}, \code{\link{niche_multistart}},
#'   \code{\link{niche_boundary_report}}
#' @export
prepare_phylo_niche_data <- function(env_occ_list,
                                    env_m_list,
                                    tree,
                                    mask_species = NULL,
                                    eta = 1,
                                    rate_mu_scale = 1,
                                    rate_log_sigma_scale = 1,
                                    rate_rho_scale = NULL,
                                    z_rho_anc_prior = 0,
                                    z_rho_anc_sigma = 1) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required for phylogenetic data preparation.\n",
         "Install it with: install.packages('ape')", call. = FALSE)
  }
  if (is.null(names(env_occ_list)) || is.null(names(env_m_list))) {
    stop("'env_occ_list' and 'env_m_list' must be named lists, ",
         "with names matching the tip labels of 'tree'.")
  }
  sp_names <- names(env_occ_list)
  tip_labels <- tree$tip.label
  missing_tips <- setdiff(sp_names, tip_labels)
  if (length(missing_tips) > 0) {
    stop("Species in 'env_occ_list' not found in tree tip labels: ",
         paste(missing_tips, collapse = ", "))
  }
  missing_sp <- setdiff(tip_labels, sp_names)
  if (length(missing_sp) > 0) {
    stop("Tree tips not found in 'env_occ_list': ", paste(missing_sp,
         collapse = ", "))
  }
  if (!identical(sort(names(env_m_list)), sort(sp_names))) {
    stop("'env_m_list' must have the same species names as 'env_occ_list'")
  }

  if (!is.null(mask_species)) {
    if (!is.character(mask_species)) {
      stop("'mask_species' must be a character vector of tip labels")
    }
    bad <- setdiff(mask_species, tip_labels)
    if (length(bad) > 0) {
      stop("'mask_species' not in tree: ", paste(bad, collapse = ", "))
    }
  }
  masked <- tip_labels %in% (mask_species %||% character(0))

  env_occ_list <- env_occ_list[tip_labels]
  env_m_list <- env_m_list[tip_labels]
  S <- length(tip_labels)
  occ_mats <- lapply(env_occ_list, as.matrix)
  m_mats <- lapply(env_m_list, as.matrix)
  P <- ncol(occ_mats[[1]])
  for (s in seq_len(S)) {
    sp <- tip_labels[s]
    if (ncol(occ_mats[[s]]) != P) {
      stop("Inconsistent number of columns in 'env_occ_list' for species: ",
           sp)
    }
    if (ncol(m_mats[[s]]) != P) {
      stop("Number of columns in 'env_m_list' does not match 'env_occ_list' ",
           "for species: ", sp)
    }
    if (!masked[s]) {
      if (nrow(occ_mats[[s]]) < 2) {
        stop("Need at least 2 occurrence points for species: ", sp,
             " (or list it in 'mask_species')")
      }
      if (nrow(m_mats[[s]]) < 2) {
        stop("Need at least 2 background points for species: ", sp)
      }
    }
  }

  # Masked species contribute zero rows to the likelihood arrays; they remain
  # in C, so the phylogenetic prior still informs their niche parameters.
  if (any(masked)) {
    occ_mats[masked] <- lapply(occ_mats[masked],
                               function(m) m[0, , drop = FALSE])
    m_mats[masked] <- lapply(m_mats[masked],
                             function(m) m[0, , drop = FALSE])
  }

  N_occ <- vapply(occ_mats, nrow, integer(1))
  N_m <- vapply(m_mats, nrow, integer(1))
  occ_start <- as.integer(c(1L, cumsum(N_occ[-S]) + 1L))
  m_start <- as.integer(c(1L, cumsum(N_m[-S]) + 1L))
  occ_all <- do.call(rbind, occ_mats)
  env_m_all <- do.call(rbind, m_mats)

  C <- ape::vcv(tree)
  C <- C[tip_labels, tip_labels]

  # Data-adaptive ancestral priors from the realised occurrences.
  mu_anc_prior <- colMeans(occ_all)
  mu_anc_sigma <- apply(occ_all, 2, stats::sd) * 10
  mu_anc_sigma <- pmax(mu_anc_sigma, 1e-06)

  t_root <- diag(C)
  t_max <- max(t_root)
  if (is.null(rate_rho_scale)) {
    rate_rho_scale <- if (is.finite(t_max) && t_max > 0) 1 / t_max else 1
  }

  list(
    S = S, P = P,
    N_occ_total = nrow(occ_all), occ = occ_all,
    occ_start = occ_start, N_occ = as.integer(N_occ),
    N_m_total = nrow(env_m_all), env_m = env_m_all,
    m_start = m_start, N_m = as.integer(N_m),
    C = C,
    mu_anc_prior = as.numeric(mu_anc_prior),
    mu_anc_sigma = as.numeric(mu_anc_sigma),
    rate_mu_scale = rate_mu_scale,
    rate_log_sigma_scale = rate_log_sigma_scale,
    rate_rho_scale = rate_rho_scale,
    z_rho_anc_prior = z_rho_anc_prior,
    z_rho_anc_sigma = z_rho_anc_sigma,
    eta = eta,
    species = tip_labels,
    masked = masked
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Draw one shared background sample for a clade
#'
#' Draws a single background (accessible-area) sample and returns it so the
#' identical rows can be passed to every species' \code{env_m_list} entry in
#' \code{\link{prepare_phylo_niche_data}()}. Sharing the sample keeps the
#' availability envelope commensurate across tips and avoids the per-species
#' Monte Carlo holes that create a spurious \eqn{\sigma \to 0} likelihood
#' mode.
#'
#' @param env_matrix Numeric matrix of environmental values over the clade
#'   accessible area (e.g. sampled cells of the biorealm/bioregion).
#' @param n Integer, number of background rows to draw (default: all rows).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A matrix with \code{n} rows sampled without replacement from
#'   \code{env_matrix}.
#' @export
shared_background_sample <- function(env_matrix, n = NULL, seed = NULL) {
  env_matrix <- as.matrix(env_matrix)
  if (is.null(n) || n >= nrow(env_matrix)) return(env_matrix)
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
  env_matrix[sample.int(nrow(env_matrix), n), , drop = FALSE]
}
