#' Project ancestral niche states to geographic space
#'
#' Takes ancestral niche parameter estimates from a phylogenetic evolutionary
#' model and projects them onto a geographic raster using the xnicher
#' habitat-suitability engine.
#'
#' @param phy An object of class \code{phylo} from \pkg{ape}, with tip labels
#'   matching the species for which ancestral states were reconstructed.
#' @param anc_mu A \code{S} x \code{P} matrix of posterior mean ancestral
#'   niche centroid positions, one row per node (tips + internal nodes).
#' @param anc_Sigma A \code{S} x \code{P} x \code{P} array of posterior mean
#'   ancestral niche covariance matrices.
#' @param env A multi-layer \code{SpatRaster} from \pkg{terra}, with one
#'   layer per environmental variable, in the same order as the niche
#'   dimensions.
#' @param node_labels Character vector of node labels (tips first, then
#'   internal nodes in the order returned by the model).
#' @param output_dir Character. Directory to write GeoTIFF output files.
#'   If \code{NULL}, returns in-memory rasters.
#' @param overwrite Logical. Whether to overwrite existing output files.
#' @param suffix Character. Suffix appended to output file names (before
#'   extension).
#' @param ... Additional arguments passed to
#'   \code{\link[xnicher]{habitat_suitability}}.
#'
#' @return A named list of SpatRaster objects (one per node), or \code{NULL}
#'   if files were written to disk. Names are node labels.
#'
#' @details
#' The function computes a habitat-suitability raster for each ancestral
#' (or terminal) niche state using the \code{habitat_suitability} function
#' from \pkg{xnicher}. The suitability is the standardized multivariate
#' normal density \eqn{S(x) = exp(-0.5 (x - mu)^\top Sigma^{-1} (x - mu))}
#' evaluated at every grid cell of the environmental raster.
#'
#' For terminal (tip) nodes, the projected niche represents the contemporary
#' niche envelope. For internal nodes, it represents the reconstructed
#' ancestral niche. These can be compared to test hypotheses about niche
#' conservatism vs. divergence.
#'
#' @export
project_ancestral_niche <- function(phy,
                                     anc_mu,
                                     anc_Sigma,
                                     env,
                                     node_labels,
                                     output_dir = NULL,
                                     overwrite = FALSE,
                                     suffix = "ancestral",
                                     ...) {
  if (!requireNamespace("xnicher", quietly = TRUE)) {
    stop(
      "Package 'xnicher' is required but not installed.\n",
      "Install it with:\n",
      '  devtools::install_github("alrobles/xnicher")',
      call. = FALSE
    )
  }

  if (!requireNamespace("terra", quietly = TRUE)) {
    stop(
      "Package 'terra' is required but not installed.",
      call. = FALSE
    )
  }

  if (is.null(output_dir)) {
    output_dir <- tempdir()
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Validate dimensions
  if (nrow(anc_mu) != length(node_labels)) {
    stop("nrow(anc_mu) must equal length(node_labels)")
  }

  if (nrow(anc_mu) != dim(anc_Sigma)[[1]]) {
    stop("dim(anc_Sigma)[1] must equal nrow(anc_mu)")
  }

  # Project each node
  rasters <- vector("list", nrow(anc_mu))
  names(rasters) <- node_labels

  for (i in seq_len(nrow(anc_mu))) {
    mu_i <- as.numeric(anc_mu[i, ])
    Sigma_i <- anc_Sigma[i, , ]

    param <- list(mu = mu_i, Sigma = Sigma_i)

    out_file <- file.path(output_dir,
                          paste0(node_labels[i], "_", suffix, ".tif"))

    if (nzchar(output_dir) && !is.null(output_dir)) {
      rasters[[i]] <- xnicher::habitat_suitability(
        param = param,
        env = env,
        output = out_file,
        overwrite = overwrite,
        ...
      )
    } else {
      rasters[[i]] <- xnicher::habitat_suitability(
        param = param,
        env = env,
        ...
      )
    }
  }

  rasters
}

#' Generate suitability raster for a single niche state
#'
#' Convenience wrapper that projects a single niche (mu, Sigma) onto a
#' geographic raster using xnicher's habitat-suitability engine.
#'
#' @param mu Numeric vector of length \code{P}, niche centroid.
#' @param Sigma Symmetric positive-definite \code{P x P} matrix.
#' @param env A multi-layer \code{SpatRaster} from \pkg{terra}, with one
#'   layer per environmental variable.
#' @param output Character. File path for output GeoTIFF. If \code{""}
#'   (default), returns an in-memory \code{SpatRaster}.
#' @param overwrite Logical. Whether to overwrite existing output file.
#' @param return_log Logical. If \code{FALSE} (default), returns suitability
#'   \eqn{S(x) \in (0, 1]}. If \code{TRUE}, returns log-suitability
#'   \eqn{\log S(x) \le 0}.
#' @param threads Integer. Number of parallel threads for the C++ kernel.
#'
#' @return A one-layer \code{SpatRaster} named \code{"suitability"} (or
#'   \code{"log_suitability"} when \code{return_log = TRUE}).
#'
#' @export
#' @importFrom RcppParallel defaultNumThreads
project_single_niche <- function(mu,
                                   Sigma,
                                   env,
                                   output = "",
                                   overwrite = FALSE,
                                   return_log = FALSE,
                                   threads = NULL) {
  if (!requireNamespace("xnicher", quietly = TRUE)) {
    stop(
      "Package 'xnicher' is required but not installed.\n",
      "Install it with:\n",
      '  devtools::install_github("alrobles/xnicher")',
      call. = FALSE
    )
  }

  # Resolve default thread count
  if (is.null(threads)) {
    threads <- RcppParallel::defaultNumThreads()
  }

  param <- list(mu = mu, Sigma = Sigma)

  xnicher::habitat_suitability(
    param = param,
    env = env,
    output = output,
    overwrite = overwrite,
    return_log = return_log,
    threads = threads
  )
}

#' Compare niche overlap between two ancestral states
#'
#' Computes pairwise niche overlap metrics (Schoener's D, Warren's I) between
#' two projected niche rasters. Useful for testing niche conservatism vs.
#' divergence along phylogenetic branches.
#'
#' @param raster1 A \code{SpatRaster} of suitability values from
#'   \code{project_single_niche()} or \code{project_ancestral_niche()}.
#' @param raster2 A \code{SpatRaster} of suitability values.
#' @param threshold Numeric. Suitability threshold for binary overlap metrics
#'   (default 0.1, corresponding to the fundamental niche boundary).
#' @param method Character. Which metrics to compute:
#'   \code{"schoener"} (Schoener's D),
#'   \code{"warren"} (Warren's I),
#'   \code{"binary"} (binary overlap = proportion of shared suitable area).
#'
#' @return A named numeric vector of overlap metrics.
#'
#' @details
#' Schoener's D measures the normalized sum of absolute differences in
#' suitability between two niches:
#' D = 1 - 0.5 * sum(|S1(x) - S2(x)|) / sum(S1(x) + S2(x))
#'
#' Warren's I measures the Hellinger distance between the two suitability
#' distributions. Both range from 0 (no overlap) to 1 (identical niches).
#'
#' @export
niche_overlap <- function(raster1,
                           raster2,
                           threshold = 0.1,
                           method = c("schoener", "warren", "binary")) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required but not installed.", call. = FALSE)
  }

  method <- match.arg(method, several.ok = TRUE)

  # Extract suitability values (ignore NA)
  s1 <- terra::values(raster1)
  s2 <- terra::values(raster2)

  if (is.null(s1) || is.null(s2)) {
    stop("Could not extract raster values.")
  }

  # Remove cells where either raster is NA
  valid <- stats::complete.cases(s1, s2)
  s1 <- s1[valid]
  s2 <- s2[valid]

  # Remove zero-suitability cells (outside niche boundary)
  nonzero <- (s1 > 0) | (s2 > 0)
  s1 <- s1[nonzero]
  s2 <- s2[nonzero]

  if (length(s1) == 0) {
    stop("No overlapping suitable cells between the two rasters.")
  }

  results <- list()

  if ("schoener" %in% method) {
    # Schoener's D
    diff_sum <- sum(abs(s1 - s2))
    total_sum <- sum(s1 + s2)
    D <- 1 - 0.5 * diff_sum / total_sum
    results$D <- D
  }

  if ("warren" %in% method) {
    # Warren's I: Hellinger distance
    sqrt_s1 <- sqrt(s1 / sum(s1))
    sqrt_s2 <- sqrt(s2 / sum(s2))
    H <- sum((sqrt_s1 - sqrt_s2)^2)
    I <- 1 - H / 2
    results$I <- I
  }

  if ("binary" %in% method) {
    # Binary overlap: fraction of shared suitable area above threshold
    b1 <- as.integer(s1 >= threshold)
    b2 <- as.integer(s2 >= threshold)
    shared <- sum(b1 * b2)
    total1 <- sum(b1)
    total2 <- sum(b2)
    overlap_frac <- shared / pmax(total1, total2)
    results$binary_overlap <- overlap_frac
  }

  unlist(results)
}
