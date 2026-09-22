#' Project supplied niche states to geographic space
#'
#' Takes niche parameter estimates (extant species or ancestral states
#' previously reconstructed by a phylogenetic evolutionary model) and
#' projects them onto a geographic raster using the xnicher
#' habitat-suitability engine. This function does not itself perform
#' ancestral reconstruction; it consumes states estimated elsewhere.
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
#' @return A named list of SpatRaster objects (one per node); names are
#'   node labels. When \code{output_dir} is \code{NULL} the rasters are
#'   in memory; otherwise they are file-backed GeoTIFFs.
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

  write_files <- !is.null(output_dir)
  if (write_files && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Validate dimensions
  if (nrow(anc_mu) != length(node_labels)) {
    stop("nrow(anc_mu) must equal length(node_labels)")
  }

  if (nrow(anc_mu) != dim(anc_Sigma)[[1]]) {
    stop("dim(anc_Sigma)[1] must equal nrow(anc_mu)")
  }

  p <- ncol(anc_mu)
  if (!isTRUE(all.equal(dim(anc_Sigma), c(nrow(anc_mu), p, p)))) {
    stop("dim(anc_Sigma) must be (nrow(anc_mu), ncol(anc_mu), ncol(anc_mu))")
  }

  # node_labels must correspond to tip + node labels of phy when available
  if (!is.null(phy)) {
    tree_nodes <- c(phy$tip.label, phy$node.label)
    unknown <- setdiff(node_labels, tree_nodes)
    if (length(phy$node.label) > 0 && length(unknown) > 0) {
      stop("node_labels not found in 'phy': ",
           paste(utils::head(unknown, 5), collapse = ", "))
    }
  }

  # Project each node
  rasters <- vector("list", nrow(anc_mu))
  names(rasters) <- node_labels

  for (i in seq_len(nrow(anc_mu))) {
    mu_i <- as.numeric(anc_mu[i, ])
    Sigma_i <- anc_Sigma[i, , ]

    param <- list(mu = mu_i, Sigma = Sigma_i)

    out_file <- if (write_files)
      file.path(output_dir, paste0(node_labels[i], "_", suffix, ".tif"))
    else ""

    rasters[[i]] <- xnicher::habitat_suitability(
      param = param,
      env = env,
      output = out_file,
      overwrite = overwrite,
      ...
    )
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
#' @param threshold Numeric. Suitability threshold for the binary metric
#'   (default 0.1). This is an arbitrary cutoff on the suitability scale,
#'   not an estimate of a fundamental-niche boundary.
#' @param method Character. Which metrics to compute:
#'   \code{"schoener"} (Schoener's D),
#'   \code{"warren"} (Warren's I),
#'   \code{"binary"} (binary overlap = proportion of shared suitable area).
#'
#' @return A named numeric vector of overlap metrics.
#'
#' @details
#' These are raster-based summaries comparing two projected suitability
#' surfaces. The two rasters must share extent, resolution, and CRS.
#'
#' Schoener's D is computed on the suitability densities normalized to
#' sum 1:
#' D = 1 - 0.5 * sum(|p1(x) - p2(x)|)
#'
#' Warren's I is the Hellinger-based similarity of the same normalized
#' densities. Both range from 0 (no overlap) to 1 (identical niches).
#' These are descriptive comparisons of projected surfaces, not a formal
#' test of niche conservatism.
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

  # rasters must share extent, resolution and alignment
  if (!isTRUE(terra::compareGeom(raster1, raster2, stopOnError = FALSE))) {
    stop("raster1 and raster2 must have identical extent, resolution, ",
         "and CRS (see terra::compareGeom).")
  }

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
    # Schoener's D on normalized densities (each raster scaled to sum 1)
    p1 <- s1 / sum(s1)
    p2 <- s2 / sum(s2)
    results$D <- 1 - 0.5 * sum(abs(p1 - p2))
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
    denom <- pmax(total1, total2)
    results$binary_overlap <- if (denom > 0) shared / denom else 0
  }

  unlist(results)
}
