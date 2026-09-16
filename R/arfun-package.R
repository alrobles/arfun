#' arfun: Ancestral Reconstruction of the Fundamental Niche
#'
#' Bayesian estimation of the fundamental niche as a multivariate normal
#' ellipsoid from occurrence data using Stan via cmdstanr, with phylogenetic
#' evolutionary models (BM, OU) and ancestral state reconstruction.
#'
#' \pkg{arfun} is a complete framework for:
#' \enumerate{
#'   \item Estimating the fundamental niche as a bivariate normal ellipsoid
#'     from occurrence data (single species)
#'   \item Fitting phylogenetic evolutionary models (BM, OU) to the niche
#'     parameters across multiple species
#'   \item Reconstructing ancestral niche states
#'   \item Projecting niches into geographic space
#' }
#'
#' @section Single-species niche estimation:
#'
#' The fundamental niche is estimated as a bivariate normal ellipsoid
#' where occurrence points follow a multivariate normal distribution
#' with mean mu (centroid) and covariance Sigma = diag(sigma) * R * diag(sigma)
#'
#' where mu is the centroid (optimal environment), sigma are the
#' tolerances (niche breadths), and R is the correlation matrix.
#'
#' The package provides:
#' \describe{
#'   \item{prepare_niche_data}{Convert occurrence data into the list
#'     format required by the Stan model.}
#'   \item{fit_niche}{Compile and fit the Bayesian presence-only niche
#'     model using cmdstanr.}
#'   \item{extract_pars}{Extract posterior summaries of the niche centroid
#'     (\code{mu}), standard deviations (\code{sigma}), covariance matrix
#'     (\code{Sigma}) and correlation matrix (\code{R_corr}).}
#' }
#'
#' @section Phylogenetic evolutionary models:
#'
#' For clade-level analyses, \pkg{arfun} provides five evolutionary models:
#'
#' \tabular{ll}{
#'   \strong{Model} \tab \strong{Description} \cr
#'   bm_constant   \tab BM on mu, sigma, rho; constant Sigma
#'     (Sharma et al. 2024 case) \cr
#'   bm_evolving  \tab BM on mu, sigma, rho; evolving Sigma
#'     (per-species covariance) \cr
#'   ou_constant  \tab OU on mu, sigma, rho; constant Sigma \cr
#'   ou_evolving  \tab OU on mu, sigma, rho; evolving Sigma \cr
#'   null         \tab Independent species, no phylogeny; each species has
#'     its own mu_s and Sigma_s
#' }
#'
#' @section Parallelisation:
#'
#' All models support \code{reduce_sum} parallelisation via:
#' \itemize{
#'   \item \code{compile_model(stan_threads = TRUE)} — compile with thread
#'     support
#'   \item \code{fit_niche(..., grainsize = n, threads = m)} — run with
#'     parallel likelihood evaluation
#' }
#'
#' For HPC workflows, pre-compile the model once and pass the \code{CmdStanModel}
#' object to avoid redundant compilation across parallel jobs.
#'
#' @name arfun
#' @aliases arfun-package
#' @author Angel Robles-Fernandez
#' @author Jorge Soberon
#' @useDynLib arfun, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @references
#' Jimenez, L., Soberon, J., Christen, J.A. and Soto, D. (2019).
#' On the problem of modeling a fundamental niche from occurrence data.
#' Ecological Modelling 397: 74-83.
#' DOI: 10.1016/j.ecolmodel.2018.11.013
#'
#' Jimenez, L. and Soberon, J. (2022).
#' Estimating the fundamental niche: accounting for the uneven availability
#' of existing climates in the calibration area.
#' Ecological Modelling 464: 109823.
#' DOI: 10.1016/j.ecolmodel.2021.109823
#'
#' Sharma, S. et al. (2024). Phylogenetic comparative methods for niche
#' evolution. In review.
NULL
