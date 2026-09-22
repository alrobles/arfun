## Internal helpers shared across the niche-evolution objectives.

`%||%` <- function(x, y) if (is.null(x)) y else x

# Canonical parameterization code: public long names map to the short
# codes used internally and by the C++ backend.
.canon_param <- function(parameterization) {
  model <- match.arg(
    parameterization,
    c("noncentered_bounded", "noncentered", "centered",
      "centered_lograte", "unb", "bnd", "ctr", "ctrlr"))
  switch(model,
         noncentered = "unb", noncentered_bounded = "bnd",
         centered = "ctr", centered_lograte = "ctrlr",
         model)
}

# TRUE for the non-centered parameterizations (raw standard-normal
# deviations mapped through chol(C)).
.is_noncentered <- function(model) model %in% c("unb", "bnd")

# Stan-style matrix element names in column-major order:
# "prefix[1,1]", "prefix[2,1]", ..., "prefix[S,P]".
.mat_names <- function(prefix, S, P) {
  as.vector(outer(seq_len(S), seq_len(P),
                  function(s, k) paste0(prefix, "[", s, ",", k, "]")))
}

# Species trait matrices implied by an unpacked parameter vector:
# non-centered codes solve the phylo map, centered codes read the
# matrices directly. Returns list(mu, log_sigma, z_rho).
.reconstruct_traits <- function(pr, data, model) {
  S <- data$S; P <- data$P
  if (!.is_noncentered(model)) {
    return(list(mu = pr$mu, log_sigma = pr$log_sigma, z_rho = pr$z_rho))
  }
  L_C <- t(chol(data$C))
  mu <- matrix(NA_real_, S, P)
  log_sigma <- matrix(NA_real_, S, P)
  for (k in seq_len(P)) {
    mu[, k] <- pr$mu_anc[k] +
      sqrt(pr$rate_mu[k]) * (L_C %*% pr$mu_raw[, k])
    log_sigma[, k] <- pr$log_sigma_anc[k] +
      sqrt(pr$rate_ls[k]) * (L_C %*% pr$log_sigma_raw[, k])
  }
  z_rho <- drop(pr$z_rho_anc +
                sqrt(pr$rate_rho) * (L_C %*% pr$z_rho_raw))
  list(mu = mu, log_sigma = log_sigma, z_rho = z_rho)
}
