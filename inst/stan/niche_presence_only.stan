// niche_presence_only.stan
// Bayesian presence-only fundamental niche model
//
// Estimates the fundamental niche as a multivariate normal (ellipsoid)
// from occurrence data only, without background correction.
//
// This is the simplest model: a standard multivariate normal likelihood.
// Equivalent to the presence-only model in the 'nicher' package
// (Jimenez et al. 2019), but with full Bayesian inference.
//
// Parallelization: uses reduce_sum to parallelize the loop over occurrence
// points when compiled with stan_threads = TRUE.
//
// Parameters:
//   mu    - niche centroid (optimal environmental conditions)
//   sigma - standard deviations (tolerance breadth per variable)
//   L_corr - Cholesky factor of correlation matrix (trade-offs among variables)
//
// The covariance matrix is decomposed as:
//   Sigma = diag(sigma) * L_corr * L_corr' * diag(sigma)
//
// Priors:
//   mu     ~ normal(mu_prior, sigma_prior)  (data-adaptive)
//   sigma  ~ exponential(sigma_rate)
//   L_corr ~ lkj_corr_cholesky(eta)

functions {
  // Partial sum for reduce_sum: computes the likelihood for a subset of occurrences
  real partial_sum_lpmf(array[] int slice_idx, int start, int end,
                        matrix occ, vector mu, matrix L_cov) {
    real sum_lp = 0;
    for (i in start:end) {
      sum_lp += multi_normal_cholesky_lpdf(occ[i]' | mu, L_cov);
    }
    return sum_lp;
  }
}

data {
  int<lower=1> N_occ;              // number of occurrence points
  int<lower=1> P;                  // number of environmental dimensions
  matrix[N_occ, P] occ;           // occurrence environmental data
  int<lower=0> grainsize;         // grainsize for reduce_sum (0 = auto)

  // Prior hyperparameters
  vector[P] mu_prior;             // prior mean for centroid
  vector<lower=0>[P] sigma_prior; // prior sd for centroid
  real<lower=0> sigma_rate;       // rate for exponential prior on sigma
  real<lower=0> eta;              // LKJ shape parameter (1 = uniform)
}

transformed data {
  array[N_occ] int occ_idx;
  for (i in 1:N_occ) occ_idx[i] = i;
}

parameters {
  vector[P] mu;                          // niche centroid
  vector<lower=0>[P] sigma;              // standard deviations
  cholesky_factor_corr[P] L_corr;        // Cholesky factor of correlation matrix
}

transformed parameters {
  // Covariance Cholesky factor: L_cov = diag(sigma) * L_corr
  // so that Sigma = L_cov * L_cov' = diag(sigma) * R * diag(sigma)
  matrix[P, P] L_cov = diag_pre_multiply(sigma, L_corr);
}

model {
  // Priors
  mu ~ normal(mu_prior, sigma_prior);
  sigma ~ exponential(sigma_rate);
  L_corr ~ lkj_corr_cholesky(eta);

  // Likelihood: parallelized via reduce_sum
  // When grainsize = 0, Stan chooses an appropriate grainsize automatically
  if (grainsize > 0) {
    target += reduce_sum(partial_sum_lpmf, occ_idx, grainsize,
                         occ, mu, L_cov);
  } else {
    for (i in 1:N_occ)
      target += multi_normal_cholesky_lpdf(occ[i]' | mu, L_cov);
  }
}

generated quantities {
  // Full covariance matrix for interpretation
  matrix[P, P] Sigma = multiply_lower_tri_self_transpose(L_cov);
  // Correlation matrix
  matrix[P, P] R_corr = multiply_lower_tri_self_transpose(L_corr);
}
