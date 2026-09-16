// niche_bm_constant_sigma.stan
// arfunEvolution: Brownian motion (BM) niche evolution with a constant Sigma.
//
// The niche centroid mu[s] evolves along the phylogeny under Brownian motion
// (BM). The niche shape -- tolerances (log sigma) and correlation (rho, on the
// Fisher-z scale) -- is shared by all species: a single, constant covariance
// matrix Sigma. This is the simplest phylogenetic model and corresponds to the
// constant-Sigma case of Sharma et al. (2024).
//
// Per-species likelihoods use the unweighted background correction of
// Jimenez & Soberon (2019):
//   log L_s = -0.5 * sum_q_occ_s - N_occ[s] * log_sum_exp(a_s)
//
// The correlation is parameterised on the Fisher-z scale z_rho = atanh(rho)
// and back-transformed via tanh, keeping rho in (-1, 1). Designed for the
// bivariate niche ellipsoid (P = 2) produced by the 'arfun' package.
//
// reduce_sum parallelises the per-species likelihood across cores when the
// model is compiled with stan_threads = TRUE. grainsize controls partitioning.

functions {
  // Lower-triangular Cholesky factor of the niche covariance for a bivariate
  // (P = 2) ellipsoid with standard deviations sigma and correlation rho:
  //   Sigma = diag(sigma) * R(rho) * diag(sigma),  R = [[1, rho], [rho, 1]]
  matrix cov_chol_from_rho(vector sigma, real rho) {
    int P = num_elements(sigma);
    matrix[P, P] L_cov = rep_matrix(0.0, P, P);
    L_cov[1, 1] = sigma[1];
    L_cov[1, 2] = 0.0;
    L_cov[2, 1] = rho * sigma[2];
    L_cov[2, 2] = sigma[2] * sqrt(1.0 - square(rho));
    return L_cov;
  }

  // Log-likelihood of a single species s under the unweighted
  // Jimenez & Soberon (2019) model with centroid mu_s and shared shape.
  real species_log_lik(int s,
                       matrix occ, array[] int occ_start, array[] int N_occ,
                       matrix env_m, array[] int m_start, array[] int N_m,
                       vector mu_s, vector sigma, real rho) {
    int P = num_elements(sigma);
    matrix[P, P] L_cov = cov_chol_from_rho(sigma, rho);

    int s_start = occ_start[s];
    int s_end = s_start + N_occ[s] - 1;
    int b_start = m_start[s];
    int b_end = b_start + N_m[s] - 1;

    // Occurrence part: sum of squared Mahalanobis distances
    matrix[N_occ[s], P] occ_s = occ[s_start:s_end, ];
    matrix[P, N_occ[s]] diff_occ = occ_s' - rep_matrix(mu_s, N_occ[s]);
    matrix[P, N_occ[s]] y_occ = mdivide_left_tri_low(L_cov, diff_occ);
    real sum_q_occ = sum(columns_dot_self(y_occ));

    // Background part: log_sum_exp of -0.5 * q_j
    matrix[N_m[s], P] env_m_s = env_m[b_start:b_end, ];
    matrix[P, N_m[s]] diff_m = env_m_s' - rep_matrix(mu_s, N_m[s]);
    matrix[P, N_m[s]] y_m = mdivide_left_tri_low(L_cov, diff_m);
    vector[N_m[s]] a = (-0.5 * columns_dot_self(y_m))';

    return -0.5 * sum_q_occ - N_occ[s] * log_sum_exp(a);
  }

  // reduce_sum partial function: sums the per-species log-likelihood over the
  // slice of species indices in slice_s.
  real partial_sum_bm_constant(array[] int slice_s, int start, int end,
                               matrix occ, array[] int occ_start,
                               array[] int N_occ,
                               matrix env_m, array[] int m_start,
                               array[] int N_m,
                               matrix mu, vector sigma, real rho) {
    real lp = 0.0;
    for (i in 1:size(slice_s)) {
      int s = slice_s[i];
      if (N_occ[s] > 0) lp += species_log_lik(s, occ, occ_start, N_occ, env_m, m_start, N_m,
                            mu[s]', sigma, rho);
    }
    return lp;
  }
}

data {
  int<lower=1> S;                    // number of species
  int<lower=1> P;                    // environmental dimensions (P = 2)

  // Concatenated occurrence data (all species stacked row-wise)
  int<lower=0> N_occ_total;
  matrix[N_occ_total, P] occ;
  array[S] int<lower=1> occ_start;   // 1-based start row index per species
  array[S] int<lower=0> N_occ;       // occurrence count per species

  // Concatenated background data
  int<lower=0> N_m_total;
  matrix[N_m_total, P] env_m;
  array[S] int<lower=1> m_start;
  array[S] int<lower=0> N_m;

  // Phylogenetic variance-covariance matrix (from ape::vcv(tree))
  matrix[S, S] C;

  // Prior hyperparameters
  vector[P] mu_anc_prior;             // prior mean for ancestral centroid
  vector<lower=0>[P] mu_anc_sigma;    // prior sd for ancestral centroid
  real<lower=0> rate_mu_scale;        // half-normal scale for BM rate on mu
  real z_rho_prior;                   // prior mean for shared z_rho
  real<lower=0> z_rho_sigma;          // prior sd for shared z_rho

  int<lower=1> grainsize;             // reduce_sum grainsize
}

transformed data {
  // Pre-compute the Cholesky of the phylogenetic VCV (once, not per iteration)
  cholesky_factor_cov[S] L_C = cholesky_decompose(C);
  array[S] int species_idx;
  for (s in 1:S) species_idx[s] = s;
}

parameters {
  // Ancestral niche centroid
  vector[P] mu_anc;

  // Shared (constant) niche shape
  vector[P] log_sigma;                // shared log tolerances
  real z_rho;                         // shared Fisher-z correlation

  // BM evolutionary rate for centroids (per dimension)
  vector<lower=0>[P] rate_mu;

  // Non-centered centroid deviations (iid standard normal)
  matrix[S, P] mu_raw;
}

transformed parameters {
  // Per-species centroids via non-centered BM:
  //   mu[, k] = mu_anc[k] + sqrt(rate_mu[k]) * L_C * mu_raw[, k]
  matrix[S, P] mu;
  // Constant niche shape (shared across species)
  vector[P] sigma = exp(log_sigma);
  real rho = tanh(z_rho);

  for (k in 1:P)
    mu[, k] = mu_anc[k] + sqrt(rate_mu[k]) * (L_C * mu_raw[, k]);
}

model {
  // --- Priors on ancestral and shared shape parameters ---
  mu_anc ~ normal(mu_anc_prior, mu_anc_sigma);
  log_sigma ~ normal(0, 3);
  z_rho ~ normal(z_rho_prior, z_rho_sigma);

  // --- Prior on BM rate (half-normal) ---
  rate_mu ~ normal(0, rate_mu_scale);

  // --- Non-centered BM deviations: iid N(0,1) ---
  for (k in 1:P)
    mu_raw[, k] ~ std_normal();

  // --- Per-species likelihood (parallelised with reduce_sum) ---
  target += reduce_sum(partial_sum_bm_constant, species_idx, grainsize,
                       occ, occ_start, N_occ, env_m, m_start, N_m,
                       mu, sigma, rho);
}

generated quantities {
  // Per-species covariance and correlation matrices (constant across species
  // here, but kept as arrays for a uniform interface) and per-species
  // log-likelihood for model comparison (LOO/WAIC).
  array[S] matrix[P, P] Sigma;
  array[S] matrix[P, P] R_corr;
  vector[S] log_lik;

  {
    matrix[P, P] L_cov = cov_chol_from_rho(sigma, rho);
    matrix[P, P] R = rep_matrix(rho, P, P);
    for (i in 1:P) R[i, i] = 1.0;

    for (s in 1:S) {
      Sigma[s] = multiply_lower_tri_self_transpose(L_cov);
      R_corr[s] = R;
      if (N_occ[s] > 0) {
        log_lik[s] = species_log_lik(s, occ, occ_start, N_occ,
                                   env_m, m_start, N_m,
                                   mu[s]', sigma, rho);
      } else {
        log_lik[s] = 0;
      }
    }
  }
}
