// niche_bm_evolving_sigma.stan
// arfunEvolution: Brownian motion (BM) niche evolution with evolving Sigma.
//
// Niche centroids mu[s], log tolerances log_sigma[s], and the correlation
// rho[s] (on the Fisher-z scale z_rho = atanh(rho)) all evolve along the
// phylogeny under independent Brownian motion (BM) processes. Each species
// therefore has its own covariance matrix Sigma[s], so the niche shape itself
// is an evolutionary trait. This is the full niche-evolution model.
//
// The correlation is kept in (-1, 1) by modelling the unconstrained
// z_rho[s] = atanh(rho[s]) as the evolving trait and back-transforming with
// tanh. Designed for the bivariate niche ellipsoid (P = 2) of the 'arfun'
// package.
//
// Per-species likelihoods use the unweighted background correction of
// Jimenez & Soberon (2019). Non-centered parameterization is used for
// efficient HMC sampling. reduce_sum parallelises the per-species likelihood.

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
  // Jimenez & Soberon (2019) model with per-species centroid and shape.
  real species_log_lik(int s,
                       matrix occ, array[] int occ_start, array[] int N_occ,
                       matrix env_m, array[] int m_start, array[] int N_m,
                       vector mu_s, vector sigma_s, real rho_s) {
    int P = num_elements(sigma_s);
    matrix[P, P] L_cov = cov_chol_from_rho(sigma_s, rho_s);

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

  // reduce_sum partial function over a slice of species.
  real partial_sum_bm_evolving(array[] int slice_s, int start, int end,
                               matrix occ, array[] int occ_start,
                               array[] int N_occ,
                               matrix env_m, array[] int m_start,
                               array[] int N_m,
                               matrix mu, matrix log_sigma, vector rho) {
    real lp = 0.0;
    for (i in 1:size(slice_s)) {
      int s = slice_s[i];
      if (N_occ[s] > 0) lp += species_log_lik(s, occ, occ_start, N_occ, env_m, m_start, N_m,
                            mu[s]', exp(log_sigma[s]'), rho[s]);
    }
    return lp;
  }
}

data {
  int<lower=1> S;                    // number of species
  int<lower=1> P;                    // environmental dimensions (P = 2)

  // Concatenated occurrence data
  int<lower=0> N_occ_total;
  matrix[N_occ_total, P] occ;
  array[S] int<lower=1> occ_start;
  array[S] int<lower=0> N_occ;

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
  real<lower=0> rate_log_sigma_scale; // half-normal scale for BM rate on log_sigma
  real z_rho_anc_prior;               // prior mean for ancestral z_rho
  real<lower=0> z_rho_anc_sigma;      // prior sd for ancestral z_rho
  real<lower=0> rate_rho_scale;       // half-normal scale for BM rate on z_rho

  int<lower=1> grainsize;             // reduce_sum grainsize
}

transformed data {
  cholesky_factor_cov[S] L_C = cholesky_decompose(C);
  array[S] int species_idx;
  for (s in 1:S) species_idx[s] = s;
}

parameters {
  // Ancestral niche parameters
  vector[P] mu_anc;
  vector[P] log_sigma_anc;
  real z_rho_anc;

  // BM evolutionary rates (per dimension)
  vector<lower=0>[P] rate_mu;
  vector<lower=0>[P] rate_log_sigma;
  real<lower=0> rate_rho;

  // Non-centered species deviations (iid standard normal)
  matrix[S, P] mu_raw;
  matrix[S, P] log_sigma_raw;
  vector[S] z_rho_raw;
}

transformed parameters {
  // Per-species centroids, log tolerances and Fisher-z correlation via
  // non-centered BM.
  matrix[S, P] mu;
  matrix[S, P] log_sigma;
  vector[S] z_rho;
  vector[S] rho;

  for (k in 1:P) {
    mu[, k] = mu_anc[k] + sqrt(rate_mu[k]) * (L_C * mu_raw[, k]);
    log_sigma[, k] = log_sigma_anc[k]
                     + sqrt(rate_log_sigma[k]) * (L_C * log_sigma_raw[, k]);
  }
  z_rho = z_rho_anc + sqrt(rate_rho) * (L_C * z_rho_raw);
  rho = tanh(z_rho);
}

model {
  // --- Ancestral niche priors ---
  mu_anc ~ normal(mu_anc_prior, mu_anc_sigma);
  log_sigma_anc ~ normal(0, 3);
  z_rho_anc ~ normal(z_rho_anc_prior, z_rho_anc_sigma);

  // --- BM rate priors (half-normal) ---
  rate_mu ~ normal(0, rate_mu_scale);
  rate_log_sigma ~ normal(0, rate_log_sigma_scale);
  rate_rho ~ normal(0, rate_rho_scale);

  // --- Non-centered deviations: iid N(0,1) ---
  for (k in 1:P) {
    mu_raw[, k] ~ std_normal();
    log_sigma_raw[, k] ~ std_normal();
  }
  z_rho_raw ~ std_normal();

  // --- Per-species likelihood (parallelised with reduce_sum) ---
  target += reduce_sum(partial_sum_bm_evolving, species_idx, grainsize,
                       occ, occ_start, N_occ, env_m, m_start, N_m,
                       mu, log_sigma, rho);
}

generated quantities {
  // Per-species covariance and correlation matrices and per-species
  // log-likelihood for model comparison (LOO/WAIC).
  array[S] matrix[P, P] Sigma;
  array[S] matrix[P, P] R_corr;
  vector[S] log_lik;

  for (s in 1:S) {
    vector[P] mu_s = mu[s]';
    vector[P] sigma_s = exp(log_sigma[s]');
    real rho_s = rho[s];
    matrix[P, P] L_cov = cov_chol_from_rho(sigma_s, rho_s);

    if (N_occ[s] > 0) {
      log_lik[s] = species_log_lik(s, occ, occ_start, N_occ,
                                 env_m, m_start, N_m,
                                 mu_s, sigma_s, rho_s);
    } else {
      log_lik[s] = 0;
    }

    Sigma[s] = multiply_lower_tri_self_transpose(L_cov);
    R_corr[s, 1, 1] = 1.0;
    R_corr[s, 2, 2] = 1.0;
    R_corr[s, 1, 2] = rho_s;
    R_corr[s, 2, 1] = rho_s;
  }
}
