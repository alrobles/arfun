// niche_null_model.stan
// arfunEvolution: null model (independent species, no phylogenetic structure).
//
// Each species has an independent niche centroid mu[s], log tolerances
// log_sigma[s], and correlation rho[s] (via z_rho = atanh(rho)). There is no
// phylogenetic structure: species are iid draws from a common prior. This is
// the baseline model against which the phylogenetic BM/OU models are compared
// (e.g. via LOO) to test whether phylogeny improves prediction.
//
// Per-species likelihoods use the unweighted background correction of
// Jimenez & Soberon (2019). Designed for the bivariate niche ellipsoid
// (P = 2) of the 'arfun' package. reduce_sum parallelises the per-species
// likelihood.

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
  // Jimenez & Soberon (2019) model.
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

    matrix[N_occ[s], P] occ_s = occ[s_start:s_end, ];
    matrix[P, N_occ[s]] diff_occ = occ_s' - rep_matrix(mu_s, N_occ[s]);
    matrix[P, N_occ[s]] y_occ = mdivide_left_tri_low(L_cov, diff_occ);
    real sum_q_occ = sum(columns_dot_self(y_occ));

    matrix[N_m[s], P] env_m_s = env_m[b_start:b_end, ];
    matrix[P, N_m[s]] diff_m = env_m_s' - rep_matrix(mu_s, N_m[s]);
    matrix[P, N_m[s]] y_m = mdivide_left_tri_low(L_cov, diff_m);
    vector[N_m[s]] a = (-0.5 * columns_dot_self(y_m))';

    return -0.5 * sum_q_occ - N_occ[s] * log_sum_exp(a);
  }

  // reduce_sum partial function over a slice of species.
  real partial_sum_null(array[] int slice_s, int start, int end,
                        matrix occ, array[] int occ_start, array[] int N_occ,
                        matrix env_m, array[] int m_start, array[] int N_m,
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

  // Prior hyperparameters (shared across species)
  vector[P] mu_prior;                 // prior mean for per-species centroid
  vector<lower=0>[P] mu_sigma;        // prior sd for per-species centroid
  real z_rho_prior;                   // prior mean for per-species z_rho
  real<lower=0> z_rho_sigma;          // prior sd for per-species z_rho

  int<lower=1> grainsize;             // reduce_sum grainsize
}

transformed data {
  array[S] int species_idx;
  for (s in 1:S) species_idx[s] = s;
}

parameters {
  // Per-species niche parameters (independent across species)
  matrix[S, P] mu;
  matrix[S, P] log_sigma;
  vector[S] z_rho;
}

transformed parameters {
  // Back-transformed per-species correlation (kept in (-1, 1))
  vector[S] rho = tanh(z_rho);
}

model {
  // --- Independent per-species priors (no phylogenetic structure) ---
  for (s in 1:S) {
    mu[s] ~ normal(mu_prior, mu_sigma);
    log_sigma[s] ~ normal(0, 3);
  }
  z_rho ~ normal(z_rho_prior, z_rho_sigma);

  // --- Per-species likelihood (parallelised with reduce_sum) ---
  target += reduce_sum(partial_sum_null, species_idx, grainsize,
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
