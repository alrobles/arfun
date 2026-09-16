// niche_ou_constant_sigma.stan
// arfunEvolution: Ornstein-Uhlenbeck (OU) niche evolution with constant Sigma.
//
// The niche centroid mu[s] evolves along the phylogeny under an
// Ornstein-Uhlenbeck (OU) process (Hansen 1997). The niche shape --
// tolerances (log sigma) and correlation (rho, Fisher-z) -- is shared by all
// species: a single, constant covariance matrix Sigma. The OU process
// introduces a selection strength alpha and an optimum (the ancestral
// centroid mu_anc), and converges to BM as alpha -> 0.
//
// Per-species likelihoods use the unweighted background correction of
// Jimenez & Soberon (2019). Designed for the bivariate niche ellipsoid
// (P = 2) of the 'arfun' package. reduce_sum parallelises the per-species
// likelihood.

functions {
  // OU variance-covariance matrix for one trait dimension (Hansen 1997):
  //   V_ij = sigma_sq * exp(-alpha * tau_ij)
  //          * (1 - exp(-2 * alpha * tanc_ij)) / (2 * alpha)
  // Uses expm1 for numerical stability when alpha is close to 0.
  matrix ou_vcv(matrix tanc, real alpha, real sigma_sq) {
    int S = rows(tanc);
    vector[S] t_root = diagonal(tanc);
    matrix[S, S] V;
    real two_alpha = 2 * alpha;
    for (i in 1:S) {
      for (j in 1:S) {
        real dist = t_root[i] + t_root[j] - 2 * tanc[i, j];
        real varanc = -0.5 * expm1(-two_alpha * tanc[i, j]) / alpha;
        V[i, j] = sigma_sq * exp(-alpha * dist) * varanc;
      }
    }
    return V;
  }

  // Lower-triangular Cholesky factor of the niche covariance for a bivariate
  // (P = 2) ellipsoid with standard deviations sigma and correlation rho.
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
  real partial_sum_ou_constant(array[] int slice_s, int start, int end,
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

  // Shared branch-length matrix from ape::vcv(tree). Diagonal is root-to-tip
  // time for each tip; off-diagonal is time from root to the MRCA of the pair.
  matrix[S, S] T_anc;

  // Prior hyperparameters
  vector[P] mu_anc_prior;             // prior mean for ancestral centroid
  vector<lower=0>[P] mu_anc_sigma;    // prior sd for ancestral centroid
  vector<lower=0>[P] alpha_mu_scale;  // half-normal scale for OU alpha on centroids
  vector<lower=0>[P] sigma_sq_mu_scale; // half-normal scale for OU sigma^2 on centroids
  real z_rho_prior;                   // prior mean for shared z_rho
  real<lower=0> z_rho_sigma;          // prior sd for shared z_rho

  int<lower=1> grainsize;             // reduce_sum grainsize
}

transformed data {
  array[S] int species_idx;
  for (s in 1:S) species_idx[s] = s;
}

parameters {
  // Ancestral niche centroid (also the OU optimum)
  vector[P] mu_anc;

  // Shared (constant) niche shape
  vector[P] log_sigma;
  real z_rho;

  // OU evolutionary parameters for centroids (per dimension)
  // alpha is bounded away from 0 to avoid the 0/0 singularity in ou_vcv.
  vector<lower=1e-6>[P] alpha_mu;
  vector<lower=0>[P] sigma_sq_mu;

  // Non-centered centroid deviations (iid standard normal)
  matrix[S, P] mu_raw;
}

transformed parameters {
  // Per-species centroids via non-centered OU.
  matrix[S, P] mu;
  // Constant niche shape (shared across species)
  vector[P] sigma = exp(log_sigma);
  real rho = tanh(z_rho);

  for (k in 1:P) {
    matrix[S, S] V_mu_k = ou_vcv(T_anc, alpha_mu[k], sigma_sq_mu[k]);
    matrix[S, S] L_V_mu_k = cholesky_decompose(V_mu_k);
    mu[, k] = mu_anc[k] + L_V_mu_k * mu_raw[, k];
  }
}

model {
  // --- Priors on ancestral and shared shape parameters ---
  mu_anc ~ normal(mu_anc_prior, mu_anc_sigma);
  log_sigma ~ normal(0, 3);
  z_rho ~ normal(z_rho_prior, z_rho_sigma);

  // --- Priors on OU parameters (half-normal) ---
  alpha_mu ~ normal(0, alpha_mu_scale);
  sigma_sq_mu ~ normal(0, sigma_sq_mu_scale);

  // --- Non-centered OU deviations: iid N(0,1) ---
  for (k in 1:P)
    mu_raw[, k] ~ std_normal();

  // --- Per-species likelihood (parallelised with reduce_sum) ---
  target += reduce_sum(partial_sum_ou_constant, species_idx, grainsize,
                       occ, occ_start, N_occ, env_m, m_start, N_m,
                       mu, sigma, rho);
}

generated quantities {
  // Per-species covariance and correlation matrices (constant across species
  // here) and per-species log-likelihood for model comparison (LOO/WAIC).
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
