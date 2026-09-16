// niche_logpost_xptr.cpp
//
// Pure-C++ ObjFun for ucminfcpp::ucminf_xptr: the exact Stan log-posterior
// of the 2D phylo niche model on the unconstrained ("math") scale.
// Ported verbatim from logpost_math() in 17_multistart_map.R (which was
// verified term-by-term against the param-less lp_debug_centered.stan
// probe: priors/jacobian ~1e-6, likelihood ~0.03 float noise).
//
//   - species ll: J&S, lower-triangular 2x2 Cholesky solve
//   - priors: mu_anc ~ N(mu_anc_prior, mu_anc_sigma),
//             log_sigma_anc ~ N(0,3), z_rho_anc ~ N(prior, sigma),
//             rates: half-normal on exp(lrate) + Jacobian sum(lrate)
//                    OR (ctrlr) lognormal: N(log(scale), 1) on lrate
//   - trait prior: unb/bnd raws ~ N(0,1) then affine through L_C;
//                  ctr/ctrlr MVN(anc, rate*C + jitter*I)
//   - rho = rho_cap * tanh(z_rho); sigma floor rejects the degenerate
//     sigma->0 spike of the J&S likelihood.
//
// Theta layout (same as unpack() in the R script; Stan decl order):
//   mu_anc[P] | log_sigma_anc[P] | lrate_mu[P] | lrate_ls[P] |
//   matA[S*P col-major] | matB[S*P] | z_rho_anc | lrate_rho | vec[S]
//   unb/bnd: matA=mu_raw, matB=log_sigma_raw, vec=z_rho_raw
//   ctr/ctrlr: matA=mu,    matB=log_sigma,     vec=z_rho
//
// Compile from R:
//   Rcpp::sourceCpp("niche_logpost_xptr.cpp")
//   xp <- make_niche_logpost_xptr(occ, occ_start, N_occ, env_m, m_start,
//         N_m, C, mu_anc_prior, mu_anc_sigma, z_rho_anc_prior,
//         z_rho_anc_sigma, rate_mu_scale, rate_ls_scale, rate_rho_scale,
//         model = "ctr", rho_cap = 0.98, sigma_floor = 0.02)
//   res <- ucminfcpp::ucminf_xptr(theta0, xp, control = list(...))

#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "ucminf_core.hpp"

namespace {

constexpr double kLog2Pi = 1.8378770664093453;  // log(2*pi)
constexpr double kRho2Eps = 1e-12;              // pmax(1-rho^2, eps)

inline double dnorm_lp(double x, double mu, double sd) {
    const double z = (x - mu) / sd;
    return -0.5 * z * z - std::log(sd) - 0.5 * kLog2Pi;
}

// Lower-triangular Cholesky of a column-major symmetric n x n matrix.
// Returns false if not positive definite.
bool chol_lower(const std::vector<double>& A, int n,
                std::vector<double>& L) {
    L.assign(static_cast<size_t>(n) * n, 0.0);
    for (int j = 0; j < n; ++j) {
        double d = A[j * n + j];
        for (int k = 0; k < j; ++k) d -= L[k * n + j] * L[k * n + j];
        if (!(d > 0.0) || !std::isfinite(d)) return false;
        L[j * n + j] = std::sqrt(d);
        for (int i = j + 1; i < n; ++i) {
            double s = A[j * n + i];
            for (int k = 0; k < j; ++k) s -= L[k * n + i] * L[k * n + j];
            L[j * n + i] = s / L[j * n + j];
        }
    }
    return true;
}

// Solve L x = b for lower-triangular column-major L (n x n).
void solve_lower(const std::vector<double>& L, int n,
                 const double* b, double* x) {
    for (int i = 0; i < n; ++i) {
        double s = b[i];
        for (int k = 0; k < i; ++k) s -= L[k * n + i] * x[k];
        x[i] = s / L[i * n + i];
    }
}

// y = L x for lower-triangular column-major L (n x n).
void matvec_lower(const std::vector<double>& L, int n,
                  const double* x, double* y) {
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int k = 0; k <= i; ++k) s += L[k * n + i] * x[k];
        y[i] = s;
    }
}

// MVN log-density: x ~ N(mean_i = mu_scalar, Sigma) via Cholesky.
double dmvn_lp(const double* x, double mu_scalar,
               const std::vector<double>& Sigma, int n,
               std::vector<double>& L, std::vector<double>& z) {
    if (!chol_lower(Sigma, n, L)) return -std::numeric_limits<double>::infinity();
    std::vector<double> d(n);
    for (int i = 0; i < n; ++i) d[i] = x[i] - mu_scalar;
    solve_lower(L, n, d.data(), z.data());
    double q = 0.0, ld = 0.0;
    for (int i = 0; i < n; ++i) {
        q += z[i] * z[i];
        ld += std::log(L[i * n + i]);
    }
    return -0.5 * q - ld - 0.5 * n * kLog2Pi;
}

// Streaming log-sum-exp of -0.5 * sum q^2 over background rows.
struct LseAcc {
    double m = -std::numeric_limits<double>::infinity();
    double s = 0.0;
    void add(double v) {
        if (v > m) { s = s * std::exp(m - v) + 1.0; m = v; }
        else       { s += std::exp(v - m); }
    }
    double get() const { return m + std::log(s); }
};

enum class ModelKind : int { UNB = 0, BND = 1, CTR = 2, CTRLR = 3 };

struct NicheState {
    // --- data --------------------------------------------------------------
    int S = 0;
    int P = 0;
    std::vector<double> occ;        // col-major (tot_occ, P)
    int                 occ_nrow = 0;
    std::vector<int>    occ_start;  // 1-based from R -> converted to 0-based
    std::vector<int>    N_occ;
    std::vector<double> env_m;      // col-major (tot_bg, P)
    int                 bg_nrow = 0;
    std::vector<int>    m_start;
    std::vector<int>    N_m;
    std::vector<double> C;          // col-major S x S
    std::vector<double> L_C;        // lower chol of C (precomputed)
    double jitter = 1e-9;

    // --- priors ------------------------------------------------------------
    std::vector<double> mu_anc_prior;
    std::vector<double> mu_anc_sigma;   // per-dim (length P)
    double z_rho_anc_prior = 0.0;
    double z_rho_anc_sigma = 1.0;
    double rate_mu_scale  = 1.0;
    double rate_ls_scale  = 1.0;
    double rate_rho_scale = 1.0;

    // --- config ------------------------------------------------------------
    ModelKind model = ModelKind::CTR;
    double rho_cap     = 0.98;
    double sigma_floor = 0.02;
    double gradstep_rel = 1e-6;
    double gradstep_abs = 1e-8;
    bool   use_central  = true;

    // --- scratch ------------------------------------------------------------
    std::vector<double> mu, log_sigma, z_rho;      // derived traits
    std::vector<double> sig_work, chol_work, z_work; // MVN scratch
};

// J&S species log-likelihood (2D). Returns -Inf on failure.
double species_ll(const NicheState& st, int s,
                  const double* mu_s,       // P
                  const double* log_sigma_s,// P
                  double rho) {
    const double sg1 = std::exp(log_sigma_s[0]);
    const double sg2 = std::exp(log_sigma_s[1]);
    // lower-triangular 2x2: L00 = sg1, L10 = rho*sg2,
    //                       L11 = sg2*sqrt(1-rho^2)
    const double L11 = sg2 * std::sqrt(std::max(1.0 - rho * rho, kRho2Eps));
    const double L10 = rho * sg2;
    if (!(sg1 > 0.0) || !(L11 > 0.0)) {
        return -std::numeric_limits<double>::infinity();
    }

    const int os = st.occ_start[s];
    const int no = st.N_occ[s];
    const int bs = st.m_start[s];
    const int nb = st.N_m[s];
    // Masked species carry no occurrence/background rows: zero likelihood.
    if (no == 0) return 0.0;

    // occurrence term: -0.5 * sum q^2 ; q = L^{-1}(x - mu)
    double q_occ = 0.0;
    for (int i = 0; i < no; ++i) {
        const double d0 = st.occ[os + i]                      - mu_s[0];
        const double d1 = st.occ[os + i + (size_t)st.occ_nrow] - mu_s[1];
        const double y0 = d0 / sg1;
        const double y1 = (d1 - L10 * y0) / L11;
        q_occ += y0 * y0 + y1 * y1;
    }

    // background normalization: -N_occ * lse(-0.5 * q_bg)
    LseAcc acc;
    for (int i = 0; i < nb; ++i) {
        const double d0 = st.env_m[bs + i]                     - mu_s[0];
        const double d1 = st.env_m[bs + i + (size_t)st.bg_nrow] - mu_s[1];
        const double y0 = d0 / sg1;
        const double y1 = (d1 - L10 * y0) / L11;
        acc.add(-0.5 * (y0 * y0 + y1 * y1));
    }
    return -0.5 * q_occ - static_cast<double>(no) * acc.get();
}

// Exact log-posterior on math scale. -Inf if outside support.
double niche_lp_eval(const NicheState& st, const double* th) {
    const int S = st.S, P = st.P;
    const bool noncentered =
        (st.model == ModelKind::UNB || st.model == ModelKind::BND);
    const bool ctrlr = (st.model == ModelKind::CTRLR);

    int i = 0;
    const double* mu_anc        = th + i; i += P;
    const double* log_sigma_anc = th + i; i += P;
    const double* lrate_mu      = th + i; i += P;
    const double* lrate_ls      = th + i; i += P;
    const double* matA          = th + i; i += S * P; // mu_raw | mu
    const double* matB          = th + i; i += S * P; // log_sigma_raw | log_sigma
    const double  z_rho_anc     = th[i++];
    const double  lrate_rho     = th[i++];
    const double* vec           = th + i; i += S;     // z_rho_raw | z_rho

    // ---- ancestor priors ----------------------------------------------------
    double lp = 0.0;
    for (int k = 0; k < P; ++k) {
        lp += dnorm_lp(mu_anc[k], st.mu_anc_prior[k], st.mu_anc_sigma[k]);
        lp += dnorm_lp(log_sigma_anc[k], 0.0, 3.0);
    }
    lp += dnorm_lp(z_rho_anc, st.z_rho_anc_prior, st.z_rho_anc_sigma);

    // ---- rate priors ---------------------------------------------------------
    if (ctrlr) {
        // declared log_rate params: lognormal prior, no Jacobian
        for (int k = 0; k < P; ++k) {
            lp += dnorm_lp(lrate_mu[k], std::log(st.rate_mu_scale), 1.0);
            lp += dnorm_lp(lrate_ls[k], std::log(st.rate_ls_scale), 1.0);
        }
        lp += dnorm_lp(lrate_rho, std::log(st.rate_rho_scale), 1.0);
    } else {
        // half-normal on rate = exp(lrate) + Jacobian +sum(lrate)
        for (int k = 0; k < P; ++k) {
            lp += dnorm_lp(std::exp(lrate_mu[k]), 0.0, st.rate_mu_scale)
                + lrate_mu[k];
            lp += dnorm_lp(std::exp(lrate_ls[k]), 0.0, st.rate_ls_scale)
                + lrate_ls[k];
        }
        lp += dnorm_lp(std::exp(lrate_rho), 0.0, st.rate_rho_scale)
            + lrate_rho;
    }

    // ---- traits --------------------------------------------------------------
    NicheState& w = const_cast<NicheState&>(st); // scratch only
    if (noncentered) {
        for (size_t j = 0; j < (size_t)S * P; ++j) {
            lp += dnorm_lp(matA[j], 0.0, 1.0);
            lp += dnorm_lp(matB[j], 0.0, 1.0);
        }
        for (int s = 0; s < S; ++s) lp += dnorm_lp(vec[s], 0.0, 1.0);
        w.mu.assign(S * P, 0.0);
        w.log_sigma.assign(S * P, 0.0);
        w.z_rho.assign(S, 0.0);
        std::vector<double> tmp(S);
        for (int k = 0; k < P; ++k) {
            matvec_lower(st.L_C, S, matA + (size_t)k * S, tmp.data());
            const double rmu = std::exp(0.5 * lrate_mu[k]);
            for (int s = 0; s < S; ++s)
                w.mu[k * S + s] = mu_anc[k] + rmu * tmp[s];
            matvec_lower(st.L_C, S, matB + (size_t)k * S, tmp.data());
            const double rls = std::exp(0.5 * lrate_ls[k]);
            for (int s = 0; s < S; ++s)
                w.log_sigma[k * S + s] = log_sigma_anc[k] + rls * tmp[s];
        }
        matvec_lower(st.L_C, S, vec, tmp.data());
        const double rr = std::exp(0.5 * lrate_rho);
        for (int s = 0; s < S; ++s)
            w.z_rho[s] = z_rho_anc + rr * tmp[s];
    } else {
        w.mu.assign(matA, matA + (size_t)S * P);
        w.log_sigma.assign(matB, matB + (size_t)S * P);
        w.z_rho.assign(vec, vec + S);
        // centered MVN priors: x_k ~ N(anc_k, rate_k * C + jitter I)
        w.sig_work.assign((size_t)S * S, 0.0);
        w.z_work.assign(S, 0.0);
        for (int k = 0; k < P; ++k) {
            const double rmu = std::exp(lrate_mu[k]);
            const double rls = std::exp(lrate_ls[k]);
            for (int j = 0; j < S * S; ++j)
                w.sig_work[j] = rmu * st.C[j] + (j % (S + 1) == 0 ? st.jitter : 0.0);
            lp += dmvn_lp(w.mu.data() + (size_t)k * S, mu_anc[k],
                          w.sig_work, S, w.chol_work, w.z_work);
            for (int j = 0; j < S * S; ++j)
                w.sig_work[j] = rls * st.C[j] + (j % (S + 1) == 0 ? st.jitter : 0.0);
            lp += dmvn_lp(w.log_sigma.data() + (size_t)k * S,
                          log_sigma_anc[k], w.sig_work, S,
                          w.chol_work, w.z_work);
        }
        const double rr = std::exp(lrate_rho);
        for (int j = 0; j < S * S; ++j)
            w.sig_work[j] = rr * st.C[j] + (j % (S + 1) == 0 ? st.jitter : 0.0);
        lp += dmvn_lp(w.z_rho.data(), z_rho_anc, w.sig_work, S,
                      w.chol_work, w.z_work);
    }
    if (!std::isfinite(lp)) return -std::numeric_limits<double>::infinity();

    // ---- regularization + likelihood ----------------------------------------
    const double cap = (st.model == ModelKind::UNB) ? 1.0 : st.rho_cap;
    for (size_t j = 0; j < w.log_sigma.size(); ++j) {
        if (std::exp(w.log_sigma[j]) < st.sigma_floor)
            return -std::numeric_limits<double>::infinity();
    }
    // mu/log_sigma are stored col-major (S x P): element (s,k) at k*S + s.
    for (int s = 0; s < S; ++s) {
        const double rho = cap * std::tanh(w.z_rho[s]);
        if (!std::isfinite(rho))
            return -std::numeric_limits<double>::infinity();
        double mus[2] = { w.mu[s], w.mu[S + s] };
        double lss[2] = { w.log_sigma[s], w.log_sigma[S + s] };
        const double ll = species_ll(st, s, mus, lss, rho);
        if (!std::isfinite(ll))
            return -std::numeric_limits<double>::infinity();
        lp += ll;
    }
    return lp;
}

} // namespace

// ---------------------------------------------------------------------------
// make_niche_logpost_xptr
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP make_niche_logpost_xptr(
    Rcpp::NumericMatrix occ,
    Rcpp::IntegerVector occ_start,
    Rcpp::IntegerVector N_occ,
    Rcpp::NumericMatrix env_m,
    Rcpp::IntegerVector m_start,
    Rcpp::IntegerVector N_m,
    Rcpp::NumericMatrix C,
    Rcpp::NumericVector mu_anc_prior,
    Rcpp::NumericVector mu_anc_sigma,
    double z_rho_anc_prior,
    double z_rho_anc_sigma,
    double rate_mu_scale,
    double rate_ls_scale,
    double rate_rho_scale,
    std::string model = "ctr",
    double rho_cap = 0.98,
    double sigma_floor = 0.02,
    double jitter = 1e-9,
    std::string grad = "central",
    Rcpp::NumericVector gradstep = Rcpp::NumericVector::create(1e-6, 1e-8)
) {
    if (gradstep.size() != 2) Rcpp::stop("`gradstep` must have length 2.");
    if (grad != "forward" && grad != "central")
        Rcpp::stop("`grad` must be 'forward' or 'central'.");
    const int P = occ.ncol();
    if (P != 2 || env_m.ncol() != 2)
        Rcpp::stop("This objective is the 2D (P=2) niche model.");
    const int S = N_occ.size();
    if (C.nrow() != S || C.ncol() != S)
        Rcpp::stop("`C` must be S x S.");
    if (mu_anc_prior.size() != P)
        Rcpp::stop("`mu_anc_prior` must have length P.");
    if (mu_anc_sigma.size() != P)
        Rcpp::stop("`mu_anc_sigma` must have length P.");

    auto st = std::make_shared<NicheState>();
    st->S = S; st->P = P;
    st->occ_nrow = occ.nrow();
    st->bg_nrow  = env_m.nrow();
    st->occ.assign(occ.begin(), occ.end());
    st->occ_start.assign(occ_start.begin(), occ_start.end());
    st->N_occ.assign(N_occ.begin(), N_occ.end());
    st->env_m.assign(env_m.begin(), env_m.end());
    st->m_start.assign(m_start.begin(), m_start.end());
    st->N_m.assign(N_m.begin(), N_m.end());
    st->C.assign(C.begin(), C.end());
    st->mu_anc_prior.assign(mu_anc_prior.begin(), mu_anc_prior.end());
    st->mu_anc_sigma.assign(mu_anc_sigma.begin(), mu_anc_sigma.end());
    st->z_rho_anc_prior = z_rho_anc_prior;
    st->z_rho_anc_sigma = z_rho_anc_sigma;
    st->rate_mu_scale   = rate_mu_scale;
    st->rate_ls_scale   = rate_ls_scale;
    st->rate_rho_scale  = rate_rho_scale;
    st->jitter          = jitter;
    st->rho_cap         = rho_cap;
    st->sigma_floor     = sigma_floor;
    st->use_central     = (grad == "central");
    st->gradstep_rel    = gradstep[0];
    st->gradstep_abs    = gradstep[1];

    if (model == "unb")       st->model = ModelKind::UNB;
    else if (model == "bnd")  st->model = ModelKind::BND;
    else if (model == "ctr")  st->model = ModelKind::CTR;
    else if (model == "ctrlr") st->model = ModelKind::CTRLR;
    else Rcpp::stop("Unknown model '%s' (unb|bnd|ctr|ctrlr).", model);

    // Precompute L_C once (needed by the non-centered transform).
    if (!chol_lower(st->C, S, st->L_C))
        Rcpp::stop("`C` is not positive definite.");

    // 1-based R indices -> 0-based row offsets.
    for (int s = 0; s < S; ++s) {
        st->occ_start[s] -= 1;
        st->m_start[s]   -= 1;
    }

    auto fn = std::make_unique<ucminf::ObjFun>(
        [st](const std::vector<double>& x,
             std::vector<double>&       g,
             double&                    f)
        {
            if (g.size() != x.size()) g.resize(x.size(), 0.0);
            const double lp0 = niche_lp_eval(*st, x.data());
            if (!std::isfinite(lp0)) {
                f = std::numeric_limits<double>::infinity();
                std::fill(g.begin(), g.end(), 0.0);
                return;
            }
            f = -lp0;

            // Finite-difference gradient over all coordinates.
            std::vector<double> xp = x;
            for (size_t j = 0; j < x.size(); ++j) {
                const double dx =
                    std::abs(x[j]) * st->gradstep_rel + st->gradstep_abs;
                xp[j] = x[j] + dx;
                const double fp = niche_lp_eval(*st, xp.data());
                if (st->use_central) {
                    xp[j] = x[j] - dx;
                    const double fm = niche_lp_eval(*st, xp.data());
                    g[j] = (std::isfinite(fp) && std::isfinite(fm))
                        ? -(fp - fm) / (2.0 * dx)
                        : (std::isfinite(fp) ? -(fp - lp0) / dx : 0.0);
                } else {
                    g[j] = std::isfinite(fp) ? -(fp - lp0) / dx : 0.0;
                }
                xp[j] = x[j];
            }
        });

    return Rcpp::XPtr<ucminf::ObjFun>(fn.release(), true);
}

// ---------------------------------------------------------------------------
// xptr_niche_logpost_eval: evaluate the compiled lp at one point
// (parity check against the R implementation; returns lp, not -lp)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double xptr_niche_logpost_eval(SEXP xptr, Rcpp::NumericVector theta) {
    Rcpp::XPtr<ucminf::ObjFun> fn(xptr);
    std::vector<double> x = Rcpp::as<std::vector<double>>(theta);
    std::vector<double> g(x.size(), 0.0);
    double f = 0.0;
    (*fn)(x, g, f);
    return -f;
}
