# arfun

**Ancestral Reconstruction of the Fundamental Niche**

`arfun` estimates species' environmental niches as bivariate normal
ellipses in environmental space — from presence-only data, normalized
against a user-supplied environmental background representing the
accessible area `M` — and fits how those ellipses evolve along a
phylogeny under Brownian motion (BM), Ornstein–Uhlenbeck (OU), or an
independent-species null model. The bivariate ellipsoid is a model-based
estimate of the fundamental niche conditional on the data and on `M`;
it is not a direct measurement of physiological tolerance.

## What it does

- **Presence-only niche estimation.** The Jiménez & Soberón likelihood
  fits an environmental ellipse `(mu, sigma, rho)` to occurrence points
  normalized against the accessible area `M`.
- **Niche evolution.** Joint models let `mu`, `sigma` and `rho` evolve
  along a phylogeny under BM or OU, so each species' ellipse borrows
  strength from its relatives.
- **Boundary diagnostics.** A multi-start MAP scan reports which
  species sit on estimation boundaries (`sigma -> 0`, `|rho| -> 1`,
  `sigma -> Inf`) before any MCMC is run — and supports *masking* poorly
  identified species so the phylogeny predicts their niche as a
  validation target.
- **Warm starts via xnicher.** `niche_warmstart()` fits each species
  independently with `xnicher` and injects the result into the joint
  optimization.
- **Virtual species.** Simulate niche evolution (BM / OU / null) inside a
  real accessible area, draw occurrences, and check recovery — the
  reverse-engineering loop used to validate the whole pipeline.

## Installation

```r
# hard dependencies (all on CRAN)
install.packages(c("ucminfcpp", "xnicher"))

# development version
remotes::install_github("alrobles/arfun")
```

Stan models compile with
[cmdstanr](https://mc-stan.org/cmdstanr/) for full Bayesian inference
(requires a working CmdStan installation; `install.packages("cmdstanr",
repos = c("https://stan-dev.r-universe.dev", getOption("repos")))`);
the MAP path needs only `ucminfcpp`.

## Quick start

```r
library(arfun)

# per-species occurrence and shared-background environmental matrices
data <- prepare_phylo_niche_data(env_occ_list, env_m_list, tree)

# multi-start MAP on the exact joint log-posterior (C++ backend)
fit <- niche_multistart(data, n_starts = 32)
fit$boundary   # per-species diagnostics

# full Bayesian fit of the BM evolving-sigma model
evo <- fit_evolution(model = "bm_evolving", niche_data = data)
```

See `vignette("virtual-species-bioregions")` for an end-to-end
simulate-and-recover benchmark inside a real One Earth bioregion.

## The models

| model | niche evolution |
|---|---|
| `bm_evolving` | BM on centroid, tolerances and correlation |
| `ou_evolving` | OU toward an ancestral optimum on all parameters |
| `bm_constant`, `ou_constant` | evolution on the centroid with a shared shape |
| `null` | independent species, no phylogenetic signal |
| `presence_only` | single-species J&S ellipse |

## Formal verification (experimental)

The optional `formal/` subproject checks selected mathematical statements in
Lean 4.33.1 with a pinned mathlib revision. Its initial 26 lemmas/theorems cover
bivariate covariance geometry, likelihood normalization, a conditional
quadrature-error bound, and matrix/algebraic identities used in BM and OU.
It does not change the R/Stan fitting pipeline or add an R installation dependency.

With Elan/Lake and Python 3.11 or newer, run from the repository root:

```bash
python3 formal/check.py --bootstrap
```

For subsequent checks, omit `--bootstrap`. The verifier builds the library,
audits theorem axioms, replays the project modules with `leanchecker`, and
checks that deliberately invalid proof controls are rejected.

See the [beginner's manual in Spanish](formal/MANUAL.es.md) and the
[theorem-to-equation registry](formal/theorems.json). These are proofs of the
stated real-number specifications, not a certification of the entire Stan
executable, MCMC convergence, or ecological assumptions.

## Citation

Jiménez, L., Soberón, J., Christen, J.A. & Soto, D. (2019). On the
problem of modeling a fundamental niche from occurrence data.
*Ecological Modelling* 397: 74–83.
