# arfun

**Ancestral Reconstruction of the Fundamental Niche**

`arfun` estimates species' fundamental niches as bivariate normal
ellipses in environmental space — from presence-only data, correcting for
the set of environments actually available to the species — and fits how
those ellipses evolve along a phylogeny under Brownian motion (BM),
Ornstein–Uhlenbeck (OU), or an independent-species null model.

## What it does

- **Presence-only niche estimation.** The Jiménez & Soberón likelihood
  fits an environmental ellipse `(mu, sigma, rho)` to occurrence points
  normalized against the accessible area `M`.
- **Niche evolution.** Joint models let `mu`, `sigma` and `rho` evolve
  along a phylogeny under BM or OU, so each species' ellipse borrows
  strength from its relatives.
- **Identifiability diagnostics.** A multi-start MAP scan reports which
  species sit on estimation boundaries (`sigma -> 0`, `|rho| -> 1`,
  `sigma -> Inf`) before any MCMC is run — and supports *masking* poorly
  identified species so the phylogeny predicts their niche as a
  validation target.
- **Virtual species.** Simulate niche evolution (BM / OU / null) inside a
  real accessible area, draw occurrences, and check recovery — the
  reverse-engineering loop used to validate the whole pipeline.

## Installation

```r
# development version
remotes::install_github("alrobles/arfun")

# fast MAP optimizer backend
remotes::install_github("alrobles/ucminfcpp")
```

Stan models compile with
[cmdstanr](https://mc-stan.org/cmdstanr/) for full Bayesian inference;
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

## Citation

Jiménez, L., Soberón, J., Christen, J.A. & Soto, D. (2019). On the
problem of modeling a fundamental niche from occurrence data.
*Ecological Modelling* 397: 74–83.
