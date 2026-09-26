# arfun 0.3.1

* Fix: `niche_warmstart()$per_species` now returns `mu2` (the second
  centroid coordinate) alongside `mu1`, `sigma1`, `sigma2`, `rho`.

# arfun 0.3.0

* `niche_logpost()`: exact log-posterior on the mathematical
  (unconstrained) scale, matching the four Stan parameterizations.
* Compiled C++ objective via `niche_logpost_xptr()` +
  `ucminfcpp::ucminf_xptr` for fast MAP estimation.
* `niche_multistart()`: Sobol/LHS multistart MAP with boundary
  diagnostics (`niche_boundary_report()`).
* `niche_warmstart()`: per-species independent fits via `xnicher`
  assembled into a joint-model start point.
* `prepare_phylo_niche_data()`: `mask_species` keeps weakly identified
  species in the phylogenetic prior while excluding them from the
  likelihood (mask-and-predict validation).
* `simulate_niche_evolution()`: BM/OU/null niche simulation on a tree.
* `shared_background_sample()`: a single shared background sample for
  all species, as required by the J&S likelihood under a common
  accessible area M.
* Vignette: virtual species evolving inside a real One Earth bioregion.
* All Stan models accept `N_occ = 0` (masked species); they still emit
  predicted niche parameters in `generated quantities`.
