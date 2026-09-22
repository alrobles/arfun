#!/usr/bin/env Rscript
# Build the cached bm_evolving MCMC fit for the empirical-data vignette.
# Runs the exact same data-prep code the vignette teaches, on the bundled
# Vulpes extdata, then saves posterior draws + summary to an .rds that ships
# in inst/extdata/vulpes/vulpes_bm_evolving_cache.rds.
#
# Run inside the xsdm container on KU HPC:
#   R_LIBS_USER=~/R/arfun_lib:~/.local/lib/R apptainer exec $SIF Rscript vulpes_mcmc_cache.R

suppressPackageStartupMessages({
  library(ape); library(terra); library(arfun)
})

EXT <- Sys.getenv("EXT_DIR", "/home/a474r867/vulpes_extdata")
OUT <- Sys.getenv("OUT_RDS",
                  "/home/a474r867/vulpes_extdata/vulpes_bm_evolving_cache.rds")

## ---- same data prep as the vignette --------------------------------------
tree <- read.tree(file.path(EXT, "vulpes_tree.nwk"))
tree$tip.label <- gsub("_", " ", tree$tip.label)

occ <- read.csv(file.path(EXT, "vulpes_occurrences.csv"))
occ$species <- as.character(occ$species)

clim <- rast(file.path(EXT, "vulpes_climate.tif"))
m <- vect(file.path(EXT, "vulpes_m_polygons.gpkg"))

pts <- vect(occ, geom = c("lon", "lat"), crs = "EPSG:4326")
env_at_occ <- extract(clim, pts, ID = FALSE)
env_occ_list <- split(as.data.frame(env_at_occ), occ$species)
env_occ_list <- lapply(env_occ_list,
                       function(x) as.matrix(x[complete.cases(x), ]))

# z-score using clade-pooled occurrence moments (documented in vignette)
all_occ <- do.call(rbind, env_occ_list)
mu0 <- colMeans(all_occ); sd0 <- apply(all_occ, 2, sd)
env_occ_list <- lapply(env_occ_list, function(x) scale(x, mu0, sd0))

# background: sample inside each species' M polygon
bg_list <- lapply(seq_len(nrow(m)), function(i) {
  clim_m <- mask(crop(clim, m[i]), m[i])
  s <- spatSample(clim_m, size = 2000, method = "regular",
                  na.rm = TRUE, warn = FALSE)
  as.matrix(s[complete.cases(s), ])
})
names(bg_list) <- m$species
bg_list <- lapply(bg_list, function(x) scale(x, mu0, sd0))

mask_sp <- names(which(sapply(env_occ_list, nrow) < 10))
cat("masked species:", paste(mask_sp, collapse = ", "), "\n")

stan_data <- prepare_phylo_niche_data(env_occ_list, bg_list, tree,
                                      mask_species = mask_sp)

## ---- fit ------------------------------------------------------------------
# fit_evolution() applies stan_only() internally: bookkeeping fields such as
# `species`/`masked` are dropped before reaching cmdstanr.
fit <- fit_evolution(stan_data, model = "bm_evolving",
                     chains = 2, iter_warmup = 1000, iter_sampling = 1000,
                     seed = 7, mod = compile_model("bm_evolving"))

summ <- fit$summary()
draws <- posterior::as_draws_df(fit$draws())

## ---- post-hoc node posterior (BM conditional mean | tips, per draw) -------
S <- stan_data$S
depth <- ape::node.depth.edgelength(tree)
mrc <- ape::mrca(tree, full = TRUE)
V <- matrix(depth[mrc], nrow = nrow(mrc))
Ti <- seq_len(S); Ni <- (S + 2):nrow(V)
W <- V[Ni, Ti] %*% solve(V[Ti, Ti])
bc <- V[Ni, Ni] - W %*% V[Ti, Ni]

blocks <- list(
  mu1 = list(tip = "mu",        k = 1, anc = "mu_anc[1]",
             rate = "rate_mu[1]"),
  mu2 = list(tip = "mu",        k = 2, anc = "mu_anc[2]",
             rate = "rate_mu[2]"),
  ls1 = list(tip = "log_sigma", k = 1, anc = "log_sigma_anc[1]",
             rate = "rate_log_sigma[1]"),
  ls2 = list(tip = "log_sigma", k = 2, anc = "log_sigma_anc[2]",
             rate = "rate_log_sigma[2]"),
  zr  = list(tip = "z_rho",     k = NULL, anc = "z_rho_anc",
             rate = "rate_rho"))

node_post <- lapply(blocks, function(b) {
  tag <- if (is.null(b$k)) sprintf("%s[%d]", b$tip, 1:S)
         else sprintf("%s[%d,%d]", b$tip, 1:S, b$k)
  tag <- tag[tag %in% names(draws)]
  stopifnot(length(tag) == S)
  yT <- as.matrix(draws[, tag, drop = FALSE])
  anc_d <- draws[[b$anc]]; q <- draws[[b$rate]]
  keep <- is.finite(anc_d) & is.finite(q) & q > 0 &
    rowSums(!is.finite(yT)) == 0
  dd <- which(keep)
  means <- t(vapply(dd, function(d)
    anc_d[d] + W %*% (yT[d, ] - anc_d[d]), numeric(length(Ni))))
  sds <- t(vapply(dd, function(d)
    sqrt(pmax(q[d], 0) * pmax(diag(bc), 0)), numeric(length(Ni))))
  list(mean = colMeans(means),
       sd = sqrt(colMeans(sds^2) + apply(means, 2, var)))
})

saveRDS(list(
  draws = draws,
  summary = summ,
  diagnostics = fit$diagnostic_summary(),
  node_posterior = node_post,
  node_ids = Ni,
  species = stan_data$species,
  masked = stan_data$masked,
  scale_center = mu0, scale_sd = sd0,
  args = list(model = "bm_evolving", chains = 2,
              iter_warmup = 1000, iter_sampling = 1000, seed = 7)
), OUT)
cat("saved:", OUT, "\n")
