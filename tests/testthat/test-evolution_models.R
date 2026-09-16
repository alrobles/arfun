#' Structural tests for the arfun Stan evolutionary models
#'
#' Verifies that each Stan model file (present in inst/stan/) has the expected
#' blocks, required data declarations, and generated quantities. The tests
#' mirror the style of test-stan_models.R (for the presence-only model) and
#' test-stan_phylo.R from nicherbayes.
#'
#' @noRd
library(testthat)
library(arfun)

test_that("all Stan evolutionary models have required blocks and declarations", {
  models <- c(
    "bm_constant" = "niche_bm_constant_sigma.stan",
    "bm_evolving" = "niche_bm_evolving_sigma.stan",
    "ou_constant" = "niche_ou_constant_sigma.stan",
    "ou_evolving" = "niche_ou_evolving_sigma.stan",
    "null"        = "niche_null_model.stan"
  )

  for (model_name in names(models)) {
    path <- model_path(model_name)

    # Read the entire file as a single string
    content <- paste(readLines(path), collapse = "\n")

    # --- Required blocks ---
    expect_true(
      grepl("\\bdata\\s*\\{", content),
      info = paste("Missing data block in", model_name)
    )
    expect_true(
      grepl("\\bparameters\\s*\\{", content),
      info = paste("Missing parameters block in", model_name)
    )
    expect_true(
      grepl("\\btransformed parameters\\s*\\{", content),
      info = paste("Missing transformed parameters block in", model_name)
    )
    expect_true(
      grepl("\\bmodel\\s*\\{", content),
      info = paste("Missing model block in", model_name)
    )
    expect_true(
      grepl("\\bgenerated quantities\\s*\\{", content),
      info = paste("Missing generated quantities block in", model_name)
    )

    # --- Required data fields ---
    required_data <- c("S", "P", "N_occ_total", "N_m_total", "occ", "env_m",
                       "occ_start", "m_start", "grainsize")
    for (field in required_data) {
      expect_true(
        grepl(paste0(field, "\\b"), content),
        info = paste("Missing data field", field, "in", model_name)
      )
    }

    # Phylogenetic matrix: models with phylogeny need C or T_anc
    phylo_models <- c("bm_constant", "bm_evolving", "ou_constant", "ou_evolving")
    if (model_name %in% phylo_models) {
      expect_true(
        grepl("matrix\\[S,\\s*S\\]\\s+C\\b", content) ||
          grepl("matrix\\[S,\\s*S\\]\\s+T_anc\\b", content),
        info = paste("Missing phylogenetic matrix (C or T_anc) in", model_name)
      )
    }

    # Null model should NOT have phylogenetic matrix
    if (model_name == "null") {
      expect_false(
        grepl("matrix\\[S,\\s*S\\]\\s+C\\b", content),
        info = "Null model should not have phylogenetic C matrix"
      )
      expect_false(
        grepl("matrix\\[S,\\s*S\\]\\s+T_anc\\b", content),
        info = "Null model should not have phylogenetic T_anc matrix"
      )
    }

    # --- Partial functions for reduce_sum ---
    expect_true(
      grepl("partial_sum", content),
      info = paste("Missing partial_sum function in", model_name)
    )

    # --- Generated quantities ---
    expect_true(
      grepl("Sigma\\[s\\]", content),
      info = paste("Missing Sigma[s] in generated quantities of", model_name)
    )
    expect_true(
      grepl("R_corr\\[s", content),
      info = paste("Missing R_corr[s] in generated quantities of", model_name)
    )
    expect_true(
      grepl("log_lik\\[s\\]", content),
      info = paste("Missing log_lik[s] in generated quantities of", model_name)
    )

    # --- Per-species likelihood ---
    expect_true(
      grepl("log_sum_exp", content),
      info = paste("Missing log_sum_exp (Jimenez likelihood) in", model_name)
    )
  }
})

test_that("model_path returns valid paths for all evolutionary models", {
  models <- c("bm_constant", "bm_evolving", "ou_constant", "ou_evolving", "null")

  for (model in models) {
    path <- model_path(model)
    expect_true(file.exists(path),
                info = paste("model_path(", model, ") returned non-existent path:", path))
    expect_true(grepl("\\.stan$", path),
                info = paste("model_path(", model, ") did not return a .stan file"))
  }
})

test_that("model names are case-sensitive and match switch mapping", {
  expected_filenames <- c(
    bm_constant = "niche_bm_constant_sigma.stan",
    bm_evolving = "niche_bm_evolving_sigma.stan",
    ou_constant = "niche_ou_constant_sigma.stan",
    ou_evolving = "niche_ou_evolving_sigma.stan",
    null        = "niche_null_model.stan"
  )

  for (model_name in names(expected_filenames)) {
    path <- model_path(model_name)
    fname <- basename(path)
    expect_equal(fname, expected_filenames[[model_name]],
                 label = paste0("model_path(", model_name, ")"),
                 expected.label = expected_filenames[[model_name]])
  }
})
