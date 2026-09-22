#' Keep only the fields a Stan model declares
#'
#' Filters a prepared data list (from
#' \code{\link{prepare_phylo_niche_data}()} or \code{\link{prepare_phylo_data}()})
#' down to exactly the variables declared in the \code{data {}} block of
#' the named Stan model. Prepare functions attach bookkeeping fields
#' (\code{species}, \code{masked}, \code{eta}, \code{T_anc}, \dots) that are
#' useful in R but rejected by \pkg{cmdstanr}'s type validation; this
#' helper drops them so the list can be passed straight to
#' \code{fit_evolution()} or \code{CmdStanModel$sample()}.
#'
#' @param data A named list of Stan data, as returned by the prepare
#'   functions.
#' @param model Character string specifying the model whose data block is
#'   parsed (e.g. \code{"bm_evolving"}). See \code{\link{model_path}()}.
#'
#' @return A named list containing only the declared fields, in the
#'   original order.
#'
#' @details The parser reads the bundled \code{.stan} file at install time
#'   of each call, so the field set always matches the shipped model --
#'   adding a data variable to a Stan file automatically whitelists it.
#'
#' @examples
#' \dontrun{
#' stan_data <- prepare_phylo_niche_data(env_occ_list, env_m_list, tree)
#' stan_data <- stan_only(stan_data, "bm_evolving")
#' fit <- fit_evolution(stan_data, model = "bm_evolving")
#' }
#'
#' @export
stan_only <- function(data, model) {
  fields <- .stan_data_fields(model)
  missing <- setdiff(fields, names(data))
  if (length(missing)) {
    warning("declared fields absent from 'data': ",
            paste(missing, collapse = ", "), call. = FALSE)
  }
  data[intersect(fields, names(data))]
}

# Variable names declared in the `data {}` block of a bundled Stan model.
.stan_data_fields <- function(model) {
  lines <- readLines(model_path(model), warn = FALSE)
  in_block <- FALSE
  fields <- character(0)
  for (ln in lines) {
    ln <- sub("//.*$", "", ln)                    # strip comments
    if (!in_block) {
      if (grepl("^\\s*data\\s*\\{", ln)) in_block <- TRUE
      next
    }
    if (grepl("\\}", ln)) break
    # declaration lines end with "name;" (possibly "name[n];")
    m <- regmatches(ln,
                    regexpr("[A-Za-z_][A-Za-z0-9_]*(\\[[^]]*\\])?\\s*;\\s*$",
                            ln))
    if (length(m)) {
      fields <- c(fields, sub("\\s*;\\s*$", "",
                              sub("\\[[^]]*\\]", "", m)))
    }
  }
  fields
}
