## R CMD check results

0 errors | 0 warnings | 2 notes

* New submission. The suggested package cmdstanr is not on CRAN but is
  resolvable via the Additional_repositories field
  (https://stan-dev.r-universe.dev) and all its uses are conditional
  (`requireNamespace`). The package fits Stan models only when cmdstanr
  is installed; MAP estimation via ucminfcpp works without it.

* The non-portable compiler flag note (-mno-omit-leaf-frame-pointer)
  comes from the local R Makeconf on this Ubuntu system, not from any
  flag set by the package.

## Downstream dependencies

None (new submission).
