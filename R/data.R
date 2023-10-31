#' Simulated Metabolite Data
#'
#' A 15x100x4 array containing 20 simulated target metabolites and 80 simulated
#' noise metabolites across 15 observations and 4 time points. Intended for use
#' in `prolong` package plotting functions or as predictors for `sim_outcome` in
#' the `prolong()` function.
#'
#' @format ## `sim_metabs`
#' A 3-D numeric array with n = 15 rows, p = 100 columns, and t = 4 slices
#'
#' @usage data(sim_metabs)

"sim_metabs"

#' Simulated Clinical Outcome Data
#'
#' A 15x4 matrix containing 20 simulated clinical outcome values. Intended for
#' use as the response variable with predictors from `sim_metabs` in
#' the `prolong()` function.
#'
#' @format ## `sim_outcome`
#' A numeric matrix with n = 15 rows and t = 4 columns
#'
#' @usage data(sim_metabs)



"sim_outcome"
