#' Confidence-region projection for integrated Emax responses
#'
#' RcppArmadillo implementation for fixed replicated doses, independent Gaussian
#' errors, and a common unknown variance. See [emax_region()], [emax_project()],
#' and the getting-started vignette for assumptions and numerical diagnostics.
#' @useDynLib emaxProjection, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"
