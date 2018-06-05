#' waldi: Location-adjusted Wald statistic for general model classes
#'
#' The waldi R package provides methods for the computation of the location-adjsuted Wald statistics derived in Di Caterina and Kosmidis <arXiv:1710.11217> for 'glm' and 'brglmFit' abjects, and corresponding methods for location-adjusted Wald confidence intervals.
#'
#' @references
#'
#' Di Caterina C, Kosmidis I (2017). Location-adjusted Wald statistic
#' for scalar parameters *arXiv*, **arXiv:1710.11217**
#'
#' @docType package
#' @name waldi-package
#' @import enrichwith
#' @importFrom stats approx coef confint.default make.link model.matrix predict qnorm spline weights
NULL

#' Generic method for computing location-adjusted Wald statistics for
#' the parameters of general models
#' @param object a fitted model object (e.g. the result of a
#'     \code{\link{glm}} call)
#' @param null a scalar or a vector of the same length as
#'     \code{coef(object)} defining the parameter values for the null
#'     hypotheses. \code{0} by default
#' @param adjust if \code{TRUE} (default) return the location-adjusted
#'     Wald statistic, if \code{FALSE} return the standard Wald
#'     statistic
#' @param which for which parameters should the statistics be
#'     computed? This can be a numeric or character vector.  The
#'     default (\code{NULL}) is to compute statistics for all
#'     parameters
#' @param parallel if \code{TRUE}, the statistics are computed in
#'     parallel using the backend provided to \code{\link{foreach}}
#' @param ... further arguments to be passed to or from other methods
#'
#' @details
#'
#' ADD TECHNICAL DETAILS HERE
#'
#' @seealso  \code{\link{waldi.glm}} \code{\link{waldi.betareg}}
#'
#' @references
#'
#' Di Caterina C, Kosmidis I (2017). Location-adjusted Wald statistic
#' for scalar parameters *arXiv*, **arXiv:1710.11217**
#'
#' @export
waldi <- function(object, null = 0, adjust = TRUE, which = NULL, parallel = FALSE, ...) {
    UseMethod("waldi")
}

if(getRversion() >= "2.15.1")  {
    utils::globalVariables(c("i"))
}
