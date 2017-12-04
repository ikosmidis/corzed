#' corzed: Location-adjusted Wald statistic for general model classes
#'
#' ADD PACKAGE DESCRIPTION HERE
#'
#' @details
#'
#' ANY OTHER DETAILS?
#'
#' @author Ioannis Kosmidis \email{i.kosmidis@ucl.ac.uk},
#'         Claudia Di Caterina \email{dicaterina@stat.unipd.it}
#'
#' @references
#'
#' Di Caterina C, Kosmidis I (2017). Location-adjusted Wald statistic
#' for scalar parameters *arXiv*, **arXiv:1710.11217**
#'
#' @docType package
#' @name corzed-package
#' @import enrichwith
NULL

#' Generic method for computing location-adjusted Wald statistics for
#' the parameters of general models
#' @param object a fitted model object (e.g. the result of a
#'     \code{\link{glm}} call)
#' @param null a scalar or a vector of the same length as
#'     \code{coef(object)} defining the parameter values for the null
#'     hypotheses. \code{0} by default
#' @param adjust if \code{TRUE} (default) compute the
#'     location-adjusted Wald statistic is returned, if \code{FALSE}
#'     return the standard Wald statistic
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
#' \code{corzed} stands for "cor"rected "zed"-statistic.
#'
#' @author Ioannis Kosmidis \email{i.kosmidis@ucl.ac.uk},
#'         Claudia Di Caterina \email{dicaterina@stat.unipd.it}
#'
#' @seealso \code{\link{corzed.glm}}
#'
#' @references
#'
#' Di Caterina C, Kosmidis I (2017). Location-adjusted Wald statistic
#' for scalar parameters *arXiv*, **arXiv:1710.11217**
#'
#' @export
corzed <- function(object, null = 0, adjust = TRUE, which = NULL, parallel = FALSE, ...) {
    UseMethod("corzed")
}
