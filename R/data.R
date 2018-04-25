#' Crying babies
#'
#' 18 matched pairs of binomial observations in disaggregated view
#' (one row per baby) from a study on the effect of lulling on the
#' crying of babies. Matching is per day and each day pair in the
#' orginial data consists of the number of babies not crying out of a
#' fixed number of control babies, and the outcome of lulling on a
#' single child. A total of 143 babies were involved in the
#' experiment. Interest is in testing the effect of lulling on the
#' crying of children
#'
#' \itemize{
#'   \item lull; a factor with levels \code{yes} and \code{no}
#'   \item day; a factor with 18 levels, one for each day of the experiment
#'   \item y; \code{0} if the baby cried, and \code{1} otherwise
#' }
#'
#' @source
#'
#' Cox, D. R. (1970) _Analysis of Binary Data_ (page 61).  London:
#' Chapman \& Hall.
#'
#' Davison, A. C. (1988) Approximate conditional inference in
#' generalized linear models.  _J. R. Statist. Soc._ B, *50*,
#' 445-461.
#'
#' @docType data
#' @keywords datasets
"babies"


