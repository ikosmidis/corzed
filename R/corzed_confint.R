#' Confidence intervals by inversion of the location-adjusted Wald statistic
#'
#' @inheritParams corzed
#'
#' @param level the confidence level required. Default is 0.95
#'
#' @param length The length of the grid on which the location-adjusted
#'     statistic is evaluated. Default is 20. See details.
#'
#' @param return_values Return the values of the statistic on the grid
#'     instead of confidence intervals? Default is \code{FALSE}.
#'
#' @seealso \code{\link{summary.glm}}
#'
#' @examples
#' clotting <- data.frame(
#'    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
#'    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
#'    lot = factor(c(rep(1, 9), rep(2, 9))))
#'
#' modML <- glm(conc ~ log(u)*lot, data = clotting, family = Gamma(link="log"))
#' corzed_confint(modML, parallel = FALSE)
#'
#' \dontrun{
#' ## Now do the same in parallel
#' library("foreach")
#' library("doMC")
#' registerDoMC(2)
#' corzed_confint(modML, parallel = TRUE)
#' }
#'
#' @export
corzed_confint <- function(object, level = 0.95, adjust = TRUE, which,
                           parallel = TRUE, numeric = TRUE,
                           length = 20, return_values = FALSE) {

    ci <- function(j) {
        stat <- function(b) {
            corzed(object, null = b, adjust = TRUE, which = j)
        }
        bs <- seq(low[j], upp[j], length = length)
        if (aliased[j]) {
            return(NAout)
        }
        vals <- sapply(bs, stat)
        if (return_values) {
            vals
        }
        else {
            sp <- spline(x = bs, y = vals)
            approx(sp$y, sp$x, xout = -cutoff)$y
        }
    }

    cis <- confint.default(object, level = level)
    len <- apply(cis, 1, diff)
    par_names <- rownames(cis)
    if (!adjust) {
        return(cis[which,])
    }
    low <- cis[, 1] - len/3
    upp <- cis[, 2] + len/3
    names(low) <- names(upp) <- names(par_names) <- par_names
    aliased <- apply(cis, 1, function(x) all(is.na(x)))
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    cutoff <- qnorm(a)
    pct <- paste(round(100 * a, 1), "%")
    out_length <- ifelse(return_values, length, 2)
    NAout <- rep(NA, out_length)
    if (missing(which)) {
        which <- seq.int(length(coef(object)))
    }
    foreach_object <- eval(as.call(c(list(quote(foreach::foreach), i = which, .combine = "rbind", .init = NULL))))
    if (parallel) {
        setup_parallel()
        out <- foreach::`%dopar%`(foreach_object, if (aliased[i]) NAout else ci(i))
    }
    else {
        out <- foreach::`%do%`(foreach_object, if (aliased[i]) NAout else ci(i))
    }
    rownames(out) <- par_names[which]
    colnames(out) <- pct
    out
}
