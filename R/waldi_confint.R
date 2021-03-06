#' Confidence intervals by inversion of the location-adjusted Wald statistic
#'
#' @inheritParams waldi
#' @inheritParams waldi.glm
#'
#' @param level the confidence level required. Default is 0.95
#'
#' @param adjust if \code{TRUE} (default) use the location-adjusted
#'     Wald statistic, if \code{FALSE} use the standard Wald
#'     statistic
#'
#' @param length The length of the grid on which the location-adjusted
#'     statistic is evaluated. Default is 20. See details.
#'
#' @param return_values Return the values of the statistic on the grid
#'     instead of confidence intervals? Default is \code{FALSE}
#'
#' @param quantiles A \code{length(which)} (or
#'     \code{length(coef(object))} if \code{which} is not specified)
#'     times 2 matrix with user-supplied quantiles to be used for the
#'     inversion of the statistics (e.g. for bootstrap studentized
#'     intervals). Default is \code{NULL}. If specified, then the
#'     value for \code{level} is ignored.
#'
#' @seealso \code{\link{summary.glm}}
#'
#' @examples
#' clotting <- data.frame(
#'    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
#'    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
#'    lot = factor(c(rep(1, 9), rep(2, 9))))
#'
#' modML <- glm(conc ~ log(u)*lot, data = clotting,
#'              family = Gamma(link="log"))
#' waldi_confint(modML, parallel = FALSE)
#'
#' \dontrun{
#' ## Now do the same in parallel using 4 cores
#' library("foreach")
#' library("doMC")
#' registerDoMC(4)
#'
#' ## Differences between the Wald and location-adjusted Wald statistic
#' data("babies", package = "waldi")
#' babies_ml0 <- glm(formula = y ~ day + lull - 1, family = binomial,
#'                   data = babies)
#' out_waldi <- waldi_confint(babies_ml0, adjust = TRUE, which = "lullyes",
#'                              parallel = TRUE, return_values = TRUE)
#' out_wald <- waldi_confint(babies_ml0, adjust = FALSE, which = "lullyes",
#'                            parallel = TRUE, return_values = TRUE)
#' ## Statistics and critical values for 95\% 2-sided intervals
#' with(out_waldi, plot(grid, value, type = "l", col = "red"))
#' with(out_wald, points(grid, value, type = "l", col = "blue"))
#' abline(a = qnorm(0.975), b = 0, lty = 2)
#' abline(a = qnorm(0.025), b = 0, lty = 2)
#' legend(x = "topright", legend = c(expression(t), expression(t^'*')),
#'        col = c("blue", "red"), lty = 1)
#' }
#'
#' @export
waldi_confint <- function(object, level = 0.95, quantiles = NULL,
                          adjust = TRUE, which,
                          parallel = FALSE, numerical = TRUE,
                          length = 5, return_values = FALSE) {
    ci <- function(j) {
        stat <- function(b) {
            waldi(object, null = b, adjust = adjust, which = j, numerical = numerical)
        }
        bs <- seq(low[j], upp[j], length = length)
        if (aliased[j]) {
            return(NAout)
        }
        vals <- sapply(bs, stat)
        if (return_values) {
            data.frame(grid = bs, value = vals, parameter = j)
        }
        else {
            sp <- spline(x = bs, y = vals)
            approx(sp$y, sp$x, xout = -cutoff[i, ])$y
        }
    }

    cis <- confint.default(object, level = level + (1 - level)/2)
    len <- apply(cis, 1, diff)
    par_names <- rownames(cis)
    if (missing(which)) {
        which <- seq.int(length(coef(object)))
    }
    npar <- length(which)
    ## if (!adjust & !return_values) {
    ##     return(cis[which,])
    ## }
    low <- cis[, 1] - len/2
    upp <- cis[, 2] + len/2

    names(low) <- names(upp) <- names(par_names) <- par_names
    aliased <- apply(cis, 1, function(x) all(is.na(x)))
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    if (no_quantiles <- is.null(quantiles)) {
        cutoff <- matrix(qnorm(a), npar, 2, byrow = TRUE)
        rownames(cutoff) <- par_names[which]
    }
    else {
        if (!isTRUE(nrow(quantiles) == npar) | !isTRUE(ncol(quantiles) == 2)) {
            stop("quantiles needs to be a ", npar, " times", 2, " matrix.")
        }
        cutoff <- quantiles
    }
    pct <- paste(round(100 * a, 1), "%")
    out_length <- ifelse(return_values, length, 2)
    NAout <- rep(NA, out_length)
    foreach_object <- eval(as.call(c(list(quote(foreach::foreach), i = which, .combine = "rbind", .init = NULL))))
    if (parallel) {
        setup_parallel()
        out <- foreach::`%dopar%`(foreach_object, if (aliased[i]) NAout else ci(i))
    }
    else {
        out <- foreach::`%do%`(foreach_object, if (aliased[i]) NAout else ci(i))
    }
    if (!return_values) {
        rownames(out) <- par_names[which]
        colnames(out) <- if (no_quantiles) pct else c("low", "upper")
    }
    out
}

