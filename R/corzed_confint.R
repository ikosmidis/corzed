#' Confidence intervals by inversion of the location-adjusted Wald statistic

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
