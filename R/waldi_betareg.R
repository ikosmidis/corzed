#' Location-adjusted Wald statistics for \code{\link{betareg}} objects
#'
#' @inheritParams waldi
#'
#' @param numerical shall numerical derivatives be used for the
#'     computation of the location-adjusted statistics? Default is \code{TRUE}. Currently unused.
#'
#' @seealso \code{\link[betareg] {betareg}}; \code{\link[betareg]{summary.betareg}}
#'
#' @examples
#'
#' @export
waldi.betareg <- function(object, null = 0, adjust = TRUE, which = NULL, parallel = FALSE, numerical = TRUE, ...) {

    adj_t_numeric <- function(j) {
        u <- numDeriv::grad(kappa, theta, j = j)
        V <- numDeriv::hessian(kappa, theta, j = j)
        a <- -t[j] * u
        a[j] <- 1 + a[j]
        t[j] - sum(a * b, na.rm = TRUE)/ses[j] +
            (sum(F * (tcrossprod(a, u)), na.rm = TRUE)/ses[j] +
             0.5 * t[j] * sum(F * V, na.rm = TRUE))/ses[j]
    }

    adj_t_analytic <- function(j) {
    }

    object <- enrich(object, with = c("auxiliary functions"))
    theta <- coef(object, model = "full")
    info <- object$auxiliary_functions$information
    ## Robustify against aliasing
    aliased <- is.na(theta)
    p_all <- length(theta)
    theta_ind <- seq.int(p_all)
    ## Compute inverse Fisher information, standard errors and t values
    F_inds <- c(theta_ind[!aliased], p_all + 1)
    F <- matrix(NA, p_all + 1, p_all + 1)

    F[F_inds, F_inds] <- solve(info(theta, phi, type = "expected")[F_inds, F_inds])
    ses <- sqrt(diag(F))[1L:p_all]
    t <- (theta - null)/ses
    theta_names <- names(theta)

    if (is.null(which)) {
        which <- theta_ind
    }
    if (is.character(which)) {
        if (all(!(which %in% theta_names)))
            stop("invalid variables selected")
        which <- match(which, theta_names)
    }

    ## if no correction return t
    if (!adjust) {
        return(t[which])
    }
    ## otherwise continue to the computation of the adjustment (need
    ## bias at the mle and information function)
    b <- object$auxiliary_functions$bias(theta)

    theta[is.na(theta)] <- 0
    kappa <- function(par, j) {
        inv_i <- matrix(NA, p_all + 1, p_all + 1)
        inv_i[F_inds, F_inds] <- solve(info(par[-(p_all + 1)], par[p_all + 1],
                                            type = "expected")[F_inds, F_inds])
        sqrt(inv_i[j, j])
    }

    if (numerical) {
        adj_t <- adj_t_numeric
    }
    else {
        stop("analytical")
    }

    foreach_object <- eval(as.call(c(list(quote(foreach::foreach), i = which, .combine = "c"))))
    if (parallel) {
        setup_parallel()
        t[which] <- foreach::`%dopar%`(foreach_object, ifelse(aliased[i], NA, adj_t(i)))
    }
    else {
        t[which] <- foreach::`%do%`(foreach_object, ifelse(aliased[i], NA, adj_t(i)))
    }
    t[which]
}
