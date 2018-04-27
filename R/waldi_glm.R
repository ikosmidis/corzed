#' Location-adjusted Wald statistics for \code{\link{glm}} and
#' \code{\link[brglm2]{brglmFit}} objects
#'
#' @inheritParams waldi
#'
#' @param numerical shall numerical derivatives be used for the
#'     computation of the location-adjusted statistics? Default is \code{TRUE}
#'
#' @seealso \code{\link{summary.glm}}
#'
#' @examples
#'
#' data("babies", package = "waldi")
#' babies_ml <- glm(formula = y ~ day + lull - 1, family = binomial, data = babies)
#' waldi(babies_ml, numerical = FALSE)
#'
#'
#' clotting <- data.frame(
#'    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
#'    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
#'    lot = factor(c(rep(1, 9), rep(2, 9))))
#'
#' modML <- glm(conc ~ log(u)*lot, data = clotting, family = Gamma(link="log"))
#' waldi(modML, parallel = FALSE)
#'
#' ## Unidentifiability
#' modMLo <- glm(conc ~ log(u)*lot + I(2*log(u)), data = clotting, family = Gamma(link="log"))
#' waldi(modMLo)
#'
#' @export
waldi.glm <- function(object, null = 0, adjust = TRUE, which = NULL, parallel = FALSE,
                       numerical = TRUE, ...) {

    br <- inherits(object, "brglmFit")

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
        Fj <- F[j, , drop = TRUE]
        prods <- sapply(d1_info, function(x) x %*% Fj)
        v2 <- function(u, v) {
             prods[, u] %*% F %*% prods[, v] - 0.5 * Fj %*% d2_info(u, v) %*% Fj
        }
        u <- -0.5 * drop(Fj %*% prods) / ses[j]
        V <- -tcrossprod(u) / ses[j] + outer(p_inds, p_inds, Vectorize(v2)) / ses[j]
        a <- -t[j] * u
        a[j] <- 1 + a[j]
        t[j] - sum(a * b, na.rm = TRUE)/ses[j] +
            (sum(F * (tcrossprod(a, u)), na.rm = TRUE)/ses[j] +
             0.5 * t[j] * sum(F * V, na.rm = TRUE))/ses[j]
    }

    object <- enrich(object, with = c("auxiliary functions", "mle of dispersion"))
    beta <- coef(object, model = "mean")
    phi <- coef(object, model = "dispersion")
    info <- object$auxiliary_functions$information
    no_dispersion <- object$family$family %in% c("poisson", "binomial")
    ## Robustify against aliasing
    aliased <- is.na(beta)
    p_all <- length(beta)
    beta_ind <- seq.int(p_all)
    ## Compute inverse Fisher information, standard errors and t values
    if (no_dispersion) {
        F_inds <- beta_ind[!aliased]
        F <- matrix(NA, p_all, p_all)
        p_inds <- seq.int(p_all)
    }
    else {
        F_inds <- c(beta_ind[!aliased], p_all + 1)
        F <- matrix(NA, p_all + 1, p_all + 1)
        p_inds <- seq.int(p_all + 1)
    }

    ## We can save the inversion here by extracting directly from the object. But let's keep it for now
    info0 <- info(beta, phi, type = "expected")[F_inds, F_inds]
    F[F_inds, F_inds] <- solve(info0)
    ses <- sqrt(diag(F))[1L:p_all]
    t <- (beta - null)/ses
        beta_names <- names(beta)

    if (is.null(which)) {
        which <- beta_ind
    }
    if (is.character(which)) {
        if (all(!(which %in% beta_names)))
            stop("invalid variables selected")
        which <- match(which, beta_names)
    }

    ## if no correction return t
    if (!adjust) {
        return(t[which])
    }
    ## otherwise continue to the computation of the adjustment (need
    ## bias at the mle and information function)
    b <- if (br) 0 else object$auxiliary_functions$bias(beta, phi)
    if (no_dispersion) {
        theta <- beta
        theta[is.na(theta)] <- 0
        kappa <- function(par, j) {
            inv_i <- matrix(NA, p_all, p_all)
            inv_i[F_inds, F_inds] <- solve(info(par, 1, type = "expected")[F_inds, F_inds])
            sqrt(inv_i[j, j])
        }
    }
    else {
        theta <- c(beta, phi)
        theta[is.na(theta)] <- 0
        kappa <- function(par, j) {
            inv_i <- matrix(NA, p_all + 1, p_all + 1)
            inv_i[F_inds, F_inds] <- solve(info(par[-(p_all + 1)], par[p_all + 1],
                                                type = "expected")[F_inds, F_inds])
            sqrt(inv_i[j, j])
        }
    }

    if (numerical) {
        adj_t <- adj_t_numeric
    }
    else {
        X <- model.matrix(object)
        family <- enrich(object$family)
        link <- enrich(make.link(object$family$link))
        etas <- predict(object, type = "link")
        mus <- predict(object, type = "response")
        w <- weights(object, type = "working")
        m <- weights(object, type = "prior")
        d1mus <- link$mu.eta(etas)
        d2mus <- link$d2mu.deta(etas)
        d3mus <- link$d3mu.deta(etas)
        vars <- family$variance(mus)
        d1vars <- family$d1variance(mus)
        d2vars <- family$d2variance(mus)
        e <- d2mus / d1mus
        edash <- (d3mus * d1mus - d2mus^2) / d1mus^2
        l <- d1vars * d1mus / vars
        ldash <- ( (d1vars * d2mus + d2vars * d1mus^2) * vars - d1vars^2 * d1mus^2) / vars^2
        dw <- w * (2 * e - l) * X
        ## First derivatives
        d1_info_db_b <- lapply(seq.int(p_all), function(param) t(X * dw[, param]) %*% X / phi)
        d2_info_dbb_b <- function(u, v) t(X * w * ((2 * e - l)^2 + (2 * edash - ldash)) * X[, u] * X[, v]) %*% X / phi
        if (no_dispersion) {
            d1_info <- d1_info_db_b
            d2_info <- d2_info_dbb_b
        }
        else {
            zeta <- - m / phi
            d2as <- family$d2afun(zeta)
            d3as <- family$d3afun(zeta)
            d4as <- family$d4afun(zeta)
            d1_info <- lapply(d1_info_db_b, function(x) rbind(cbind(x, 0), 0))
            d1_info_dp <- -info0 / phi
            d1_info_dp[p_all + 1, p_all + 1] <- sum(m^2 * d3as)/(2 * phi^6) - 2 * sum(m^2 * d2as) / phi^5
            d1_info <- c(d1_info, list(d1_info_dp))
            d2_info_dpp <- 2 * info0 / phi^2
            d2_info_dpp[p_all + 1, p_all + 1] <- 10 * sum(m^2 * d2as) / phi^6 - 5 * sum(m^2 * d3as) / phi^7 + 0.5 * sum(m^2 * d4as) / phi^8
            d2_info <- function(u, v) {
                if (u <= p_all & v <= p_all) {
                    out <- rbind(cbind(d2_info_dbb_b(u, v), 0), 0)
                }
                if (u == p_all + 1 & v == p_all + 1) {
                    out <- d2_info_dpp
                }
                if (u <= p_all & v == p_all + 1) {
                    out <- -d1_info[[u]] / phi
                }
                if (v <= p_all & u == p_all + 1) {
                    out <- -d1_info[[v]] / phi
                }
                out
            }
        }

        ## Second derivatives

        adj_t <- adj_t_analytic

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
