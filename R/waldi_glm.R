#' Location-adjusted Wald statistics for \code{\link{glm}} objects
#'
#' @inheritParams waldi
#'
#' @param numeric shall numerical derivatives be used for the
#'     computation of the location-adjusted statistics? Default is \code{TRUE}
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
#' waldi(modML, parallel = FALSE)
#'
#' ## Unidentifiability
#' modMLo <- glm(conc ~ log(u)*lot + I(2*log(u)), data = clotting, family = Gamma(link="log"))
#' waldi(modMLo)
#'
#' @export
waldi.glm <- function(object, null = 0, adjust = TRUE, which = NULL, parallel = FALSE,
                       numeric = TRUE, ...) {

    adj_t_numeric <- function(j) {
        u <- numDeriv::grad(kappa, theta, j = j)
        V <- numDeriv::hessian(kappa, theta, j = j)
        a <- -t[j] * u
        a[j] <- 1 + a[j]
        t[j] - sum(a * b, na.rm = TRUE)/ses[j] +
            (sum(F * (tcrossprod(a, u)), na.rm = TRUE)/ses[j] +
             0.5 * t[j] * sum(F * V, na.rm = TRUE))/ses[j]
    }

    ## Can be optimised further to save on some matrix multiplications
    ## by vectorising some of the operations.
    ## Needs work to avoid errors when there are unidentified parameters (see test_analytic.R)
    adj_t_analytic <- function(j) {
        d1mus <- link$mu.eta(etas)
        d2mus <- link$d2mu.deta(etas)
        d3mus <- link$d3mu.deta(etas)
        vars <- family$variance(mus)
        d1vars <- family$d1variance(mus)
        d2vars <- family$d2variance(mus)
        u_analytical <- function(j) {
            e <- d2mus/d1mus
            ## derivative of log variance wrt eta
            l <- d1vars*d1mus/vars
            dw <- w*(2*e - l)*X
            if(no_dispersion) {
            dinfo <- lapply(seq.int(p_all), function(cov) t(X*dw[, cov])%*%X )
            dse <- sapply(seq.int(p_all), function(cov) {
                - 0.5*(diag(F%*%dinfo[[cov]]%*%F))[j]/ses[j]
            })
        }
        else {
            zeta <- - m/phi
            d2as <- family$d2afun(zeta)
            d3as <- family$d3afun(zeta)
            dinfo_beta11 <- lapply(seq.int(p_all), function(cov) t(X*dw[, cov])%*%X/phi)
            dinfo_beta <- lapply(seq.int(p_all), function(cov) cbind(rbind(dinfo_beta11[[cov]], rep(0, p_all)), rep(0, p_all + 1)))
            dse_beta <- sapply(seq.int(p_all), function(cov) {
                - 0.5*(diag(F%*%dinfo_beta[[cov]]%*%F))[j]/ses[j]
            })
            dinfo_phi11 <- - t(X)%*%diag(w)%*%X/phi^2
            dinfo_phi22 <- sum(m^2*d3as)/(2*phi^6) - 2*sum(m^2*d2as)/phi^5
            dinfo_phi <- cbind(rbind(dinfo_phi11, rep(0, p_all)), c(rep(0, p_all), dinfo_phi22))
            dse_phi <- - 0.5*(diag(F%*%dinfo_phi%*%F))[j]/ses[j]
            dse <- as.vector(c(dse_beta, dse_phi))
        }
            dse
        }
        V_analytical <- function(j) {
            e <- d2mus/d1mus
            edash <- (d3mus*d1mus - d2mus^2)/d1mus^2
            ## derivative of log variance wrt eta
            l <- d1vars*d1mus/vars
            ldash <- ((d1vars*d2mus + d2vars*d1mus^2)*vars - d1vars^2*d1mus^2)/vars^2
            if(no_dispersion) {
                dinfo <- function(u) {
                    dw <- w*(2*e - l)*X[, u]
                t(X*dw)%*%X
                }
                d2info <- function(u, v) {
                    d2w <- w*(2*e - l)^2*X[, u]*X[, v] +
                        w*(2*edash - ldash)*X[, u]*X[, v]
                    t(X*d2w)%*%X
                }
                der2 <- function(u, v) {
                    - 0.25 * diag(F%*%dinfo(u)%*%F)[j] * diag(F%*%dinfo(v)%*%F)[j] / ses[j]^3 +
                        0.5 * diag(F%*%dinfo(u)%*%F%*%dinfo(v)%*%F)[j] / ses[j] +
                        0.5*diag(F%*%dinfo(v)%*%F%*%dinfo(u)%*%F)[j]/ses[j] -
                        0.5*diag(F%*%d2info(u, v)%*%F)[j]/ses[j]
                }
                d2se <- outer(seq.int(p_all), seq.int(p_all), Vectorize(der2))
            }
            else {
                zeta <- - m/phi
                d2as <- family$d2afun(zeta)
                d3as <- family$d3afun(zeta)
                d4as <- family$d4afun(zeta)
                dinfo_beta11 <- function(u) {
                    dw <- w*(2*e - l)*X[, u]
                    t(X * dw)%*%X/phi
                }
                dinfo_beta <- function(u) {
                    cbind(rbind(dinfo_beta11(u), rep(0, p_all)), rep(0, p_all + 1))
                }
                d2info_beta11 <- function(u, v) {
                    d2w <- w*(2*e - l)^2*X[, u]*X[, v] + w*(2*edash - ldash)*X[, u]*X[, v]
                    t(X*d2w)%*%X/phi
                }
                d2info_beta <- function(u, v) {
                    cbind(rbind(d2info_beta11(u, v), rep(0, p_all)), rep(0, p_all + 1))
                }
                d2info_betaphi11 <- function(u) {
                    dw <- w*(2*e - l)*X[, u]
                    - t(X * dw)%*%X/phi^2
                }
                d2info_betaphi <- function(u) {
                    cbind(rbind(d2info_betaphi11(u), rep(0, p_all)), rep(0, p_all + 1))
                }
                dinfo_phi11 <- - t(X)%*%diag(w)%*%X/phi^2
                dinfo_phi22 <- sum(m^2*d3as)/(2*phi^6) - 2*sum(m^2*d2as)/phi^5
                dinfo_phi <- cbind(rbind(dinfo_phi11, rep(0, p_all)), c(rep(0, p_all), dinfo_phi22))
                d2info_phi11 <- 2*t(X)%*%diag(w)%*%X/phi^3
                d2info_phi22 <- 10*sum(m^2*d2as)/phi^6 - 5*sum(m^2*d3as)/phi^7 +
                    0.5*sum(m^2*d4as)/phi^8
                d2info_phi <- cbind(rbind(d2info_phi11, rep(0, p_all)), c(rep(0, p_all),
                                                                          d2info_phi22))
                dinfo <- function(u) {
                    if(is.element(u, 1:p_all)) dinfo <- dinfo_beta(u)
                    if(u == p_all + 1) dinfo <- dinfo_phi
                    return(dinfo)
                }
                d2info <- function(u, v) {
                    if((is.element(u, 1:p_all)) & (is.element(v, 1:p_all))) {
                        d2info <- d2info_beta(u, v)
                    }
                    if((is.element(u, 1:p_all)) & (v == p_all + 1)) {
                        d2info <- d2info_betaphi(u)
                    }
                    if((is.element(v, 1:p_all)) & (u == p_all + 1)) {
                        d2info <- d2info_betaphi(v)
                    }
                    if((u == p_all + 1) & (v == p_all + 1)) {
                        d2info <- d2info_phi
                    }
                    return(d2info)
                }
                der2 <- function(u, v) {
                    - 0.25*diag(F%*%dinfo(u)%*%F)[j] * diag(F%*%dinfo(v)%*%F)[j]/ses[j]^3 +
                        0.5*diag(F%*%dinfo(u)%*%F%*%dinfo(v)%*%F)[j]/ses[j] +
                        0.5*diag(F%*%dinfo(v)%*%F%*%dinfo(u)%*%F)[j]/ses[j] -
                        0.5*diag(F%*%d2info(u, v)%*%F)[j]/ses[j]
                }
                d2se <- outer(seq.int(p_all + 1), seq.int(p_all + 1), Vectorize(der2))
            }
            d2se
        }
        a <- -t[j] * u_analytical(j)
        a[j] <- 1 + a[j]
        t[j] - sum(a * b, na.rm = TRUE)/ses[j] +
            (sum(F * (tcrossprod(a, u_analytical(j))), na.rm = TRUE)/ses[j] +
             0.5 * t[j] * sum(F * V_analytical(j), na.rm = TRUE))/ses[j]
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
    }
    else {
        F_inds <- c(beta_ind[!aliased], p_all + 1)
        F <- matrix(NA, p_all + 1, p_all + 1)
    }
    F[F_inds, F_inds] <- solve(info(beta, phi, type = "expected")[F_inds, F_inds])
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
    b <- object$auxiliary_functions$bias(beta, phi)
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

    if (numeric) {
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


