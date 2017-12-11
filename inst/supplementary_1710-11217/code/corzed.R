## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE!  Provided "as is",
## licensed under GPL2 or higher.  NO WARRANTY OF FITNESS FOR ANY
## PURPOSE!
##
## Ioannis Kosmidis [aut], i.kosmidis@ucl.ac.uk
## Claudia Di Caterina [ctb], dicaterina@stat.unipd.it
##
## 26 October 2017
##
## This code still has lots of rough edges, some of which are
## indicated by the embedded comments.
##
## The intention is to include a more polished implementation as part
## or a full R package


#' Location-adjusted Wald statistics
#'
#' @param object an object of class \code{\link{glm}}
#' @param null a numeric vector of the individual null hypotheses to
#'     be tested. Default is that each parameter is zero
#' @param correction should the correction be applied? Default is
#'     \code{TRUE}. If \code{TRUE} then the location-adjusted z
#'     statistic is computed
#' @param ncores how many computing nodes should be used for the
#'     calculation of the adjusted Wald statistics (one node per
#'     parameter). Default is 1
#' @param ... placeholder for extra arguments. Currently not used
#'
#'
#' @author Ioannis Kosmidis [aut, cre] \email{i.kosmidis@ucl.ac.uk},
#'     Claudia Di Caterina [ctb] \email{dicaterina@stat.unipd.it}
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
#' corzed(modML)
#'
#' ## Aliasing
#' modMLo <- glm(conc ~ log(u)*lot + I(2*log(u)), data = clotting, family = Gamma(link="log"))
#' corzed(modMLo)
#'
#' @import enrichwith
#' @import parallel
corzed.glm <- function(object, null = 0, correction = TRUE, ncores = 1, what = NULL,
                       analytical = FALSE, ...) {
    object <- enrich(object, with = c("auxiliary functions", "mle of dispersion"))
    beta <- coef(object, model = "mean")
    phi <- coef(object, model = "dispersion")
    info <- object$auxiliary_functions$information
    no_dispersion <- object$family$family %in% c("poisson", "binomial")
    ## Robustify against aliasing
    pivots <- object$qr$pivot
    p <- object$rank
    p_all <- length(beta)
    beta_ind <- pivots[seq.int(p)]
    ## Compute inverse Fisher information, standard errors and t values
    if (no_dispersion) {
        F_inds <- beta_ind
        F <- matrix(NA, p_all, p_all)
    }
    else {
        F_inds <- c(beta_ind, p_all + 1)
        F <- matrix(NA, p_all + 1, p_all + 1)
    }
    F[F_inds, F_inds] <- solve(info(beta, phi, type = "expected")[F_inds, F_inds])
    ses <- sqrt(diag(F))[1L:p_all]
    t <- (beta - null)/ses
    ## If no correction return t
    if (!correction) {
        if (is.null(what)) {
            return(t)
        }
        else {
            return(t[what])
        }
    }
    ## Otherwise continue to the computation of the adjustment (need
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
    adj_t <- function(j) {
        ## u and V can be replaced with analytical derivatives if necessary
        u <- numDeriv::grad(kappa, theta, j = j)
        V <- numDeriv::hessian(kappa, theta, j = j)
        a <- -t[j] * u
        a[j] <- 1 + a[j]
        t[j] - sum(a * b, na.rm = TRUE)/ses[j] +
            (sum(F * (tcrossprod(a, u)), na.rm = TRUE)/ses[j] +
             0.5 * t[j] * sum(F * V, na.rm = TRUE))/ses[j]
    }
    adj_t_analytical <- function(j) {
      # ADD the analytical adjustment here
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
            - 0.25*diag(F%*%dinfo(u)%*%F)[j]*diag(F%*%dinfo(v)%*%F)[j]/ses[j]^3 +
              0.5*diag(F%*%dinfo(u)%*%F%*%dinfo(v)%*%F)[j]/ses[j] +
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

    ## Adjust the t statistics
    if (is.null(what)) {
      if(analytical) {
        t[beta_ind] <- parallel::mcmapply(adj_t_analytical, beta_ind, mc.cores = ncores)
      }
      else {
        t[beta_ind] <- parallel::mcmapply(adj_t, beta_ind, mc.cores = ncores)
      }
      return(t)
    }
    else {
        if (what > 0 & what <= p_all) {
          if(analytical) return(adj_t_analytical(what))
          else return(adj_t(what))
        }
        else {
            stop("invalid what")
        }
    }
}


corzed.betareg <- function(object, null = 0, correction = TRUE, ncores = 1,
                           what = NULL, ...) {
    if (object$control$hessian) {
        object <- update(object, hessian = FALSE, start = coef(object))
    }
    object <- enrich(object, with = "auxiliary functions")
    theta <- coef(object, model = "full")
    p_all <- length(theta)
    info <- object$auxiliary_functions$information
    F <- solve(info(theta))
    ses <- sqrt(diag(F))[1L:p_all]
    t <- (theta - null)/ses
    ## If no correction return t
    ## Otherwise continue to the computation of the adjustment (need
    ## bias at the mle and information function)
    b <- object$auxiliary_functions$bias(theta)
    kappa <- function(par, j) {
        inv_i <- solve(info(par))
        sqrt(inv_i[j, j])
    }
    adj_t <- function(j) {
        ## u and V can be replaced with analytical derivatives if necessary
        u <- numDeriv::grad(kappa, theta, j = j)
        V <- numDeriv::hessian(kappa, theta, j = j)
        a <- -t[j] * u
        a[j] <- 1 + a[j]
        t[j] - sum(a * b, na.rm = TRUE)/ses[j] +
            (sum(F * (tcrossprod(a, u)), na.rm = TRUE)/ses[j] +
             0.5 * t[j] * sum(F * V, na.rm = TRUE))/ses[j]
    }
    ## Adjust the t statistics
    if (is.null(what)) {
        if (correction) {
            return(parallel::mcmapply(adj_t, 1L:p_all, mc.cores = ncores))
        }
        else {
            t
        }
    }
    else {
        if (what > 0 & what <= p_all) {
            if (correction) {
                adj_t(what)
            }
            else {
                t[what]
            }
        }
        else {
            stop("invalid what")
        }
    }
}


corzed <- function(object, ...) {
    UseMethod("corzed")
}

corzed_ci <- function(object, level = 0.95, correction = TRUE, ncores = 1, what = NULL, length = 20, return_values = FALSE) {
    coefs <- coef(object)
    ses <- sqrt(diag(vcov(object)))
    low <- coefs - 5 * ses
    upp <- coefs + 5 * ses
    p_all <- length(coefs)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    cutoff <- qnorm(a)
    pct <- paste(round(100 * a, 1), "%")
    ci <- function(j) {
        stat <- function(b) {
            corzed(object, what = j, correction = correction, null = b)
        }
        bs <- seq(low[j], upp[j], length = length)
        vals <- sapply(bs, stat)
        if (return_values) {
            vals
        }
        else {
            sp <- spline(x = bs, y = vals)
            approx(sp$y, sp$x, xout = -cutoff)$y
        }
    }
    if (is.null(what)) {
        out <- parallel::mcmapply(ci, 1L:p_all, mc.cores = ncores)
        if (!return_values) {
            out <- t(simplify2array(out))
            rownames(out) <- names(coefs)
            colnames(out) <- pct
        }
    }
    else {
        if (what > 0 & what <= p_all) {
            out <- ci(what)
            if (!return_values) {
                out <- matrix(out, ncol = 2)
                rownames(out) <- names(coefs)[what]
                colnames(out) <- pct
            }
        }
        else {
            stop("invalid what")
        }
    }
    out
}
