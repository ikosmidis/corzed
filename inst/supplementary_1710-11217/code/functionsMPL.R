## Distributed as part of the supplementary material for the the paper
## ``Improving the accuracy of likelihood-based inference in
## meta-analysis and meta-regression''
##
## Authors: Ioannis Kosmidis, Annamaria Guolo, Cristiano Varin
## Date: 15 October 2016
## Licence: GPL 2 or greater

## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE!  Provided "as is".
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!

## Returns maximum likelihood (ML), maximum penalized likelihood
## (MPL), bias-corrected (BC) and DerSimonian & Laird (DL) estimates
## for the parameters of the a meta-regression model (see Sections
## 3.1, 3.2 and 6 of the main text)
##
## Arguments:
## object: metaLik object
## method: the estimation method to be used (default is "all")
BiasFit <- function(object,
                    method = c("all", "MPL", "BC", "DL"), ...) {
    method <- match.arg(method)
    parstau2 <- object$mle
    p <- ncol(object$X)
    tau2 <- parstau2[p + 1]
    zeta <- tau2
    pars <- c(object$mle[seq.int(p)], zeta)
    b <- with(object, bias(pars, X, sigma2, y))
    BC <- MPL <- DL <- rep(NA, p + 1)
    if (method == "BC" | method == "all") {
        BC <- pars - c(numeric(p), b)
    }
    if (method == "MPL" | method == "all") {
        MPL <- optim(par = pars,
                    fn = penloglik,
                    gr = penscores,
                    X = object$X, y = object$y, sigma2 = object$sigma2,
                    control = list(fnscale = -1), method = "BFGS")
        MPL <- MPL$par
    }
    if (method == "DL" | method == "all") {
        DL <- c(drop(object$DL), with(object, {tau2 <- object$tau2.DL}))
    }
    ML <- optim(par = pars,
              fn = loglik,
              gr = scores,
              X = object$X, y = object$y, sigma2 = object$sigma2,
              control = list(fnscale = -1), method = "BFGS")
  list(ML = ML$par,
       MPL = MPL,
       BC = BC,
       DL = DL)
}


## Performs individual tests for the parameters of a meta-regression
## model (see Section 3.3 of the main text)

## Arguments:
## object: metaLik object
## what: the index of the parameter to test for
## null: the null hypothesis (the alternative is always not the null)
## type: likelihood based or penalized likelihood based tests (defaults to likelihood)?
## optMethod: the optimization method to be used for the profiling
lrtest <- function(object,
                   what = 1,
                   null = 0,
                   type = c("loglik", "penloglik"),
                   optMethod = "Nelder-Mead",
                   alternative = c("two.sided", "less", "greater"),
                   ...) {
    alternative <- match.arg(alternative)
    X <- object$X
    y <- object$y
    p <- ncol(X)
    sigma2 <- object$sigma2
    parstau2 <- pars <- object$mle
    type <- match.arg(type)
    interval <- c(0, 1000)
    IF <- switch(type,
                 "loglik" = loglik,
               "penloglik" = penloglik)
    DIF <- switch(type,
                  "loglik" = scores,
                  "penloglik" = penscores)
    profzeta <- function(zeta,
                         ystar, Xstar, X,
                         penalty = TRUE) {
        tau2 <- zeta
        w <- 1/(sigma2 + tau2)
        beta <- coef(lm.wfit(x = Xstar, y = ystar, w = w))
        res <- ystar - Xstar %*% beta
        W.X <- sqrt(w)*X
        XWX <- crossprod(W.X)
        0.5*sum(log(w)) - 0.5*sum(w*res^2) - penalty*(0.5*log(det(XWX)))
    }
    fun <- function(range) {
        params <- numeric(p)
        profvalues <- resFull$value
        params[what] <- range
        ynew <- y - as.matrix(X[, what]) %*% range
        Xnew <- as.matrix(X[, -what])
        zetaConstr <- optimize(profzeta,
                               interval = interval,
                               ystar = ynew,
                               Xstar = Xnew,
                               X = X,
                               penalty = (type == "penloglik"),
                               maximum = TRUE)$maximum
        if (NCOL(Xnew) > 0) {
            tau2Constr <- zetaConstr
            wConstr <- 1/(sigma2 + tau2Constr)
            params[-what] <- coef(lm.wfit(x = Xnew, y = ynew, w = wConstr))
        }
        profvalues <- IF(c(params, zetaConstr), X = X, y = y,
                         sigma2 = sigma2)
        list(parameterValues = params, profileValues = profvalues)
    }
    ## Get the full maximum (penalized) likelihood estimates
    resFull <- optim(par = pars,
                     fn = IF,
                     gr = DIF,
                     X = X, y = y, sigma2 = sigma2,
                      control = list(fnscale = -1),
                     method = optMethod,
                     hessian = TRUE)
    resRestricted <- fun(null)
    statistic <- sign(resFull$par[what] - null)*sqrt(2*(resFull$value - resRestricted$profileValues))
    if (is.na(statistic)) {
        statistic <- 0
    }
    pvalue <- switch(alternative,
                     "two.sided" = 2*pnorm(-abs(statistic)),
                     "less" = pnorm(statistic, lower.tail = TRUE),
                     "greater" = pnorm(statistic, lower.tail = FALSE))
    list(statistic = statistic,
         pvalue = pvalue)
}

## Implementation Zeng and Lin' double resampling
## it assumes a call to metaLik
##
## Arguments:
## beta0: null hypothesis
## m: metaLik object
## B: number of bootstrap samples (output B^2 bootstrap estimates)
## myseed: seed to be used for random variate generation
## ci: return *also* confidence interval?
## alpha: confidence interval level (only if pval = FALSE)
double.resampling <- function(beta0, m, B = 1000, myseed = 123, ci = FALSE, alpha = 0.05,
                              alternative = c("two.sided", "less", "greater", "all"))
{
    alternative <- match.arg(alternative) ## Currently inactive
    ## take care of the random seed
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        seed.keep <- get(".Random.seed", envir = .GlobalEnv,
                         inherits = FALSE)
        on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
    }
    set.seed(myseed)
    ## extract info
    beta.DL <- c(m$DL) ## assumes meta-analysis
    tau2.DL <- m$tau2.DL
    sigma2 <- m$sigma2
    K <- m$K

    beta.resampled <- matrix(NA, nrow = B, ncol = B)
    ## some recurrent computations
    S1 <- sum(1/sigma2)
    S2 <- sum(1/sigma2^2)
    denom <- S1-S2/S1
    compute.tau2 <- function(betak, betahat){
        eps <- t(betak) - c(betahat)
        tmp <- eps^2 %*% ( 1/sigma2 ) - (K-1)
        pmax(0, tmp/denom)
    }
    compute.betahat <- function(betak, tau2){
        t(betak) %*% (1 / (sigma2 + tau2) ) / sum(1 / (sigma2 + tau2) )
    }

    ## first run
    betak <- replicate(B, beta.DL + rnorm(K) * sqrt(sigma2 + tau2.DL) )
    betahat <- compute.betahat(betak, 0.0)
    tau2hat <- compute.tau2(betak, betahat)

    ## second run
    for(b in 1:B){
        betak <- replicate(B, beta.DL + rnorm(K) * sqrt(sigma2 + tau2hat[b]))
        beta.resampled[b,] <- compute.betahat(betak, tau2hat[b])
    }
    beta.resampled <- c(beta.resampled)
    pval <- switch(alternative,
                   "two.sided" = 2 * min( mean( beta.resampled > beta0 ), mean( beta.resampled < beta0 )),
                   "less" = mean( beta.resampled < beta0 ),
                   "greater" = mean( beta.resampled > beta0 ),
                   "all" = c(two.sided = 2 * min( mean( beta.resampled > beta0 ), mean( beta.resampled < beta0 )),
                             less = mean( beta.resampled < beta0 ),
                             greater = mean( beta.resampled > beta0 )))
    if(ci){
        ci <- switch(alternative,
                     "two.sided" = quantile( beta.resampled, c(alpha / 2, 1 - alpha / 2) ),
                     "less" = c(-Inf, quantile(beta.resampled, 1 - alpha)),
                     "greater" = c(quantile(beta.resampled, alpha), Inf),
                     "all" = list(two.sided = quantile( beta.resampled, c(alpha / 2, 1 - alpha / 2) ),
                                  less = c(-Inf, quantile(beta.resampled, 1 - alpha)),
                                  greater = c(quantile(beta.resampled, alpha), Inf)))

        return( list(pval = pval, ci = ci) )
    }
    else
        return(pval)
}

## Calculates p-values and statistics for a single parameter and
## various alternatives
##
## Arguments:
## y: vector of summary statistics
## X: the model matrix for the fixed effects of tha random effects meta-regression model
## sigma2: vector of summary variances
## what: what parameter should perform_tests test for?
## null: the null hypothesis
## B: the number of bootstrap samples (output B^2 bootstrap estimates)
## silent: should warnings and errors be reported?
## maxiter: maximum number of iterations to be passed to the metatest::metatest function
## seed: seed to be used for random variate generation
perform_tests <- function(y, X, sigma2, what = 1, null = 0, B = 1000, silent = TRUE, maxiter = 1000, seed = 123) {
    require(metaLik)
    require(metatest)
    fit_metaLik <- try(metaLik(y ~ -1 + X, sigma2 = sigma2), silent = silent)
    ## Skovgaard
    test_metaLik <- try(test.metaLik(fit_metaLik, param = what, value = null, print = FALSE), silent = silent)
    stat_skovgaard <- test_metaLik$rskov
    pvalue_g_skovgaard <- pnorm(stat_skovgaard, lower.tail = FALSE)
    pvalue_d_skovgaard <- 2 * pnorm(-abs(stat_skovgaard))
    pvalue_l_skovgaard <- pnorm(stat_skovgaard, lower.tail = TRUE)
    ## Bartlett: The following has no effect if the null is zero. If
    ## the null is non-zero, it will result in correct p-values for
    ## meta-analysis but incorrect for meta-regression.
    test_metatest <- try(metatest(I(y - null*X[,what]) ~ -1 + X, variance = sigma2, npermut = 0,
                                  maxiter = maxiter),
                         silent = silent)
    if (inherits(fit_metaLik, "try-error") | inherits(test_metaLik, "try-error") | inherits(test_metatest, "try-error")) {
        res <- matrix(NA, 8, 4)
    }
    else {
        if (ncol(X) == 1) {
            ## NOTE: Definition of signed root based on bartlett
            ## corrected statistic
            ##
            ## If Bartlett corrected statistic ends up negative set it
            ## to 0 to agree with how metatest computes p-values
            sbart <- test_metatest$bartLLR[what]
            sbart <- ifelse(sbart < 0, 0, sbart)
            stat_bartlett <- sign(fit_metaLik$beta.mle[what] - null)*sqrt(sbart)
            pvalue_g_bartlett <- pnorm(stat_bartlett, lower.tail = FALSE)
            pvalue_d_bartlett <- 2 * pnorm(-abs(stat_bartlett))
            pvalue_l_bartlett <- pnorm(stat_bartlett, lower.tail = TRUE)
        }
        else {
            stat_bartlett <- pvalue_g_bartlett <- pvalue_d_bartlett <- pvalue_l_bartlett <- NA
        }

        ## Likelihood ratio
        test_LR <- try(lrtest(fit_metaLik, what = what, type = "loglik",
                              null = null,
                              optMethod = "BFGS"), silent = silent)
        if (inherits(test_LR, "try-error")) {
            test_LR <- try(lrtest(fit_metaLik, what = what, type = "loglik",
                                  null = null,
                                  optMethod = "Nelder-Mead"), silent = silent)
        }

        stat_LR <- test_LR$statistic
        pvalue_g_LR <- pnorm(stat_LR, lower.tail = FALSE)
        pvalue_d_LR <- 2 * pnorm(-abs(stat_LR))
        pvalue_l_LR <- pnorm(stat_LR, lower.tail = TRUE)

        ## Penalized likelihood ratio
        test_PLR <- try(lrtest(fit_metaLik, what = what, type = "penloglik",
                               null = null,
                               optMethod = "BFGS"), silent = silent)
        if (inherits(test_PLR, "try-error")) {
            test_PLR <- try(lrtest(fit_metaLik, what = what, type = "penloglik",
                                   null = null,
                                   optMethod = "Nelder-Mead"), silent = silent)
        }
        stat_PLR <- test_PLR$statistic
        pvalue_g_PLR <- pnorm(stat_PLR, lower.tail = FALSE)
        pvalue_d_PLR <- 2 * pnorm(-abs(stat_PLR))
        pvalue_l_PLR <- pnorm(stat_PLR, lower.tail = TRUE)


        ## Double resampling
        stat_db <- NA
        if (ncol(X) == 1 & all(X == 1)) {
            pvalue_db <- double.resampling(beta0 = null, m = fit_metaLik, B = B, alternative = "all", ci = FALSE, myseed = seed)
            pvalue_d_db <- pvalue_db["two.sided"]
            pvalue_l_db <- pvalue_db["less"]
            pvalue_g_db <- pvalue_db["greater"]
        }
        else {
            pvalue_d_db <- pvalue_g_db <- pvalue_l_db <- NA
        }

        ## Wald based on ML
        stat_waldML <- (fit_metaLik$beta.mle[what] - null)/sqrt(fit_metaLik$vcov[what, what])
        pvalue_g_waldML <- pnorm(stat_waldML, lower.tail = FALSE)
        pvalue_d_waldML <- 2 * pnorm(-abs(stat_waldML))
        pvalue_l_waldML <- pnorm(stat_waldML, lower.tail = TRUE)

        ## Wald based on DL
        stat_DL <- (fit_metaLik$DL[what] - null)/sqrt(fit_metaLik$vcov.DL[what, what])
        pvalue_g_DL <- pnorm(stat_DL, lower.tail = FALSE)
        pvalue_d_DL <- 2 * pnorm(-abs(stat_DL))
        pvalue_l_DL <- pnorm(stat_DL, lower.tail = TRUE)

        ## DL with t Knapp and Hartnung
        pvalue_g_KH <- pt(stat_DL, -diff(dim(X)), lower.tail = FALSE)
        pvalue_d_KH <- 2 * pt(-abs(stat_DL), -diff(dim(X)))
        pvalue_l_KH <- pt(stat_DL, -diff(dim(X)), lower.tail = TRUE)


        statistics <- c(stat_LR,
                        stat_PLR,
                        stat_skovgaard,
                        stat_bartlett,
                        stat_db,
                        stat_waldML,
                        stat_DL,
                        stat_DL)

        pvalues_d <- c(pvalue_d_LR,
                       pvalue_d_PLR,
                       pvalue_d_skovgaard,
                       pvalue_d_bartlett,
                       pvalue_d_db,
                       pvalue_d_waldML,
                       pvalue_d_DL,
                       pvalue_d_KH)


        pvalues_g <- c(pvalue_g_LR,
                       pvalue_g_PLR,
                       pvalue_g_skovgaard,
                       pvalue_g_bartlett,
                       pvalue_g_db,
                       pvalue_g_waldML,
                       pvalue_g_DL,
                       pvalue_g_KH)

        pvalues_l <- c(pvalue_l_LR,
                       pvalue_l_PLR,
                       pvalue_l_skovgaard,
                       pvalue_l_bartlett,
                       pvalue_l_db,
                       pvalue_l_waldML,
                       pvalue_l_DL,
                       pvalue_l_KH)

        res <- cbind(statistics, pvalues_d, pvalues_g, pvalues_l)
    }
    dimnames(res) <- list(c("LR",
                            "PLR",
                            "Skovgaard",
                            "Bartlett",
                            "ZL",
                            "WaldML",
                            "DL",
                            "KH"),
                          c("statistics", "pvalues_d", "pvalues_g", "pvalues_l"))

    res
}


## Helper functions for BiasFit and lrtest
fitfun <- function(pars,
                   X,
                   sigma2,
                   y,
                   info = FALSE) {
  p <- ncol(X)
  zeta <- pars[p + 1]
  tau2 <- zeta
  betas <- pars[seq.int(p)]
  w <- 1/(sigma2 + tau2)
  etas <- X%*%betas
  R <- drop(y - etas)
  if (info) {
    W.X <- sqrt(w)*X
    XWX <- crossprod(W.X)
    XWXinv <- try(solve(XWX), silent = TRUE)
    if (inherits(XWXinv, "try-error"))
        XWXinv <- matrix(NA, p, p)
    h <- rowSums((X %*% XWXinv) * X) * w
    expinfo <- obsinfo <- vcove <- matrix(0, p + 1, p + 1)
    expinfo[seq.int(p), seq.int(p)] <- obsinfo[seq.int(p), seq.int(p)] <- XWX
    expinfo[p + 1, p + 1] <- sum(w^2)/2
    vcove[seq.int(p), seq.int(p)] <- XWXinv
    vcove[p + 1, p + 1] <- 2/sum(w^2)
    obsinfo[p + 1, p + 1] <- -sum(w^2)/2 + sum(w^3*R^2)
    obsinfo[seq.int(p), p + 1] <- colSums(w^2*X*R)
    obsinfo[p + 1, seq.int(p)] <- colSums(w^2*X*R)
  }
  else {
    XWX <- h <- expinfo <- obsinfo <- vcove <- NULL
  }
  zstate <- betas/sqrt(diag(vcove)[1:p])
  res <- list(sigma2 = sigma2,
              y = y,
              X = X,
              tau2 = tau2,
              betas = betas,
              R = R,
              XWX = XWX,
              expinfo = expinfo,
              vcove = vcove,
              w = w,
              h = h,
              zstate = zstate)
  c(res)
}

loglik <- function(pars,
                   X,
                   sigma2,
                   y,
                   fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = FALSE)
  }
  with(fit, {
    if (tau2 < 0) {
      return(NA)
    }
    else {
      0.5*sum(log(w)) - 0.5*sum(R^2*w)
    }
  })
}

penloglik <- function(pars,
                      X,
                      sigma2,
                      y,
                      fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = TRUE)
  }
  with(fit, {
    if (tau2 < 0) {
      return(NA)
    }
    else {
        0.5*sum(log(w)) - 0.5*sum(R^2*w) - 0.5*log(det(XWX))
    }
  })
}

scores <- function(pars,
                   X,
                   sigma2,
                   y,
                   fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = TRUE)
  }
  with(fit, {
    c(colSums(w*X*R), sum(R^2*w^2 - w)/2)
  })
}

penscores <- function(pars,
                      X,
                      sigma2,
                      y,
                      fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = TRUE)
  }
  with(fit, {
    c(colSums(w*X*R), sum(R^2*w^2 - (1 - h)*w)/2)
  })
}

bias <- function(pars, X, sigma2, y, fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = TRUE)
  }
  with(fit, {
           sum(w*h)/sum(w^2)
  })
}
