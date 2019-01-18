######################  Meta-analysis simulation study with BR #######################
## Specify path; make sure that path has a directory named results
path <- "."

library("metaLik")
library("metatest")
library("pracma")
library("boot")

library("doMC")
library("dplyr")
library("plyr")

## Expected information
expinfo <- function(pars, X, sigma2) {
  p <- ncol(X)
  w <- 1/(sigma2 + pars[p + 1])
  W.X <- sqrt(w)*X
  XWX <- crossprod(W.X)
  expinfo <- matrix(0, p + 1, p + 1)	
  expinfo[seq.int(p), seq.int(p)] <- XWX
  expinfo[p + 1, p + 1] <- sum(w^2)/2
  expinfo
}

## Std errors of betahat
sterrors <- function(pars, X, sigma2) {
  p <- ncol(X)
  w <- 1/(sigma2 + pars[p + 1])
  W.X <- sqrt(w)*X
  XWX <- crossprod(W.X)
  sterrors <- try(sqrt(diag(solve(XWX))), silent = TRUE)
  if (inherits(sterrors, "try-error"))	sterrors <- rep(NA, p)	 
  sterrors
}

## Log-likelihood
loglik <- function(pars, X, sigma2, y, fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = FALSE)
  }
  with(fit, {if (psi < 0) {return(NA)} else {0.5*sum(log(w)) - 
      0.5*sum(R^2*w)}})
}

## Penalized log-likelihood
penloglik <- function(pars, X, sigma2, y, fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = TRUE)
  }
  with(fit, 0.5*sum(log(w)) - 0.5*sum(R^2*w) - 0.5*log(det(XWX)))
  with(fit, {if (psi < 0) {return(NA)} else {0.5*sum(log(w)) - 
      0.5*sum(R^2*w) - 0.5*log(det(XWX))}})
}

## Score function for beta and psi
scores <- function(pars, X, sigma2, y, fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = TRUE)
  }
  with(fit, {
    c(colSums(w*X*R), sum(R^2*w^2 - w)/2)
  })
}

penscores <- function(pars, X, sigma2, y, fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = TRUE)
  }
  with(fit, {
    c(colSums(w*X*R), sum(R^2*w^2 - (1 - h)*w)/2)
  })
}

## Score function for psi
scores_psi <- function(psi, betas, X, sigma2, y) {
  w <- 1/(sigma2 + psi)
  R <- drop(y - X%*%betas)
  sum(R^2*w^2 - w)/2
}

## Penalized score functions for psi
penscores_psi <- function(psi, betas, X, sigma2, y, fit = NULL) {
  w <- 1/(sigma2 + psi)
  R <- drop(y - X%*%betas)
  W.X <- sqrt(w)*X
  XWX <- crossprod(W.X)
  XWXinv <- try(solve(XWX), silent = TRUE)
  if (inherits(XWXinv, "try-error"))	XWXinv <- matrix(NA, ncol(X), ncol(X))
  h <- rowSums((X %*% XWXinv) * X) * w
  sum(R^2*w^2 - (1 - h)*w)/2
}

## Optimization
estimation_algorithm <- function(startval, X, y, sigma2, fn, gr){
  ss.start <- Sys.time()
  m <- 0
  p <- ncol(X)
  prev_theta <- rep(.Machine$integer.max, p + 1)
  temp_theta <- startval
  while (any(abs(gr(prev_theta, X, sigma2, y)) > 1e-6)){
    prev_theta <- temp_theta
    w <- 1/(sigma2 + temp_theta[p + 1])
    W.X <- sqrt(w)*X
    XWX <- crossprod(W.X)
    XWy <- t(X) %*% diag(w) %*% y
    XWXinv <- solve(XWX)
    temp_theta[1:p] <- as.numeric(XWXinv %*% XWy)
    psival <- uniroot(fn, c(0, 5*startval[p + 1]), betas = temp_theta[1:p],
                      X = X, sigma2 = sigma2, y = y, check.conv = TRUE, 
                      extendInt = "upX")$root
    temp_theta[p + 1] <- psival
    m <- m + 1
    if (m > 100)  break
  } 
  ss.finish <- Sys.time()
  secs <- ss.finish - ss.start
  attributes(secs) <- NULL
  list(par = temp_theta, iter = m, secs = secs)
}

## Computes quantities for location adjustment
fitfun <- function(pars, X, sigma2, y, info = TRUE) {
  p <- ncol(X)
  psi <- pars[p + 1]
  betas <- pars[seq.int(p)]
  w <- 1/(sigma2 + psi)
  etas <- X%*%betas
  R <- drop(y - etas)
  if (info) {
    W.X <- sqrt(w)*X
    XWX <- crossprod(W.X)
    XWXinv <- try(solve(XWX), silent = TRUE)
    if (inherits(XWXinv, "try-error"))	XWXinv <- matrix(NA, p, p)
    h <- rowSums((X %*% XWXinv) * X) * w
    expinfo <- obsinfo <- vcove <- matrix(0, p + 1, p + 1)
    expinfo[seq.int(p), seq.int(p)] <- obsinfo[seq.int(p), seq.int(p)] <- XWX
    expinfo[p + 1, p + 1] <- sum(w^2)/2
    vcove[seq.int(p), seq.int(p)] <- XWXinv
    vcove[p + 1, p + 1] <- 2/sum(w^2)
    obsinfo[p + 1, p + 1] <- -sum(w^2)/2 + sum(w^3 * R^2)
    obsinfo[seq.int(p), p + 1] <- colSums(w^2 * X * R)
    obsinfo[p + 1, seq.int(p)] <- colSums(w^2 * X * R)
  }
  else {
    XWX <- h <- expinfo <- obsinfo <- vcove <- NULL
  }
  bias_psi <- - sum(w*h)/sum(w^2)
  
  res <- list(sigma2 = sigma2,
              y = y,
              X = X,
              p = p,
              K = nrow(X),
              psi = psi,
              betas = betas,
              coefs = c(betas, psi),
              betahat = betas,
              psihat = psi,
              R = R,
              XWX = XWX,
              expinfo = expinfo,
              vcovMat = vcove,
              W = diag(w),
              w = w,
              h = h,
              se = sqrt(diag(vcove)),
              bias = bias_psi)
  c(res)
}

## Finds ML and/or BR estimates and SEs
BiasFit <- function(object, method = c("all", "ML", "BR")) {
  method <- match.arg(method)
  thetahat <- object$mle
  p <- ncol(object$X)
  psihat <- thetahat[p + 1]
  zeta <- psihat
  pars <- c(object$DL, object$tau2.DL)
  ML <- BR <- rep(NA, p + 1)
  iterML <- iterBR <- timeML <- timeBR <- NA
  
  if (method == "ML" | method == "all") {
    fML <- try(estimation_algorithm(pars, X = object$X, y = object$y, 
                                    sigma2 = object$sigma2, fn = scores_psi, 
                                    gr = scores), silent = TRUE)
    if (inherits(fML, "try-error")) {
      ss.start <- Sys.time()
      oML <- optim(par = pars,fn = loglik,gr = scores,X = object$X, 
                   sigma2 = object$sigma2, y = object$y,
                   control = list(fnscale = -1), method = "BFGS")
      ML <- oML$par
      ss.finish <- Sys.time()
      secs <- ss.finish - ss.start
      attributes(secs) <- NULL
      seML <- sterrors(ML, object$X, object$sigma2)
      iterML <- unname(oML$counts[2])
      timeML <-  secs
    } else {
      ML <- fML$par
      #seML <- sterrors(ML, object$X, object$sigma2)
      iterML <- fML$iter
      timeML <- fML$secs[[1]]
    }
  }
  if (method == "BR" | method == "all") {
    fBR <- try(estimation_algorithm(pars, X = object$X, y = object$y, 
                                    sigma2 = object$sigma2, fn = penscores_psi, 
                                    gr = penscores), 
               silent = TRUE)
    if (inherits(fBR, "try-error")) {
      ss.start <- Sys.time()
      oBR <- optim(par = pars, fn = penloglik, gr = penscores, X = object$X, 
                   sigma2 = object$sigma2, y = object$y,
                   control = list(fnscale = -1), method = "BFGS")
      BR <- oBR$par
      ss.finish <- Sys.time()
      secs <- ss.finish - ss.start
      attributes(secs) <- NULL
      seBR <- sterrors(BR, object$X, object$sigma2)
      iterBR <- unname(oBR$counts[2])
      timeBR <-  secs
    } else {
      BR <- fBR$par
      iterBR <- fBR$iter
      timeBR <- fBR$secs[[1]]
    }
  }
  list(ML = ML, BR = BR, iters = c(iterML, iterBR), 
       secs = c(timeML, timeBR))
}

### Implementation of t*

## The derivatives of the standard errors wrt to the parameter arg
dse <- function(pars, X, sigma2, y, which = 1, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  with(fit, {
    Wstar <- diag(-(sigma2 + psihat)^(-2))
    vstar <- - sum((sigma2 + psihat)^(-3))
    
    dinfoB <- lapply(seq.int(p), function(var) matrix(0, nrow = p + 1, 
                       ncol = p + 1))
    dinfoP <- lapply(seq.int(1), function(var) cbind(rbind(t(X) %*% 
                       Wstar%*%X, rep(0, p)), c(rep(0, p), vstar)))
    
    dinfo <- c(dinfoB, dinfoP)
    
    sapply(seq.int(p + 1), function(var) {
      - 0.5*diag(vcovMat %*% dinfo[[var]] %*% vcovMat)[which]/se[which]
    })
  })
}

## The second derivatives of the standard errors wrt to the parameter arg
ddse <- function(pars, X, sigma2, y, which = 1, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  with(fit, {
    Wstar <- diag(-(sigma2 + psihat)^(-2))
    Wstarstar <- diag(2*(sigma2 + psihat)^(-3))
    vstar <- - sum((sigma2 + psihat)^(-3))
    vstarstar <- 3*sum((sigma2 + psihat)^(-4))
    
    dinfoB <- lapply(seq.int(p), function(var) matrix(0, nrow = p + 1, 
                                                      ncol = p + 1))
    dinfoP <- lapply(seq.int(1), function(var) cbind(rbind(t(X) %*% 
                         Wstar %*% X, rep(0, p)), c(rep(0, p), vstar)))
    
    dinfo <- function(s) {
      c(dinfoB, dinfoP)[[s]]
    }
    
    ddinfoPP <- cbind(rbind(t(X) %*% Wstarstar %*% X, rep(0, p)), 
                      c(rep(0, p), vstarstar))
    
    ddinfo <- function(s,t) {
      if((is.element(s, 1:p)) & (is.element(t, 1:p))) ddinfo <- 
          matrix(0, p + 1, p + 1)
      if((is.element(s, 1:p)) & (t == p + 1)) ddinfo <- 
          matrix(0, p + 1, p + 1)
      if((s == p + 1) & (is.element(t, 1:p))) ddinfo <- 
          matrix(0, p + 1, p + 1)
      if((s == p + 1) & (t == p + 1)) ddinfo <- ddinfoPP
      return(ddinfo)
    }
    
    der <- function(s, t) {
      - 0.25 * diag(vcovMat %*% dinfo(s) %*% vcovMat)[which] *
        diag(vcovMat %*% dinfo(t) %*% vcovMat)[which]/se[which]^3 +
        0.5 * diag(vcovMat %*% dinfo(s) %*% vcovMat %*% dinfo(t) 
             %*% vcovMat)[which]/se[which] + 0.5 * diag(vcovMat 
                 %*% dinfo(t) %*% vcovMat %*% dinfo(s) %*% 
                     vcovMat)[which]/se[which] - 0.5 * diag(vcovMat 
                        %*% ddinfo(s, t) %*% vcovMat)[which]/se[which]
    }
    outer(seq.int(p + 1), seq.int(p + 1), Vectorize(der))
  })
}

tau <- function(pars, X, sigma2, y, null = 0, which = 1, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  with(fit, {
    (coefs[which] - null) * se[which]^(-1)
  })
}

dtau <- function(pars, X, sigma2, y, null = 0, which = 1, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  tau.val <- tau(null = null, which = which, fit = fit)
  dse.val <- dse(which = which, fit = fit)
  with(fit, {
    evec <- rep(0, p + 1)
    evec[which] <- 1
    (evec - tau.val * dse.val)/se[which]
  })
}

ddtau <- function(pars, X, sigma2, y, null = 0, which = 1, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  tau.val <- tau(null = null, which = which, fit = fit)
  dtau.val <- dtau(null = null, which = which, fit = fit)
  dse.val <- dse(which = which, fit = fit)
  ddse.val <- ddse(which = which, fit = fit)
  with(fit, {
    -(dtau.val %*% t(dse.val) + dse.val %*% t(dtau.val) +
        tau.val * ddse.val)/se[which]
  })
}

waldi_meta <- function(pars, X, sigma2, y, null = 0, which = 1, 
                       method = c("ML", "BR"), adjust = TRUE, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  tau <- tau(null = null, which = which, fit = fit)
  nam <- names(tau)
  if (!adjust) {
    tau
  }
  else {
    with(fit, {
      traces <- sum(diag(vcovMat %*% ddtau(null = null, which = which, 
                                           fit = fit)))
      if (method == "BR") tauc <- tau - 0.5 * traces
      else {
        bias_psi <- bias
        biasvec <- c(rep(0, p), bias_psi)
        dtau <-  dtau(null = null, which = which, fit = fit)
        tauc <- tau - (dtau %*% biasvec + 0.5 * traces)
      }
      
      tauc <- as.double(tauc)
      names(tauc) <- nam
      tauc
    })
  }
}

## Implementation Zeng and Lin's double resampling
double.resampling <- function(beta0, m, B = 1000, myseed = 123, 
                              ci = FALSE, alpha = 0.05,
                              alternative = c("two.sided", "less", 
                                              "greater", "all")) {
  alternative <- match.arg(alternative) ## Currently inactive
  ## take care of the random seed
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = 
             FALSE)) {
    seed.keep <- get(".Random.seed", envir = .GlobalEnv,
                     inherits = FALSE)
    on.exit(assign(".Random.seed", seed.keep, envir = 
                     .GlobalEnv))
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
  denom <- S1 - S2/S1
  compute.tau2 <- function(betak, betahat) {
    eps <- t(betak) - c(betahat)
    tmp <- eps^2 %*% (1/sigma2 ) - (K - 1)
    pmax(0, tmp/denom)
  }
  compute.betahat <- function(betak, tau2) {
    t(betak) %*% (1/(sigma2 + tau2) )/sum(1/(sigma2 + 
                                               tau2) )
  }
  ## first run
  betak <- replicate(B, beta.DL + rnorm(K) * sqrt(sigma2 + 
                                                    tau2.DL) )
  betahat <- compute.betahat(betak, 0.0)
  tau2hat <- compute.tau2(betak, betahat)
  ## second run
  for (b in 1:B) {
    betak <- replicate(B, beta.DL + rnorm(K) * sqrt(sigma2 + 
                                                      tau2hat[b]))
    beta.resampled[b, ] <- compute.betahat(betak, tau2hat[b])
  }
  beta.resampled <- c(beta.resampled)
  pval <- switch(alternative,
                 "two.sided" = 2 * min(mean(beta.resampled > 
                                              beta0), 
                                       mean(beta.resampled < 
                                              beta0)),
                 "less" = mean(beta.resampled < beta0),
                 "greater" = mean(beta.resampled > beta0),
                 "all" = c(two.sided = 2 * min (mean(beta.resampled > 
                                                       beta0), 
                                                mean(beta.resampled < 
                                                       beta0)),
                           less = mean(beta.resampled < beta0),
                           greater = mean(beta.resampled > beta0)))
  if(ci) {
    ci <- switch(alternative,
                 "two.sided" = quantile(beta.resampled, c(alpha/2, 1 - 
                                                            alpha/2)),
                 "less" = c(-Inf, quantile(beta.resampled, 1 - alpha)),
                 "greater" = c(quantile(beta.resampled, alpha), Inf),
                 "all" = list(two.sided = quantile(beta.resampled, 
                                                   c(alpha/2, 1 - alpha/2)),
                              less = c(-Inf, quantile(beta.resampled, 
                                                      1 - alpha)),
                              greater = c(quantile(beta.resampled, 
                                                   alpha), Inf)))
    return(list(pval = pval, ci = ci))
  }
  else
    return(pval)
}


## Calculates p-values and statistics for a single parameter
perform_tests <- function(y, X, sigma2, what = 1, null = 0, Kval,
                          B = 1000, silent = TRUE, maxiter = 1000, 
                          seed = 123) {
  
  fit_metaLik <- try(metaLik(y ~ -1 + X, sigma2 = sigma2), 
                     silent = silent)
  
  if (inherits(fit_metaLik, "try-error")) res <- matrix(NA, nrow = 3, ncol = 4)
  else {
    ## Double resampling - Zeng & Lin
    stat_db <- NA
    if (ncol(X) == 1 & all(X == 1) & (K < 100)) {
      pvalue_db <- double.resampling(beta0 = null, m = fit_metaLik, B = B, 
                          alternative = "all", ci = FALSE, myseed = seed)
      pvalue_d_db <- pvalue_db["two.sided"]
      pvalue_l_db <- pvalue_db["less"]
      pvalue_g_db <- pvalue_db["greater"]
    }
    else {
      pvalue_d_db <- pvalue_l_db <- pvalue_g_db <- NA
    }
    ## Wald based on DL
    stat_DL <- (fit_metaLik$DL[what] - null)/sqrt(fit_metaLik$vcov.DL[what, what])
    pvalue_g_DL <- pnorm(stat_DL, lower.tail = FALSE)
    pvalue_d_DL <- 2 * pnorm(- abs(stat_DL))
    pvalue_l_DL <- pnorm(stat_DL, lower.tail = TRUE)
    ## DL with t Knapp and Hartnung
    pvalue_g_KH <- pt(stat_DL, - diff(dim(X)), lower.tail = FALSE)
    pvalue_d_KH <- 2 * pt(- abs(stat_DL), - diff(dim(X)))
    pvalue_l_KH <- pt(stat_DL, - diff(dim(X)), lower.tail = TRUE)
    statistics <- c(stat_db, stat_DL, stat_DL)
    pvalues_d <- c(pvalue_d_db, pvalue_d_DL, pvalue_d_KH)
    pvalues_g <- c(pvalue_g_db, pvalue_g_DL, pvalue_g_KH)
    pvalues_l <- c(pvalue_l_db, pvalue_l_DL, pvalue_l_KH)
    res <- cbind(statistics, pvalues_d, pvalues_g, pvalues_l)
  }
  dimnames(res) <- list(c("ZL", "DL", "KH"), c("statistics", "pval_2sided", 
                                               "pval_right", "pval_left"))
  res
}

## as in Zeng and Lin (2015)
generate.sigma2s <- function(K){
  sigma2s <- 0.25 * rchisq(K, df = 1)
  ok <- (sigma2s > 0.009 & sigma2s < 0.6)
  while (sum(ok) < K) {
    tmp <- 0.25 * rchisq(K - sum(ok), df = 1)
    sigma2s[!ok] <- tmp
    ok <- (sigma2s > 0.009 & sigma2s < 0.6)
  }
  sigma2s
}

simulate.BG <- function(beta = 0.5, var.re = 0.01,
                        sigma2s = runif(10, 0, 1)) {
  N <- length(sigma2s)
  list(K = N, truebeta = beta, truepsi = var.re, 
       betahat = beta + rnorm(N)*sqrt(var.re + sigma2s), 
       betahat.se = sqrt(sigma2s), X = matrix(1, nrow = N))
}

## sample sizes to consider
Ks <- c(seq(5, 50, by = 5), 100, 200)

## variance components to consider
truepsis <- seq(0, 0.1, length = 11)
## other simulation constants
truebeta <- 0.5

## compute X matrices and generate sigma2
sigma2 <- Xmat <- as.list(rep(0, length(Ks)))
## generate sigmas and set Xs
set.seed(123)
for (j in seq_along(Ks)) {
  sigma2[[j]] <- generate.sigma2s(Ks[j])
  Xmat[[j]] <- matrix(1, nrow = Ks[j])
}

## Specify number of cores, simulation size, bootstrap replications
registerDoMC(40)
nsimu <- 10000
nboot_ZL <- 1000

res <- list()

### fixed K

for(K in c(10, 20)) {
  cat(paste("K = ", K), "\n")
  
  for(psi in truepsis) {
    cat(paste("psi = ", psi), "\n")
    
    set.seed(123)
    simu_data_mat <- matrix(NA, nrow = K, ncol = nsimu)
    
    for (i in seq.int(nsimu)) {
      dat <- simulate.BG(beta = truebeta, var.re = psi, sigma2s = 
                           sigma2[[which(Ks == K)]])
      simu_data_mat[, i] <- dat$betahat
    }
    
    simu_data <- data.frame(t(simu_data_mat))
    simu_data$sample <- seq.int(nsimu)
    
    ## Bootstrap seeds
    simu_data$seeds <- sample(seq.int(nsimu * 100), nsimu, 
                              replace = FALSE)
    
    attr(simu_data, "beta") <- truebeta
    attr(simu_data, "K") <- K
    attr(simu_data, "psi") <- psi
    attr(simu_data, "sigma2") <- sigma2[[which(Ks == K)]]
    attr(simu_data, "Nsim") <- nsimu
    
    temp_res <- dlply(simu_data, ~ sample, function(response) {
      temp_data <- unname(unlist(response[-(K + 1:2)]))
      
      ## Preliminary metaLik fit
      temp_fit <- metaLik(temp_data ~ 1, sigma2 = sigma2[[which(Ks == K)]])
      ## more accurate fit via Biasfit
      bias_fit <- BiasFit(temp_fit, method = "all")
      
      ## quantities for location adjustment
      # ML
      fitML <- fitfun(bias_fit$ML, temp_fit$X, temp_fit$sigma2, temp_fit$y)
      # BR
      fitBR <- fitfun(bias_fit$BR, temp_fit$X, temp_fit$sigma2, temp_fit$y)
      
      # Wald statistics
      zstat_ml <- waldi_meta(null = truebeta, which = 1, adjust = FALSE, 
                             method = "ML", fit = fitML)
      zstat_br <- waldi_meta(null = truebeta, which = 1, adjust = FALSE, 
                             method = "BR", fit = fitBR)
      zstat_ml_cor <- waldi_meta(null = truebeta, which = 1, adjust = TRUE, 
                                 method = "ML", fit = fitML)
      zstat_br_cor <- waldi_meta(null = truebeta, which = 1, adjust = TRUE, 
                                 method = "BR", fit = fitBR)
      # ZL and DL statistics
      stat_values <- perform_tests(temp_fit$y, temp_fit$X, Kval = K,
                                   temp_fit$sigma2, null = truebeta, 
                                   what = 1, B = nboot_ZL, silent = FALSE, 
                                   maxiter = 10000, seed = response$seeds)
      
      if (response$sample %% 100 == 0) cat(response$sample, "\n")
      stats <- data.frame(name = c("ml", "ml_cor", "br", "br_cor", "DL", "KH"),
                          value = c(zstat_ml, zstat_ml_cor, zstat_br, zstat_br_cor,
                                    stat_values["DL", 1], stat_values["KH", 1]))
      pvalues <- data.frame(statistic = c(rep(c("ZL", "KH"), times = 3)),
                            pvalue = c(stat_values["ZL", 2], stat_values["KH", 2],
                                       stat_values["ZL", 3], stat_values["KH", 3],
                                       stat_values["ZL", 4], stat_values["KH", 4]),
                            type = rep(c("pval_2sided", "pval_right", "pval_left"), 
                                         each = 2))
      list(id_sample = response$sample, sample = temp_data, stats = stats, 
           pvalues = pvalues, K = attr(simu_data, "K"), 
           psi = attr(simu_data, "psi"), beta = attr(simu_data, "beta"), 
           sigma2 = attr(simu_data, "sigma2"))                                       
    }, .parallel = TRUE)
    
    res <- c(res, temp_res)
    save.image(paste0("~/metaK=", K, "psi=", psi, ".rda"))
    
  }   
}

## Summary of results
res_stats <- ldply(res, function(x) data.frame(x$stats, K = x$K, psi = x$psi))
res_pvals <- ldply(res, function(x) data.frame(x$pvalues, K = x$K, psi = x$psi))

typeI_statistics <- ddply(res_stats, ~ name + K + psi, function(x) {
  levels <- c(0.1, 1, 2.5, 5)/100
  p_value_2sided <- 2 * pnorm(-abs(x$value))
  p_value_left <- pnorm(x$value)
  p_value_right <- 1 - pnorm(x$value)
  rate_2sided <- sapply(levels, function(alpha) mean(p_value_2sided < alpha))
  rate_left <- sapply(levels, function(alpha) mean(p_value_left < alpha))
  rate_right <- sapply(levels, function(alpha) mean(p_value_right < alpha))
  out <- data.frame(
    test = rep(c("2sided", "left", "right"), each = length(levels)),
    typeI = c(rate_2sided, rate_left, rate_right),
    level = rep(levels, times = 3))
  out
})

typeI_pvalues <-
  ddply(res_pvals, ~ statistic + K + psi, function(x) {
    levels <- c(0.1, 1, 2.5, 5)/100
    rate_2sided <- sapply(levels, function(alpha) 
      mean(x$pvalue[x$type == 'pval_2sided'] < alpha))
    rate_left <- sapply(levels, function(alpha) 
      mean(x$pvalue[x$type == 'pval_left'] < alpha))
    rate_right <- sapply(levels, function(alpha) 
      mean(x$pvalue[x$type == 'pval_right'] < alpha))
    out <- data.frame(
      test = rep(c("2sided", "left", "right"), each = length(levels)),
      typeI = c(rate_2sided, rate_left, rate_right),
      level = rep(levels, times = 3))
    out
  })


names(typeI_statistics) <- names(typeI_pvalues)

typeIdf <- rbind(typeI_statistics, typeI_pvalues)

typeIdf <- typeIdf %>% 
  filter(is.element(statistic, c("ml", "ml_cor", "br", "br_cor",
                                 "DL", "ZL")), test == "2sided", level == 0.05) %>%
  mutate(level_chr = paste(level * 100, "~symbol('\045')"),
         cov = (1 - typeI) * 100,
         cov_chr = paste((1 - level) * 100, "~symbol('\045')"),
         upper = typeI - qnorm(1 - 0.01/2) * sqrt(typeI * (1 - typeI)/nsimu),
         lower = typeI + qnorm(1 - 0.01/2) * sqrt(typeI * (1 - typeI)/nsimu))

typeIdf$statistic <- factor(typeIdf$statistic, levels = c("ml", "ml_cor", 
                                                          "br", "br_cor", "DL", "ZL"), ordered = TRUE)
typeIdf$Klab <- paste0("K = ", typeIdf$K)

typeIdf$Klab <- factor(typeIdf$Klab, levels = c("K = 10", "K = 20"), labels = 
                         c("italic(K) == 10", "italic(K) == 20"), ordered = TRUE)

cov_df_K <- typeIdf

## fixed psi

## Expected information
expinfo <- function(pars, X, sigma2) {
  p <- ncol(X)
  w <- 1/(sigma2 + pars[p + 1])
  W.X <- sqrt(w)*X
  XWX <- crossprod(W.X)
  expinfo <- matrix(0, p + 1, p + 1)	
  expinfo[seq.int(p), seq.int(p)] <- XWX
  expinfo[p + 1, p + 1] <- sum(w^2)/2
  expinfo
}

## Std errors of betahat
sterrors <- function(pars, X, sigma2) {
  p <- ncol(X)
  w <- 1/(sigma2 + pars[p + 1])
  W.X <- sqrt(w)*X
  XWX <- crossprod(W.X)
  sterrors <- try(sqrt(diag(solve(XWX))), silent = TRUE)
  if (inherits(sterrors, "try-error"))	sterrors <- rep(NA, p)	 
  sterrors
}

## Log-likelihood
loglik <- function(pars, X, sigma2, y, fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = FALSE)
  }
  with(fit, {if (psi < 0) {return(NA)} else {0.5*sum(log(w)) - 
      0.5*sum(R^2*w)}})
}

## Penalized log-likelihood
penloglik <- function(pars, X, sigma2, y, fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = TRUE)
  }
  with(fit, 0.5*sum(log(w)) - 0.5*sum(R^2*w) - 0.5*log(det(XWX)))
  with(fit, {if (psi < 0) {return(NA)} else {0.5*sum(log(w)) - 
      0.5*sum(R^2*w) - 0.5*log(det(XWX))}})
}

## Score function for beta and psi
scores <- function(pars, X, sigma2, y, fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = TRUE)
  }
  with(fit, {
    c(colSums(w*X*R), sum(R^2*w^2 - w)/2)
  })
}

penscores <- function(pars, X, sigma2, y, fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = TRUE)
  }
  with(fit, {
    c(colSums(w*X*R), sum(R^2*w^2 - (1 - h)*w)/2)
  })
}

## Score function for psi
scores_psi <- function(psi, betas, X, sigma2, y) {
  w <- 1/(sigma2 + psi)
  R <- drop(y - X%*%betas)
  sum(R^2*w^2 - w)/2
}

## Penalized score functions for psi
penscores_psi <- function(psi, betas, X, sigma2, y, fit = NULL) {
  w <- 1/(sigma2 + psi)
  R <- drop(y - X%*%betas)
  W.X <- sqrt(w)*X
  XWX <- crossprod(W.X)
  XWXinv <- try(solve(XWX), silent = TRUE)
  if (inherits(XWXinv, "try-error"))	XWXinv <- matrix(NA, ncol(X), ncol(X))
  h <- rowSums((X %*% XWXinv) * X) * w
  sum(R^2*w^2 - (1 - h)*w)/2
}

## Optimization
estimation_algorithm <- function(startval, X, y, sigma2, fn, gr){
  ss.start <- Sys.time()
  m <- 0
  p <- ncol(X)
  prev_theta <- rep(.Machine$integer.max, p + 1)
  temp_theta <- startval
  while (any(abs(gr(prev_theta, X, sigma2, y)) > 1e-6)){
    prev_theta <- temp_theta
    w <- 1/(sigma2 + temp_theta[p + 1])
    W.X <- sqrt(w)*X
    XWX <- crossprod(W.X)
    XWy <- t(X) %*% diag(w) %*% y
    XWXinv <- solve(XWX)
    temp_theta[1:p] <- as.numeric(XWXinv %*% XWy)
    psival <- uniroot(fn, c(0, 5*startval[p + 1]), betas = temp_theta[1:p],
                      X = X, sigma2 = sigma2, y = y, check.conv = TRUE, 
                      extendInt = "upX")$root
    temp_theta[p + 1] <- psival
    m <- m + 1
    if (m > 100)  break
  } 
  ss.finish <- Sys.time()
  secs <- ss.finish - ss.start
  attributes(secs) <- NULL
  list(par = temp_theta, iter = m, secs = secs)
}

## Computes quantities for location adjustment
fitfun <- function(pars, X, sigma2, y, info = TRUE) {
  p <- ncol(X)
  psi <- pars[p + 1]
  betas <- pars[seq.int(p)]
  w <- 1/(sigma2 + psi)
  etas <- X%*%betas
  R <- drop(y - etas)
  if (info) {
    W.X <- sqrt(w)*X
    XWX <- crossprod(W.X)
    XWXinv <- try(solve(XWX), silent = TRUE)
    if (inherits(XWXinv, "try-error"))	XWXinv <- matrix(NA, p, p)
    h <- rowSums((X %*% XWXinv) * X) * w
    expinfo <- obsinfo <- vcove <- matrix(0, p + 1, p + 1)
    expinfo[seq.int(p), seq.int(p)] <- obsinfo[seq.int(p), seq.int(p)] <- XWX
    expinfo[p + 1, p + 1] <- sum(w^2)/2
    vcove[seq.int(p), seq.int(p)] <- XWXinv
    vcove[p + 1, p + 1] <- 2/sum(w^2)
    obsinfo[p + 1, p + 1] <- -sum(w^2)/2 + sum(w^3 * R^2)
    obsinfo[seq.int(p), p + 1] <- colSums(w^2 * X * R)
    obsinfo[p + 1, seq.int(p)] <- colSums(w^2 * X * R)
  }
  else {
    XWX <- h <- expinfo <- obsinfo <- vcove <- NULL
  }
  bias_psi <- - sum(w*h)/sum(w^2)
  
  res <- list(sigma2 = sigma2,
              y = y,
              X = X,
              p = p,
              K = nrow(X),
              psi = psi,
              betas = betas,
              coefs = c(betas, psi),
              betahat = betas,
              psihat = psi,
              R = R,
              XWX = XWX,
              expinfo = expinfo,
              vcovMat = vcove,
              W = diag(w),
              w = w,
              h = h,
              se = sqrt(diag(vcove)),
              bias = bias_psi)
  c(res)
}

## Finds ML and/or BR estimates and SEs
BiasFit <- function(object, method = c("all", "ML", "BR")) {
  method <- match.arg(method)
  thetahat <- object$mle
  p <- ncol(object$X)
  psihat <- thetahat[p + 1]
  zeta <- psihat
  pars <- c(object$DL, object$tau2.DL)
  ML <- BR <- rep(NA, p + 1)
  iterML <- iterBR <- timeML <- timeBR <- NA
  
  if (method == "ML" | method == "all") {
    fML <- try(estimation_algorithm(pars, X = object$X, y = object$y, 
                                    sigma2 = object$sigma2, fn = scores_psi, 
                                    gr = scores), silent = TRUE)
    if (inherits(fML, "try-error")) {
      ss.start <- Sys.time()
      oML <- optim(par = pars,fn = loglik,gr = scores,X = object$X, 
                   sigma2 = object$sigma2, y = object$y,
                   control = list(fnscale = -1), method = "BFGS")
      ML <- oML$par
      ss.finish <- Sys.time()
      secs <- ss.finish - ss.start
      attributes(secs) <- NULL
      seML <- sterrors(ML, object$X, object$sigma2)
      iterML <- unname(oML$counts[2])
      timeML <-  secs
    } else {
      ML <- fML$par
      #seML <- sterrors(ML, object$X, object$sigma2)
      iterML <- fML$iter
      timeML <- fML$secs[[1]]
    }
  }
  if (method == "BR" | method == "all") {
    fBR <- try(estimation_algorithm(pars, X = object$X, y = object$y, 
                                    sigma2 = object$sigma2, fn = penscores_psi, 
                                    gr = penscores), 
               silent = TRUE)
    if (inherits(fBR, "try-error")) {
      ss.start <- Sys.time()
      oBR <- optim(par = pars, fn = penloglik, gr = penscores, X = object$X, 
                   sigma2 = object$sigma2, y = object$y,
                   control = list(fnscale = -1), method = "BFGS")
      BR <- oBR$par
      ss.finish <- Sys.time()
      secs <- ss.finish - ss.start
      attributes(secs) <- NULL
      seBR <- sterrors(BR, object$X, object$sigma2)
      iterBR <- unname(oBR$counts[2])
      timeBR <-  secs
    } else {
      BR <- fBR$par
      iterBR <- fBR$iter
      timeBR <- fBR$secs[[1]]
    }
  }
  list(ML = ML, BR = BR, iters = c(iterML, iterBR), 
       secs = c(timeML, timeBR))
}

### Implementation of t*

## The derivatives of the standard errors wrt to the parameter arg
dse <- function(pars, X, sigma2, y, which = 1, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  with(fit, {
    Wstar <- diag(-(sigma2 + psihat)^(-2))
    vstar <- - sum((sigma2 + psihat)^(-3))
    
    dinfoB <- lapply(seq.int(p), function(var) matrix(0, nrow = p + 1, 
                                                      ncol = p + 1))
    dinfoP <- lapply(seq.int(1), function(var) cbind(rbind(t(X) %*% 
                         Wstar%*%X, rep(0, p)), c(rep(0, p), vstar)))
    
    dinfo <- c(dinfoB, dinfoP)
    
    sapply(seq.int(p + 1), function(var) {
      - 0.5*diag(vcovMat %*% dinfo[[var]] %*% vcovMat)[which]/se[which]
    })
  })
}

## The second derivatives of the standard errors wrt to the parameter arg
ddse <- function(pars, X, sigma2, y, which = 1, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  with(fit, {
    Wstar <- diag(-(sigma2 + psihat)^(-2))
    Wstarstar <- diag(2*(sigma2 + psihat)^(-3))
    vstar <- - sum((sigma2 + psihat)^(-3))
    vstarstar <- 3*sum((sigma2 + psihat)^(-4))
    
    dinfoB <- lapply(seq.int(p), function(var) matrix(0, nrow = p + 1, 
                                                      ncol = p + 1))
    dinfoP <- lapply(seq.int(1), function(var) cbind(rbind(t(X) %*% 
                            Wstar %*% X, rep(0, p)), c(rep(0, p), vstar)))
    
    dinfo <- function(s) {
      c(dinfoB, dinfoP)[[s]]
    }
    
    ddinfoPP <- cbind(rbind(t(X) %*% Wstarstar %*% X, rep(0, p)), 
                      c(rep(0, p), vstarstar))
    
    ddinfo <- function(s,t) {
      if((is.element(s, 1:p)) & (is.element(t, 1:p))) ddinfo <- 
          matrix(0, p + 1, p + 1)
      if((is.element(s, 1:p)) & (t == p + 1)) ddinfo <- 
          matrix(0, p + 1, p + 1)
      if((s == p + 1) & (is.element(t, 1:p))) ddinfo <- 
          matrix(0, p + 1, p + 1)
      if((s == p + 1) & (t == p + 1)) ddinfo <- ddinfoPP
      return(ddinfo)
    }
    
    der <- function(s, t) {
      - 0.25 * diag(vcovMat %*% dinfo(s) %*% vcovMat)[which] *
        diag(vcovMat %*% dinfo(t) %*% vcovMat)[which]/se[which]^3 +
        0.5 * diag(vcovMat %*% dinfo(s) %*% vcovMat %*% dinfo(t) 
          %*% vcovMat)[which]/se[which] + 0.5 * diag(vcovMat 
              %*% dinfo(t) %*% vcovMat %*% dinfo(s) %*% 
                  vcovMat)[which]/se[which] - 0.5 * diag(vcovMat 
                      %*% ddinfo(s, t) %*% vcovMat)[which]/se[which]
    }
    outer(seq.int(p + 1), seq.int(p + 1), Vectorize(der))
  })
}

tau <- function(pars, X, sigma2, y, null = 0, which = 1, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  with(fit, {
    (coefs[which] - null) * se[which]^(-1)
  })
}

dtau <- function(pars, X, sigma2, y, null = 0, which = 1, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  tau.val <- tau(null = null, which = which, fit = fit)
  dse.val <- dse(which = which, fit = fit)
  with(fit, {
    evec <- rep(0, p + 1)
    evec[which] <- 1
    (evec - tau.val * dse.val)/se[which]
  })
}

ddtau <- function(pars, X, sigma2, y, null = 0, which = 1, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  tau.val <- tau(null = null, which = which, fit = fit)
  dtau.val <- dtau(null = null, which = which, fit = fit)
  dse.val <- dse(which = which, fit = fit)
  ddse.val <- ddse(which = which, fit = fit)
  with(fit, {
    -(dtau.val %*% t(dse.val) + dse.val %*% t(dtau.val) +
        tau.val * ddse.val)/se[which]
  })
}

waldi_meta <- function(pars, X, sigma2, y, null = 0, which = 1, 
                       method = c("ML", "BR"), adjust = TRUE, fit = NULL) {
  if (is.null(fit)) fit <- fitfun(pars, X, sigma2, y)
  tau <- tau(null = null, which = which, fit = fit)
  nam <- names(tau)
  if (!adjust) {
    tau
  }
  else {
    with(fit, {
      traces <- sum(diag(vcovMat %*% ddtau(null = null, which = which, 
                                           fit = fit)))
      if (method == "BR") tauc <- tau - 0.5 * traces
      else {
        bias_psi <- bias
        biasvec <- c(rep(0, p), bias_psi)
        dtau <-  dtau(null = null, which = which, fit = fit)
        tauc <- tau - (dtau %*% biasvec + 0.5 * traces)
      }
      
      tauc <- as.double(tauc)
      names(tauc) <- nam
      tauc
    })
  }
}

## Implementation Zeng and Lin's double resampling
double.resampling <- function(beta0, m, B = 1000, myseed = 123, 
                              ci = FALSE, alpha = 0.05,
                              alternative = c("two.sided", "less", 
                                              "greater", "all")) {
  alternative <- match.arg(alternative) ## Currently inactive
  ## take care of the random seed
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = 
             FALSE)) {
    seed.keep <- get(".Random.seed", envir = .GlobalEnv,
                     inherits = FALSE)
    on.exit(assign(".Random.seed", seed.keep, envir = 
                     .GlobalEnv))
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
  denom <- S1 - S2/S1
  compute.tau2 <- function(betak, betahat) {
    eps <- t(betak) - c(betahat)
    tmp <- eps^2 %*% (1/sigma2 ) - (K - 1)
    pmax(0, tmp/denom)
  }
  compute.betahat <- function(betak, tau2) {
    t(betak) %*% (1/(sigma2 + tau2) )/sum(1/(sigma2 + 
                                               tau2) )
  }
  ## first run
  betak <- replicate(B, beta.DL + rnorm(K) * sqrt(sigma2 + 
                                                    tau2.DL) )
  betahat <- compute.betahat(betak, 0.0)
  tau2hat <- compute.tau2(betak, betahat)
  ## second run
  for (b in 1:B) {
    betak <- replicate(B, beta.DL + rnorm(K) * sqrt(sigma2 + 
                                                      tau2hat[b]))
    beta.resampled[b, ] <- compute.betahat(betak, tau2hat[b])
  }
  beta.resampled <- c(beta.resampled)
  pval <- switch(alternative,
                 "two.sided" = 2 * min(mean(beta.resampled > 
                                              beta0), 
                                       mean(beta.resampled < 
                                              beta0)),
                 "less" = mean(beta.resampled < beta0),
                 "greater" = mean(beta.resampled > beta0),
                 "all" = c(two.sided = 2 * min (mean(beta.resampled > 
                                                       beta0), 
                                                mean(beta.resampled < 
                                                       beta0)),
                           less = mean(beta.resampled < beta0),
                           greater = mean(beta.resampled > beta0)))
  if(ci) {
    ci <- switch(alternative,
                 "two.sided" = quantile(beta.resampled, c(alpha/2, 1 - 
                                                            alpha/2)),
                 "less" = c(-Inf, quantile(beta.resampled, 1 - alpha)),
                 "greater" = c(quantile(beta.resampled, alpha), Inf),
                 "all" = list(two.sided = quantile(beta.resampled, 
                                                   c(alpha/2, 1 - alpha/2)),
                              less = c(-Inf, quantile(beta.resampled, 
                                                      1 - alpha)),
                              greater = c(quantile(beta.resampled, 
                                                   alpha), Inf)))
    return(list(pval = pval, ci = ci))
  }
  else
    return(pval)
}


## Calculates p-values and statistics for a single parameter
perform_tests <- function(y, X, sigma2, what = 1, null = 0, Kval,
                          B = 1000, silent = TRUE, maxiter = 1000, 
                          seed = 123) {
  
  fit_metaLik <- try(metaLik(y ~ -1 + X, sigma2 = sigma2), 
                     silent = silent)
  
  if (inherits(fit_metaLik, "try-error")) res <- matrix(NA, nrow = 3, 
                                                        ncol = 4)
  else {
    ## Double resampling - Zeng & Lin
    stat_db <- NA
    if (ncol(X) == 1 & all(X == 1) & (K < 100)) {
      pvalue_db <- double.resampling(beta0 = null, m = fit_metaLik, B = B, 
                          alternative = "all", ci = FALSE, myseed = seed)
      pvalue_d_db <- pvalue_db["two.sided"]
      pvalue_l_db <- pvalue_db["less"]
      pvalue_g_db <- pvalue_db["greater"]
    }
    else {
      pvalue_d_db <- pvalue_l_db <- pvalue_g_db <- NA
    }
    ## Wald based on DL
    stat_DL <- (fit_metaLik$DL[what] - null)/sqrt(fit_metaLik$vcov.DL[what, what])
    pvalue_g_DL <- pnorm(stat_DL, lower.tail = FALSE)
    pvalue_d_DL <- 2 * pnorm(- abs(stat_DL))
    pvalue_l_DL <- pnorm(stat_DL, lower.tail = TRUE)
    ## DL with t Knapp and Hartnung
    pvalue_g_KH <- pt(stat_DL, - diff(dim(X)), lower.tail = FALSE)
    pvalue_d_KH <- 2 * pt(- abs(stat_DL), - diff(dim(X)))
    pvalue_l_KH <- pt(stat_DL, - diff(dim(X)), lower.tail = TRUE)
    statistics <- c(stat_db, stat_DL, stat_DL)
    pvalues_d <- c(pvalue_d_db, pvalue_d_DL, pvalue_d_KH)
    pvalues_g <- c(pvalue_g_db, pvalue_g_DL, pvalue_g_KH)
    pvalues_l <- c(pvalue_l_db, pvalue_l_DL, pvalue_l_KH)
    res <- cbind(statistics, pvalues_d, pvalues_g, pvalues_l)
  }
  dimnames(res) <- list(c("ZL", "DL", "KH"), c("statistics", "pval_2sided", 
                                               "pval_right", "pval_left"))
  res
}

## as in Zeng and Lin (2015)
generate.sigma2s <- function(K){
  sigma2s <- 0.25 * rchisq(K, df = 1)
  ok <- (sigma2s > 0.009 & sigma2s < 0.6)
  while (sum(ok) < K) {
    tmp <- 0.25 * rchisq(K - sum(ok), df = 1)
    sigma2s[!ok] <- tmp
    ok <- (sigma2s > 0.009 & sigma2s < 0.6)
  }
  sigma2s
}

simulate.BG <- function(beta = 0.5, var.re = 0.01,
                        sigma2s = runif(10, 0, 1)) {
  N <- length(sigma2s)
  list(K = N, truebeta = beta, truepsi = var.re, 
       betahat = beta + rnorm(N)*sqrt(var.re + sigma2s), 
       betahat.se = sqrt(sigma2s), X = matrix(1, nrow = N))
}

## sample sizes to consider
Ks <- c(seq(5, 50, by = 5), 100, 200)

## variance components to consider
truepsis <- seq(0, 0.1, length = 11)
## other simulation constants
truebeta <- 0.5

## compute X matrices and generate sigma2
sigma2 <- Xmat <- as.list(rep(0, length(Ks)))
## generate sigmas and set Xs
set.seed(123)
for (j in seq_along(Ks)) {
  sigma2[[j]] <- generate.sigma2s(Ks[j])
  Xmat[[j]] <- matrix(1, nrow = Ks[j])
}

## Specify number of cores, simulation size, bootstrap replications
registerDoMC(40)
nsimu <- 10000
nboot_ZL <- 1000

res <- list()

for(psi in c(0.03, 0.07)) {
  cat(paste("psi = ", psi), "\n")
  
  for(K in Ks) {
    cat(paste("K = ", K), "\n")
    
    set.seed(123)
    simu_data_mat <- matrix(NA, nrow = K, ncol = nsimu)
    
    for (i in seq.int(nsimu)) {
      dat <- simulate.BG(beta = truebeta, var.re = psi, sigma2s = 
                           sigma2[[which(Ks == K)]])
      simu_data_mat[, i] <- dat$betahat
    }
    
    simu_data <- data.frame(t(simu_data_mat))
    simu_data$sample <- seq.int(nsimu)
    
    ## Bootstrap seeds
    simu_data$seeds <- sample(seq.int(nsimu * 100), nsimu, 
                              replace = FALSE)
    
    attr(simu_data, "beta") <- truebeta
    attr(simu_data, "K") <- K
    attr(simu_data, "psi") <- psi
    attr(simu_data, "sigma2") <- sigma2[[which(Ks == K)]]
    attr(simu_data, "Nsim") <- nsimu
    
    temp_res <- dlply(simu_data, ~ sample, function(response) {
      temp_data <- unname(unlist(response[-(K + 1:2)]))
      
      ## Preliminary metaLik fit
      temp_fit <- metaLik(temp_data ~ 1, sigma2 = sigma2[[which(Ks == K)]])
      ## more accurate fit via Biasfit
      bias_fit <- BiasFit(temp_fit, method = "all")
      
      ## quantities for location adjustment
      # ML
      fitML <- fitfun(bias_fit$ML, temp_fit$X, temp_fit$sigma2, temp_fit$y)
      # BR
      fitBR <- fitfun(bias_fit$BR, temp_fit$X, temp_fit$sigma2, temp_fit$y)
      
      # Wald statistics
      zstat_ml <- waldi_meta(null = truebeta, which = 1, adjust = FALSE, 
                             method = "ML", fit = fitML)
      zstat_br <- waldi_meta(null = truebeta, which = 1, adjust = FALSE, 
                             method = "BR", fit = fitBR)
      zstat_ml_cor <- waldi_meta(null = truebeta, which = 1, adjust = TRUE, 
                                 method = "ML", fit = fitML)
      zstat_br_cor <- waldi_meta(null = truebeta, which = 1, adjust = TRUE, 
                                 method = "BR", fit = fitBR)
      # ZL and DL statistics
      stat_values <- perform_tests(temp_fit$y, temp_fit$X, Kval = K,
                                   temp_fit$sigma2, null = truebeta, 
                                   what = 1, B = nboot_ZL, silent = FALSE, 
                                   maxiter = 10000, seed = response$seeds)
      
      if (response$sample %% 100 == 0) cat(response$sample, "\n")
      stats <- data.frame(name = c("ml", "ml_cor", "br", "br_cor", "DL", "KH"),
                          value = c(zstat_ml, zstat_ml_cor, zstat_br, zstat_br_cor,
                                    stat_values["DL", 1], stat_values["KH", 1]))
      pvalues <- data.frame(statistic = c(rep(c("ZL", "KH"), times = 3)),
                            pvalue = c(stat_values["ZL", 2], stat_values["KH", 2],
                                       stat_values["ZL", 3], stat_values["KH", 3],
                                       stat_values["ZL", 4], stat_values["KH", 4]),
                            type = rep(c("pval_2sided", "pval_right", "pval_left"), 
                                       each = 2))
      list(id_sample = response$sample, sample = temp_data, stats = stats, 
           pvalues = pvalues, K = attr(simu_data, "K"), 
           psi = attr(simu_data, "psi"), beta = attr(simu_data, "beta"), 
           sigma2 = attr(simu_data, "sigma2"))                                       
    }, .parallel = TRUE)
    
    res <- c(res, temp_res)
    save.image(paste0("~/meta_psi=", psi, "K=", K, ".rda"))
    
  }   
}

## Summary of results
res_stats <- ldply(res, function(x) data.frame(x$stats, K = x$K, psi = x$psi))
res_pvals <- ldply(res, function(x) data.frame(x$pvalues, K = x$K, psi = x$psi))

typeI_statistics <- ddply(res_stats, ~ name + K + psi, function(x) {
  levels <- c(0.1, 1, 2.5, 5)/100
  p_value_2sided <- 2 * pnorm(-abs(x$value))
  p_value_left <- pnorm(x$value)
  p_value_right <- 1 - pnorm(x$value)
  rate_2sided <- sapply(levels, function(alpha) mean(p_value_2sided < alpha))
  rate_left <- sapply(levels, function(alpha) mean(p_value_left < alpha))
  rate_right <- sapply(levels, function(alpha) mean(p_value_right < alpha))
  out <- data.frame(
    test = rep(c("2sided", "left", "right"), each = length(levels)),
    typeI = c(rate_2sided, rate_left, rate_right),
    level = rep(levels, times = 3))
  out
})

typeI_pvalues <-
  ddply(res_pvals, ~ statistic + K + psi, function(x) {
    levels <- c(0.1, 1, 2.5, 5)/100
    rate_2sided <- sapply(levels, function(alpha) 
      mean(x$pvalue[x$type == 'pval_2sided'] < alpha))
    rate_left <- sapply(levels, function(alpha) 
      mean(x$pvalue[x$type == 'pval_left'] < alpha))
    rate_right <- sapply(levels, function(alpha) 
      mean(x$pvalue[x$type == 'pval_right'] < alpha))
    out <- data.frame(
      test = rep(c("2sided", "left", "right"), each = length(levels)),
      typeI = c(rate_2sided, rate_left, rate_right),
      level = rep(levels, times = 3))
    out
  })


names(typeI_statistics) <- names(typeI_pvalues)

typeIdf <- rbind(typeI_statistics, typeI_pvalues)

typeIdf <- typeIdf %>% 
  filter(is.element(statistic, c("ml", "ml_cor", "br", "br_cor",
                                 "DL", "ZL")), test == "2sided", level == 0.05) %>%
  mutate(level_chr = paste(level * 100, "~symbol('\045')"),
         cov = (1 - typeI) * 100,
         cov_chr = paste((1 - level) * 100, "~symbol('\045')"),
         upper = typeI - qnorm(1 - 0.01/2) * sqrt(typeI * (1 - typeI)/nsimu),
         lower = typeI + qnorm(1 - 0.01/2) * sqrt(typeI * (1 - typeI)/nsimu))

typeIdf$statistic <- factor(typeIdf$statistic, levels = c("ml", "ml_cor", 
                                                          "br", "br_cor", "DL", "ZL"), ordered = TRUE)
typeIdf$psilab <- paste0("psi = ", typeIdf$psi)

typeIdf$psilab <- factor(typeIdf$psilab, levels = c("psi = 0.03", "psi = 0.07"), labels =
                           c("psi == 0.03", "psi == 0.07"), ordered = TRUE)

cov_df_psi <- typeIdf

save(cov_df_K, cov_df_psi, file = paste(results_path, file = paste(path, 
         "results/meta_analysis_simulation.rda", sep = "/")))


