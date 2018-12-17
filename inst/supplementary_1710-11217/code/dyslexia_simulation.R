## devtools::install_github("ikosmidis/waldi")

library("waldi")
library("betareg")
library("enrichwith")
library("parallel")
library("dplyr")
library("doMC")
library("boot")

path <- '.'       # path; make sure that path has a directory named results
group <- 1        # group
n_cores <- 24     # number of cores
n_groups <- 1     # number of groups (same as in parallel command)
nsimu <- 50000    # number of simulated samples
R <- 500          # bootstrap size
eps <- 1e-12      # replace 0, 1 observations with eps or 1 - eps

## Split up the simulation in groups
groups <- split(seq.int(nsimu), rep(seq.int(n_groups), each = nsimu/n_groups))
## Seed
set.seed(123)

todf <- function(coverlist, sample) {
    ret <- lapply(names(coverlist), function(x) {
        co <- coverlist[[x]]
        out <- data.frame(parameter = rep(rownames(co), ncol(co)),
                          level = rep(as.numeric(colnames(co)), each = nrow(co)),
                          cover = c(co),
                          statistic = x,
                          sample = sample,
                          stringsAsFactors = FALSE)
    })
    do.call("rbind", ret)
}

dyslexia_bootstrap <- function(full_fit, R = 100, ncores = 1, levels = c(0.90, 0.95, 0.99), eps = 1e-12) {
  simulate_dysl <- get_simulate_function(full_fit)
  obs_data <- model.frame(full_fit$model)
  full_coef <- coef(full_fit)
  ncoef <- length(full_coef)
  coef_names <- names(full_coef)
  zstats <- function(data, psi) {
      temp <- obs_data
      ## replace boundary observations with something small
      data[abs(1 - data) < eps] <- 1 - eps
      data[abs(data) < eps] <- eps
      temp$accuracy <- data
      temp_fit <- try(betareg(accuracy ~ dyslexia * iq | dyslexia + iq, data = temp,
                              start = coef(full_fit)),
                      silent = TRUE)
      temp_fit_BR <- try(betareg(accuracy ~ dyslexia * iq | dyslexia + iq, data = temp,
                                 start = coef(full_fit), type = "BR"),
                         silent = TRUE)
      if (inherits(temp_fit, "try-error") | inherits(temp_fit_BR, "try-error")) {
        rep(NA, 2 * length(coef(full_fit)))
    }
    else {
      temp_coefs <- coef(temp_fit)
      zstat <- (temp_coefs - psi)/sqrt(diag(vcov(temp_fit))) # p x 1 vector
      zstat_cor <- waldi(temp_fit, null = psi, adjust = TRUE)
      c(zstat, zstat_cor)
    }
  }
  generate_dysl <- function(data, mle) {
      simulate_dysl(mle)
  }
  conv <- boot(obs_data$accuracy, statistic = zstats,
               R = R, sim = "parametric", ran.gen = generate_dysl,
               mle = full_coef, psi = full_coef,
               parallel = "multicore", ncpus = ncores)
  conv <- list(zstat = conv$t[, 1:ncoef],
               zstat_cor = conv$t[, ((ncoef + 1):(2*ncoef))])
  conv <- lapply(conv, function(x) {
      if (!is.matrix(x)) {
          x <- matrix(x, nrow = 1)
      }
      dimnames(x) <- list(NULL, coef_names)
      x
  })
  ## Bootstrap estimates of scale of the various statistics
  ## (na.rm = TRUE to protect from boundary observations)
  ses  <- apply(conv$zstat, 2, sd, na.rm = TRUE)
  ses_cor <- apply(conv$zstat_cor, 2, sd, na.rm = TRUE)
  ## Quantiles for studentized level % intervals
  ## (na.rm = TRUE to protect from boundary observations)
  a <- (1 - levels)/2
  a <- c(a, 1 - a)
  statistics <- c("zstat", "zstat_cor")
  quantiles <- lapply(statistics, function(stat) {
      out <- t(apply(conv[[stat]], 2, quantile, probs = a, na.rm = TRUE))
      colnames(out) <- a
      out
  })
  names(quantiles) <-  statistics
  list(ses = ses, ses_cor = ses_cor, quantiles = quantiles)
}


## Directly from ?betareg and as in Smythson and Verkuilen
data("ReadingSkills", package = "betareg")
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
                   data = ReadingSkills)

simulate_rs_beta <- get_simulate_function(rs_beta)

## Simulation study
simu_data <- simulate_rs_beta(nsim = nsimu)

## Get number of observations and number of parameters
n_obs <- nrow(simu_data)
n_par <- length(coef(rs_beta))

simu_data <- data.frame(t(simu_data))

levels  <-  c(0.90, 0.95, 0.99)
a <- (1 - levels)/2
a <- c(a, 1 - a)

coefs <- coef(rs_beta)
cnames <- names(coefs)

t_start <- proc.time()
res <- mclapply(groups[[group]], function(k) {
    temp_data <- ReadingSkills
    cdat <- unlist(simu_data[k, ])
    cdat[abs(1 - cdat) < eps] <- 1 - eps
    cdat[abs(cdat) < eps] <- eps
    temp_data$accuracy <- cdat
    temp_fit_ml <- try(betareg(accuracy ~ dyslexia * iq | dyslexia + iq, data = temp_data,
                               start = coef(rs_beta)),
                       silent = TRUE)
    temp_fit_br <- try(betareg(accuracy ~ dyslexia * iq | dyslexia + iq, data = temp_data,
                               start = coef(rs_beta), type = "BR"),
                       silent = TRUE)
    ## If any of the fits fails then return NAs
    if (inherits(temp_fit_ml, "try-error") | inherits(temp_fit_br, "try-error")) {
        cover <- as.list(numeric(12))
        cover <- lapply(1:12, function(x) {
            matrix(NA, length(coefs), length(levels), dimnames = list(cnames, levels))
        })
        names(cover) <- c("ml", "br", "ml_cor", "br_cor", "ml_stud", "ml_cor_stud", "br_stud", "br_cor_stud", "ml_cor_ses", "ml_cor_ses_cor", "br_cor_ses", "br_cor_ses_cor")
    }
    else {
        zstat_ml_cor <- waldi(temp_fit_ml, null = coefs, adjust = TRUE)
        temp_coefs <- coef(temp_fit_ml)
        zstat_ml <- (temp_coefs - coefs)/sqrt(diag(vcov(temp_fit_ml)))
        zstat_br <- (coef(temp_fit_br) - coefs)/sqrt(diag(vcov(temp_fit_br)))
        zstat_br_cor <- waldi(temp_fit_br, null = coefs, adjust = TRUE)
        ## Same seed for all bootstraps
        set.seed(111)
        boot_results <- try(dyslexia_bootstrap(temp_fit_ml, R = R, ncores = 1, eps = eps), silent = TRUE)
        ## Same seed for all bootstraps
        set.seed(111)
        boot_results_br <- try(dyslexia_bootstrap(temp_fit_br, R = R, ncores = 1, eps = eps), silent = TRUE)
        ## scale adjusted versions of location adjusted wald statistics (which are infisible though)
        ## Adjust the scale of zstat_mle_cor using the bootstrap estimate of the sd of zstat_mle
        zstat_ml_cor_ses <- zstat_ml_cor/boot_results$ses
        ## Adjust the scale of zstat_mle_cor using the bootstrap estimate of the sd of zstat_mle_cor
        zstat_ml_cor_ses_cor <- zstat_ml_cor/boot_results$ses_cor
        ## Adjust the scale of zstat_br_cor using the bootstrap estimate of the sd of zstat_br
        zstat_br_cor_ses <- zstat_br_cor/boot_results_br$ses
        ## Adjust the scale of zstat_br_cor using the bootstrap estimate of the sd of zstat_br_cor
        zstat_br_cor_ses_cor <- zstat_br_cor/boot_results_br$ses_cor
        stats_normal <- cbind(ml = zstat_ml,
                              br = zstat_br,
                              ml_cor = zstat_ml_cor,
                              br_cor = zstat_br_cor)
        stats_scaled <- cbind(ml_cor_ses = zstat_ml_cor_ses,
                              ml_cor_ses_cor = zstat_ml_cor_ses_cor,
                              br_cor_ses = zstat_br_cor_ses,
                              br_cor_ses_cor = zstat_br_cor_ses_cor)
        stats_stud_ml <- stats_normal[, c("ml", "ml_cor")]
        colnames(stats_stud_ml) <- c("ml_stud", "ml_cor_stud")
        stats_stud_br <- stats_normal[, c("br", "br_cor")]
        colnames(stats_stud_br) <- c("br_stud", "br_cor_stud")
        normq <- qnorm(1 - (1-levels)/2)
        ## Normal
        cover_normal <- lapply(colnames(stats_normal), function(interval) {
            out <- sapply(normq, function(qq) abs(stats_normal[, interval]) < qq)
            colnames(out) <- levels
            out
        })
        names(cover_normal) <- colnames(stats_normal)
        ## Scaled
        cover_scaled <- lapply(colnames(stats_scaled), function(interval) {
            out <- sapply(normq, function(qq) abs(stats_scaled[, interval]) < qq)
            colnames(out) <- levels
            out
        })
        names(cover_scaled) <- colnames(stats_scaled)
        ## Studentized
        quantiles_ml <- boot_results$quantiles
        quantiles_br <- boot_results_br$quantiles
        cover_stud_ml <- cover_stud_br <- cover_normal[1:2]
        names(cover_stud_ml) <- colnames(stats_stud_ml)
        names(cover_stud_br) <- colnames(stats_stud_br)
        for (j in levels) {
            low <- as.character((1 - j)/2)
            upp <- as.character(1 - (1 - j)/2)
            lev <- as.character(j)
            ## ml
            qlow <- quantiles_ml$zstat[cnames, low]
            qupp <- quantiles_ml$zstat[cnames, upp]
            st <- stats_stud_ml[cnames, "ml_stud"]
            cover_stud_ml$ml_stud[cnames, lev] <- (st > qlow) & (st < qupp)
            ## ml_cor``xs
            qlow <- quantiles_ml$zstat_cor[cnames, low]
            qupp <- quantiles_ml$zstat_cor[cnames, upp]
            st <- stats_stud_ml[cnames, "ml_cor_stud"]
            cover_stud_ml$ml_cor_stud[cnames, lev] <- (st > qlow) & (st < qupp)
            ## br
            qlow <- quantiles_br$zstat[cnames, low]
            qupp <- quantiles_br$zstat[cnames, upp]
            st <- stats_stud_br[cnames, "br_stud"]
            cover_stud_br$br_stud[cnames, lev] <- (st > qlow) & (st < qupp)
            ## br_cor
            qlow <- quantiles_br$zstat_cor[cnames, low]
            qupp <- quantiles_br$zstat_cor[cnames, upp]
            st <- stats_stud_br[cnames, "br_cor_stud"]
            cover_stud_br$br_cor_stud[cnames, lev] <- (st > qlow) & (st < qupp)
        }
        cover_stud <- c(cover_stud_ml, cover_stud_br)
        cover <- c(cover_normal, cover_stud, cover_scaled)
    }
    if (k %% 1 == 0) {
        cat("group:", group, "|", n_groups, "sample:", k, "|", nsimu, "||", format(Sys.time()), "\n")
    }
    todf(cover, sample = k)
}, mc.cores = n_cores)
res <- do.call("rbind", res)
t_end <- proc.time()

to_be_saved <- c("nsimu", "t_start", "t_end",
                 "rs_beta", "simulate_rs_beta",
                 "n_obs", "n_par",
                 "dyslexia_bootstrap",
                 "res")
if (group == 1) {
    to_be_saved  <- c(to_be_saved, "simu_data")
}

save(list = to_be_saved,
     file = paste(path, paste0("results/dyslexia_simulation.rda"), sep = "/"))
