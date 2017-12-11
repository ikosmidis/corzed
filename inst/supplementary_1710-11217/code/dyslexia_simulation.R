## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE!  Provided "as is",
## licensed under GPL2 or higher.  NO WARRANTY OF FITNESS FOR ANY
## PURPOSE!
##
## Ioannis Kosmidis [aut], i.kosmidis@ucl.ac.uk
##
## 26 October 2017
##
## This code still has lots of rough edges, some of which are
## indicated by the embedded comments.

## Specify path; make sure that path has a directory named results
path <- "."
source(paste(path, "corzed.R", sep = "/"))

library("betareg")
library("enrichwith")
library("plyr")
library("dplyr")
library("doMC")
library("lmtest")
library("brglm2")
library("ggplot2")

## Specify number of cores, simulation size and set seed
registerDoMC(3)
nsimu <- 50000
set.seed(123)


## Directly from ?betareg and as in Smythson and Verkuilen
data("ReadingSkills", package = "betareg")
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
                   data = ReadingSkills)
summary(rs_beta)

## data("GasolineYield", package = "betareg")
## rs_beta <- betareg(yield ~ batch + temp, data = GasolineYield)

simulate_rs_beta <- get_simulate_function(rs_beta)

## Simulation study
simu_data <- simulate_rs_beta(nsim = nsimu)

## Get number of observations and number of parameters
n_obs <- nrow(simu_data)
n_par <- length(coef(rs_beta))

simu_data <- data.frame(t(simu_data))
simu_data$sample <- seq.int(nsimu)

res <- ddply(simu_data, ~ sample, function(response) {
    coefs <- coef(rs_beta)
    temp_data <- ReadingSkills
    temp_data$accuracy <- unlist(response[-(n_obs + 1)])
    ## temp_data <- GasolineYield
    ## temp_data$yield <- unlist(response[-(n_obs + 1)])
    temp_fit <- try(update(rs_beta, data = temp_data))
    if (inherits(temp_fit, "try-error")) {
        zstat_cor <- zstat_mle <- zstat_mle_obs <- zstat_bc <- rep(NA, length(coef(rs_beta)))
    }
    else {
        zstat_cor <- corzed(temp_fit, null = coefs, correction = TRUE)
        temp_fit_hess <- update(temp_fit, hessian = TRUE, start = coef(temp_fit))
        ## zstat_cor_obs <- corzed(temp_fit_hess, null = coefs, correction = TRUE, use_observed = TRUE)
        temp_coefs <- coef(temp_fit)
        ## MLE + expected
        zstat_mle <- (temp_coefs - coefs)/sqrt(diag(vcov(temp_fit)))
        ## MLE + observed
        zstat_mle_obs <- (temp_coefs - coefs)/sqrt(diag(vcov(temp_fit_hess)))
        temp_fit_BC <- update(temp_fit, type = "BC")
        zstat_bc <- (coef(temp_fit_BC) - coefs)/sqrt(diag(vcov(temp_fit_BC)))
    }
    if (response$sample %% 100 == 0)
        cat(response$sample, "\n")
    data.frame(
        statistic = rep(c("mle", "mle+obs", "bc", "cor"), each = n_par),
        value = c(zstat_mle, zstat_mle_obs, zstat_bc, zstat_cor),
        parameter = rep(seq.int(n_par), times = 4))
}, .parallel = TRUE)

save(nsimu, rs_beta, simulate_rs_beta, simu_data, n_obs, n_par, simu_data, res,
     file = paste(path, "results/dyslexia_simulation.rda", sep = "/"))

