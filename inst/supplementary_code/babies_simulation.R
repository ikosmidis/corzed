## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE!  Provided "as is",
## licensed under GPL2 or higher.  NO WARRANTY OF FITNESS FOR ANY
## PURPOSE!
##
## Ioannis Kosmidis [aut], i.kosmidis@ucl.ac.uk
##
## 26 October 2017

## Specify path; make sure that path has a directory named results
path <- "."
source(paste(path, "corzed.R", sep = "/"))

require("enrichwith")
require("survival")
require("plyr")
require("dplyr")
require("doMC")
require("lmtest")
require("brglm2")
require("ggplot2")

## Specify number of cores, simulation size and set seed
registerDoMC(18)
nsimu <- 50000
set.seed(123)

data("babies", package = "cond")

## clogit understands only 0-1 so expand
babies_expand <- ddply(babies, ~ lull + day, function(z) {
    data.frame(y = rep(c(0, 1), c(z$r2, z$r1)))
})


babies_ml <- glm(formula = y ~ day + lull - 1,
                  family = binomial, data = babies_expand)
babies_cond <- clogit(y ~ strata(day) + lull, data = babies_expand)

## Simulation study
simulate_babies <- get_simulate_function(babies_ml)
simu_data <- simulate_babies(coefficients = c(coef(babies_ml)[-19], 0),
                             nsim = nsimu)

## Get number of observations and number of parameters
n_obs <- nrow(simu_data)
n_par <- length(coef(babies_ml))

simu_data <- data.frame(t(simu_data))
simu_data$sample <- seq.int(nsimu)

## Perform simulation
res <- ddply(simu_data, ~ sample, function(response) {
    temp_data <- babies_expand
    temp_data$y <- unlist(response[-(n_obs + 1)])
    temp_ml <- update(babies_ml, data = temp_data)
    temp_cond <- update(babies_cond, data = temp_data)
    ## Wald
    zstat_cor <- corzed(temp_ml, null = 0, correction = TRUE, what = 19)
    mle <- coeftest(temp_ml)["lullyes", ]
    mcle <- coeftest(temp_cond)["lullyes", ]
    zstat_mle <- mle["z value"]
    zstat_cond <- mcle["z value"]
    is_inf <- is.infinite(update(temp_ml, method = "detect_separation")$beta["lullyes"])
    ## LR
    r <- lrtest(update(temp_ml, . ~ . - lull), temp_ml)
    temp_cond_sum <- summary(temp_cond)
    rc <- temp_cond_sum$logtest[1]
    scorec <- temp_cond_sum$sctest[1]
    r <- unname(sign(mle["Estimate"]) * sqrt(r$Chisq[2]))
    rc <- unname(sign(mcle["Estimate"]) * sqrt(rc))
    scorec <- unname(sign(mcle["Estimate"]) * sqrt(scorec))
    if (response$sample %% 100 == 0)
        cat(response$sample, "\n")
    data.frame(statistic = rep(c("mle", "cond", "cor", "r", "rc", "scorec"), each = 1),
               value = c(zstat_mle, zstat_cond, zstat_cor, r, rc, scorec),
               parameter = rep(1, times = 6),
               infinite = is_inf)
}, .parallel = TRUE)


## Count and replace cor with zero when gamma is infinite
sum(res$infinite)/nlevels(res$statistic)
## [1] 13
res$value <- ifelse(res$infinite & res$statistic == "cor", 0, res$value)

save(nsimu, babies_expand, babies_ml, babies_cond, simulate_babies, simu_data, n_obs, n_par, res,
     file = paste(path, "results/babies_simulation.rda", sep = "/"))
