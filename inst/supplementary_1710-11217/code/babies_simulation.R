## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE!  Provided "as is",
## licensed under GPL2 or higher.  NO WARRANTY OF FITNESS FOR ANY
## PURPOSE!
##
## Ioannis Kosmidis [aut], i.kosmidis@ucl.ac.uk
##
## 18 April 2018

## Specify path; make sure that path has a directory named results
path <- "."

require("waldi")
require("brglm2")
require("enrichwith")
require("survival")
require("plyr")
require("doMC")
require("lmtest")
require("brglm2")
require("boot")
require("cond")

## Bootstrap function

babies_bootstrap <- function(full_fit, R = 100, eps = 1e-04, ncores = 1) {
    simulate_babies <- get_simulate_function(full_fit)
    obs_data <- model.frame(full_fit$model)
    null_fit <- update(full_fit, . ~ . - lull, data = obs_data)
    full_coef <- coef(full_fit)
    null_coef <- c(coef(null_fit), 0)
    obs_z <- summary(full_fit)$coefficients["lullyes", 3]
    z_statistic <- function(data, psi) {
        temp <- obs_data
        temp$y <- data
        temp_fit <- update(full_fit, data = temp)
        est <- summary(temp_fit)$coef["lullyes", c(1, 2)]
        is_inf <- is.infinite(update(temp_fit, method = "detect_separation")$beta["lullyes"])
        c(z = (est[1] - psi)/est[2], infinite = is_inf)
    }
    generate_babies <- function(data, mle) {
        simulate_babies(mle)[[1]]
    }
    conv <- boot(obs_data$y, statistic = z_statistic,
                 R = R, sim = "parametric", ran.gen = generate_babies,
                 mle = full_coef, psi = full_coef["lullyes"],
                 parallel = "multicore", ncpus = ncores)
    prep <- boot(obs_data$y, statistic = z_statistic,
                 R = R, sim = "parametric", ran.gen = generate_babies,
                 mle = null_coef, psi = 0,
                 parallel = "multicore", ncpus = ncores)
    list(conv = (sum(abs(conv$t[, 1]) >= abs(obs_z)) + eps)/(R + 2 * eps), # 2 * min(mean(conv$t <= obs_z), mean(conv$t > obs_z)),
         prep = (sum(abs(prep$t) >= abs(obs_z)) + eps)/(R + 2 * eps), # 2 * min(mean(prep$t <= obs_z), mean(prep$t > obs_z))
         conv_left = (sum(conv$t[, 1] <= obs_z) + eps)/(R + 2 * eps),
         prep_left = (sum(prep$t[, 1] <= obs_z) + eps)/(R + 2 * eps),
         conv_right = (sum(conv$t[, 1] >= obs_z) + eps)/(R + 2 * eps),
         prep_right = (sum(prep$t[, 1] >= obs_z) + eps)/(R + 2 * eps),
         conv_inf = sum(conv$t[, 2]),
         prep_inf = sum(prep$t[, 2]))
}

## Specify number of cores, simulation size and set seed
registerDoMC(20)
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
## Bootstrap seeds
simu_data$seeds <- sample(seq.int(nsimu * 100), nsimu, replace = FALSE)

## Perform simulation
res <- ddply(simu_data, ~ sample, function(response) {
    temp_data <- babies_expand
    temp_data$y <- unlist(response[-(n_obs + 1:2)])
    temp_ml <- update(babies_ml, data = temp_data)
    temp_rb <- update(babies_ml, data = temp_data, method = "brglmFit")
    temp_cond <- update(babies_cond, data = temp_data)
    ## Wald
    zstat_cor <- waldi(temp_ml, null = 0, adjust = TRUE, which = 19)
    zstat_rb_cor <- waldi(temp_rb, null = 0, adjust = TRUE, which = 19)
    mle <- coeftest(temp_ml)["lullyes", ]
    rbe <- coeftest(temp_rb)["lullyes", ]
    mcle <- coeftest(temp_cond)["lullyes", ]
    zstat_mle <- mle["z value"]
    zstat_rbe <- rbe["z value"]
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
    set.seed(temp_data$seed)
    boots <- unlist(babies_bootstrap(temp_ml, R = 1000, ncores = 1, eps = 0.1/1000))
    if (response$sample %% 10 == 0)
        cat(response$sample, "\n")
    data.frame(name = c("mle", "rbe", "cond", "cor", "cor_rb", "r", "rc", "scorec",
                        "boot_conv_2sided", "boot_conv_right", "boot_conv_left",
                        "boot_prep_2sided", "boot_prep_right", "boot_prep_left",
                        "conv_inf", "prep_inf"),
               value = unname(c(zstat_mle, zstat_rbe, zstat_cond, zstat_cor, zstat_rb_cor, r, rc, scorec,
                                qnorm(boots[c("conv", "conv_right", "conv_left", "prep", "prep_right", "prep_left")]),
                                boots[c("conv_inf", "prep_inf")])),
               type = c(rep("statistic", 8), rep("bootstrap_statistic", 6), rep("summary", 2)),
               parameter = rep(1, times = 16),
               infinite = unname(is_inf))
}, .parallel = TRUE)


## Count and replace mle and cor with zero when gamma is infinite
sum(res$infinite)/nlevels(res$name)
## [1] 13
res$value <- ifelse(res$infinite & (res$name %in% c("mle", "cor")), 0, res$value)

save(nsimu, babies_expand, babies_ml, babies_cond, simulate_babies, simu_data, n_obs, n_par, res, babies_bootstrap,
     file = paste(path, "results/babies_simulation.rda", sep = "/"))





