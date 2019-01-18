## devtools::install_github("ikosmidis/waldi")
## Specify path; make sure that path has a directory named results
path <- "."

library("waldi")
library("enrichwith")
library("brglm2")
library("plyr")
library("dplyr")
library("doMC")

## Specify number of cores, simulation size and set seed
registerDoMC(20)
nsimu <- 1000
nboot <- 1000
eps <- 1e-04
set.seed(123)

## Bootstrap function
clotting_bootstrap <- function(full_fit, R = 100, ncores = 1) {
  simulate_clot <- get_simulate_function(full_fit)
  obs_data <- full_fit$data
  full_coef <- coef(full_fit)
  ncoef <- length(full_coef)
  
  zstats_ml <- function(data, psi) {
    temp <- obs_data
    temp$conc <- unlist(data)
    temp_fit <- try(update(full_fit, data = temp))
    if (inherits(temp_fit, "try-error")) {
      zstat_ml <- zstat_ml_cor <- rep(NA, ncoef)
      c(zstat_ml, zstat_ml_cor)
    }
    else {
      temp_coefs <- coef(temp_fit)
      zstat_ml <- waldi(temp_fit, null = psi, adjust = FALSE, numerical = FALSE)
      zstat_ml_cor <- waldi(temp_fit, null = psi, adjust  = TRUE, numerical = FALSE)
      c(zstat_ml, zstat_ml_cor)
    }
  }
  generate_clot <- function(data, mle) {
    simulate_clot(mle)
  }
  
  zstats_mom <- function(data, psi) {
    temp <- obs_data
    temp$conc <- unlist(data)
    temp_fit <- try(update(full_fit, data = temp))
    if (inherits(temp_fit, "try-error")) {
      zstat_mom <- rep(NA, ncoef)
    }
    else {
      temp_coefs <- coef(temp_fit)
      (temp_coefs - psi)/sqrt(diag(vcov(temp_fit)))
    }
  }
  generate_clot_mom <- function(data, mle) {
    simulate_clot(mle, dispersion = summary(full_fit)$dispersion)
  }
  
  boot_ml <- boot(obs_data$conc, statistic = zstats_ml,
                  R = R, sim = "parametric", ran.gen = generate_clot,
                  mle = full_coef, psi = full_coef,
                  parallel = "multicore", ncpus = ncores)
  boot_mom <- boot(obs_data$conc, statistic = zstats_mom, R = R,
                   sim = "parametric", ran.gen = generate_clot_mom,
                   mle = full_coef, psi = full_coef,
                   parallel = "multicore", ncpus = ncores)
  
  data.frame(zstat_ml = boot_ml$t[, 1:ncoef], zstat_ml_cor = 
               boot_ml$t[, ((ncoef + 1):(2*ncoef))], zstat_mom = boot_mom$t)
}

## The clotting data set
clotting <- data.frame(
  conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
  u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
  lot = factor(c(rep(1, 9), rep(2, 9))))

## The maximum likelihood fit of the gamma regression model
model_fit <- glm(conc ~ log(u)*lot, data = clotting, family = Gamma(link = "log"))
summary(model_fit)

## Simulate at the maximum likelihood estimates of beta and phi
clotting_simulate <- get_simulate_function(model_fit)
simu_data <- clotting_simulate(nsim = nsimu)

## Get number of observations and number of parameters
n_obs <- nrow(simu_data)
n_par <- length(coef(model_fit))

simu_data <- data.frame(t(simu_data))
simu_data$sample <- seq.int(nsimu)
## Bootstrap seeds
simu_data$seeds <- sample(seq.int(nsimu * 100), nsimu, replace = FALSE)

start <- seq(0, 9, 1) * 100 + 1
end <- seq(1, 10, 1) * 100

res <- list()

## Perform simulation
for (j in seq_along(start)) {
  cdat <- simu_data[start[j]:end[j], ]
  temp_res <- dlply(cdat, ~ sample, function(response) {
    temp_data <- clotting
    temp_data$conc <- unlist(response[-(n_obs + 1:2)])
    temp_ml_fit <- update(model_fit, data = temp_data)
    temp_br_fit <- update(model_fit, data = temp_data, method = "brglmFit")
    zstat_ml <- waldi(temp_ml_fit, null = coef(model_fit), adjust = FALSE, 
                      numerical = FALSE)
    zstat_rb <- waldi(temp_br_fit, null = coef(model_fit), adjust = FALSE, 
                      numerical = FALSE)
    zstat_ml_cor <- waldi(temp_ml_fit, null = coef(model_fit), adjust = TRUE, 
                          numerical = FALSE)
    zstat_rb_cor <- waldi(temp_br_fit, null = coef(model_fit), adjust = TRUE, 
                          numerical = FALSE)
    zstat_mom <- (coef(temp_ml_fit) - coef(model_fit))/sqrt(diag(vcov(temp_ml_fit)))
    zstats_temp <- c(zstat_ml, zstat_ml_cor, zstat_mom)
    ### bootstrap
    set.seed(temp_data$seed)
    boots_conv <- clotting_bootstrap(temp_ml_fit, R = nboot, ncores = 1)
    conv_mat <- rbind(data.matrix(boots_conv), zstats_temp)
    p_2sided_conv <- apply(conv_mat, 2, function(col) (sum(abs(col[-(nboot + 1)]) >= 
                                                             abs(col[nboot + 1])) + eps)/(nboot + 2 * eps))
    p_left_conv <- apply(conv_mat, 2, function(col) (sum(col[-(nboot + 1)] <= 
                                                           col[nboot + 1]) + eps)/(nboot + 2 * eps))
    p_right_conv <- apply(conv_mat, 2, function(col) (sum(col[-(nboot + 1)] >= 
                                                            col[nboot + 1]) + eps)/(nboot + 2 * eps))
    
    ### location- and scale-adjusted Wald statistics
    ## variance estimated by conventional bootstrap
    ses_ml_conv <- apply(boots_conv[, 1:n_par], 2, sd)
    ses_ml_cor_conv <- apply(boots_conv[, (n_par + 1):(2*n_par)], 2, sd)
    
    zstat_loc_sc_ml_conv <- zstat_ml_cor/ses_ml_conv   ## p x 1 vector
    zstat_loc_sc_ml_cor_conv <- zstat_ml_cor/ses_ml_cor_conv
    
    if (response$sample %% 100 == 0)
      cat(response$sample, "\n")
    stats <- data.frame(name = rep(c("ml", "ml_cor", "rb", "rb_cor", "mom", 
                                     "sc_ml_conv", "sc_ml_cor_conv"), each = n_par),
                        value = c(zstat_ml, zstat_ml_cor, zstat_rb, zstat_rb_cor, zstat_mom, 
                                  zstat_loc_sc_ml_conv, zstat_loc_sc_ml_cor_conv),
                        parameter = rep(1:n_par, times = 7))
    boot_pvalues <- data.frame(statistic = c(rep(rep(c("ml", "ml_cor", "mom"), 
                                                     each = n_par), times = 3)),
                               value = c(p_2sided_conv, p_right_conv, p_left_conv),
                               type = rep(c("boot_conv_2sided", "boot_conv_right", 
                                            "boot_conv_left"), each = 3*n_par),
                               parameter = rep(rep(1:n_par, times = 3), times = 3))
    list(stats = stats, boot_pvalues = boot_pvalues)                                           
  }, .parallel = TRUE)
  res <- c(res, temp_res)
  save.image(paste0("~/clotting_simu", end[j], ".rda"))
}

## Summary of results
res_statistics <- ldply(res, function(x) x$stats)
res_pvalues <- ldply(res, function(x) x$boot_pvalues)

## Type I error rates
typeI_statistics <- ddply(res_statistics, ~ name + parameter, function(x) {
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

typeI_pvalues <- ddply(res_pvalues, ~ statistic + parameter, function(x) {
  levels <- c(0.1, 1, 2.5, 5)/100
  rate_2sided <- sapply(levels, function(alpha) mean(x$value[x$type == 
                                                               'boot_conv_2sided'] < alpha))
  rate_left <- sapply(levels, function(alpha) mean(x$value[x$type == 
                                                             'boot_conv_left'] < alpha))
  rate_right <- sapply(levels, function(alpha) mean(x$value[x$type == 
                                                              'boot_conv_right'] < alpha))
  out <- data.frame(
    test = rep(c("2sided", "left", "right"), each = length(levels)),
    typeI = c(rate_2sided, rate_left, rate_right),
    level = rep(levels, times = 3))
  out
})

names(typeI_statistics) <- names(typeI_pvalues)
levels(typeI_pvalues$statistic) <- c("ml_boot", "ml_cor_boot", "mom_boot")

typeI <- rbind(typeI_statistics, typeI_pvalues)

typeI <- typeI %>%
  filter(test != "right") %>%
  mutate(test = recode(test,
                       "2sided" = "H[1]: beta[italic(j)] != beta[paste(italic(j), 0)]",
                       "left" = "H[1]: beta[italic(j)] < beta[paste(italic(j), 0)]",
                       "right" = "H[1]: beta[italic(j)] > beta[paste(italic(j), 0)]"),
         level_chr = paste(level*100, "~symbol('\045')"),
         upper = typeI - qnorm(1 - 0.01/2)*sqrt(typeI*(1-typeI)/nsimu),
         lower = typeI + qnorm(1 - 0.01/2)*sqrt(typeI*(1-typeI)/nsimu))

save(typeI, file = paste(path, "results/clotting_simulation.rda", sep = "/"))