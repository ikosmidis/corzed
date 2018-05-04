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

## devtools::install_github("ikosmidis/waldi")
library("waldi")
library("enrichwith")
library("brglm2")
library("plyr")
library("dplyr")
library("doMC")


## Specify path; make sure that path has a directory named results
path <- "."

## Specify number of cores, simulation size and set seed
registerDoMC(3)
nsimu <- 50000
set.seed(123)

## The clotting data set
clotting <- data.frame(
  conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
  u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
  lot = factor(c(rep(1, 9), rep(2, 9))))

## The maximum likelihood fit of the gamma regression model
model_fit <- glm(conc ~ log(u)*lot, data = clotting, family = Gamma(link = "log"))

## Simulate at the maximum likelihood estimates of beta and phi
clotting_simulate <- get_simulate_function(model_fit)
simu_data <- clotting_simulate(nsim = nsimu)

## Get number of observations and number of parameters
n_obs <- nrow(simu_data)
n_par <- length(coef(model_fit))

simu_data <- data.frame(t(simu_data))
simu_data$sample <- seq.int(nsimu)

res <- ddply(simu_data, ~ sample, function(response) {
    temp_data <- clotting
    temp_data$conc <- unlist(response[-(n_obs + 1)])
    temp_ml_fit <- update(model_fit, data = temp_data)
    temp_br_fit <- update(model_fit, data = temp_data, method = "brglmFit")
    zstat_ml <- waldi(temp_ml_fit, null = coef(model_fit), adjust = FALSE, numerical = FALSE)
    zstat_ml_cor <- waldi(temp_ml_fit, null = coef(model_fit), adjust  = TRUE, numerical = FALSE)
    zstat_rb <- waldi(temp_br_fit, null = coef(model_fit), adjust  = FALSE, numerical = FALSE)
    zstat_rb_cor <- waldi(temp_br_fit, null = coef(model_fit), adjust  = TRUE, numerical = FALSE)
    zstat_mom <- (coef(temp_ml_fit) - coef(model_fit))/sqrt(diag(vcov(temp_ml_fit)))
    if (response$sample %% 100 == 0)
        cat(response$sample, "\n")
    data.frame(statistic = rep(c("ml", "ml_cor", "rb", "rb_cor", "mom"), each = n_par),
               value = c(zstat_ml, zstat_ml_cor, zstat_rb, zstat_rb_cor, zstat_mom),
               parameter = rep(1:4, times = 5))
}, .parallel = TRUE)

save(nsimu, clotting, model_fit, clotting_simulate, simu_data, n_obs, n_par, simu_data, res,
     file = paste(path, "results/clotting_simulation.rda", sep = "/"))

