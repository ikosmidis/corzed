## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#", fig.width=7, fig.height=4, 
                      out.width="100%", dpi = 200)
options(width = 80)

## ----path----------------------------------------------------------------
path <- "."
code_path <- paste(path, "code", sep = "/")
results_path <- paste(path, "results", sep = "/")
lesions_path <- paste(path, "lesion data", sep = "/")

## ------------------------------------------------------------------------
dir(code_path)
dir(results_path)
dir(lesions_path)

## ----waldi---------------------------------------------------------------
waldi_version <- try(packageVersion("waldi"), silent = TRUE)
if (inherits(waldi_version, "try-error")) {
    devtools::install_github("ikosmidis/waldi")
}

## ----libraries, message = FALSE, warning = FALSE-------------------------
library("waldi")
library("oro.nifti")
library("boot")
library("plyr")
library("plotrix")
library("dplyr")
library("survival")
library("cond")
library("lmtest")
library("betareg")
library("enrichwith")
library("brglm2")
library("ggplot2")
library("gridExtra")
library("colorspace")

## ----dyslexia------------------------------------------------------------
data("ReadingSkills", package = "betareg")
## maximum likelihood estimates and corresponding 95\% Wald confidence intervals
rs_beta_ml <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
                      data = ReadingSkills, type = "ML", hessian = FALSE)
rs_summary_ml <- coef(summary(rs_beta_ml))
rs_ml_estimates <- do.call("rbind", lapply(rs_summary_ml,
                                           function(z) z[, c("Estimate", "Std. Error")]))
rs_ml_cis <- confint(rs_beta_ml)
## bias corrected fit and corresponding 95\% Wald confidence intervals
rs_beta_br <- update(rs_beta_ml, type = "BR")
rs_summary_br <- coef(summary(rs_beta_br))
rs_br_estimates <- do.call("rbind", lapply(rs_summary_br,
                                           function(z) z[, c("Estimate", "Std. Error")]))
rs_br_cis <- confint(rs_beta_br)
round(cbind(rs_ml_estimates, rs_br_estimates, rs_ml_cis, rs_br_cis), 3)

## ------------------------------------------------------------------------
load(paste(results_path, paste0("dyslexia_simulation.rda"), sep = "/"))
rs_coverage <- results %>%
    filter(parameter %in% c("dyslexia", "iq", "dyslexia:iq", "(phi)_dyslexia", 
                            "(phi)_iq")) %>%
    filter(statistic %in% c("ml", "br", "ml_cor", "br_cor",
                            "ml_stud", "ml_cor_stud", "br_stud", "br_cor_stud",
                            "ml_cor_ses_cor", "br_cor_ses_cor")) %>%
    mutate(parameter = recode(parameter,
                              "dyslexia" = 2, "iq" = 3, "dyslexia:iq" = 4,
                              "(phi)_dyslexia" = 6, "(phi)_iq" = 7)) %>%
    mutate(level = 100 * level) %>%
    group_by(level, statistic, parameter) %>%
    summarize(coverage = round(mean(cover, na.rm = TRUE) * 100, 1)) %>%
    as.data.frame() %>%
    reshape(idvar = c("statistic", "parameter"), v.names = "coverage",
            timevar = "level",
            direction = "wide")
rs_coverage %>% filter(statistic %in% c("ml", "br")) %>%
    select(statistic, parameter, coverage.90, coverage.95, coverage.99)

## ---- warning = TRUE-----------------------------------------------------
rs_cor_ml_cis <- waldi_confint(rs_beta_ml, level = 0.95, adjust = TRUE)
interpolation <- waldi_confint(rs_beta_ml, level = 0.95,
                               which = rownames(rs_cor_ml_cis),
                               adjust = TRUE,
                               return_values = TRUE,
                               length = 20)
intervals <- data.frame(low = rs_cor_ml_cis[, 1],
                        upp = rs_cor_ml_cis[, 2],
                        parameter = rownames(rs_cor_ml_cis))
interpolation <- interpolation %>%
    filter(!(parameter %in% c("(Intercept)", "(phi)_(Intercept)"))) %>%
    mutate(parameter = recode(parameter,
                              "dyslexia" = "beta[2]",
                              "iq" = "beta[3]",
                              "dyslexia:iq" = "beta[4]",
                              "(phi)_dyslexia" = "gamma[2]",
                              "(phi)_iq" = "gamma[3]"))
intervals <- intervals %>%
    filter(!(parameter %in% c("(Intercept)", "(phi)_(Intercept)"))) %>%
    mutate(parameter = recode(parameter,
                              "dyslexia" = "beta[2]",
                              "iq" = "beta[3]",
                              "dyslexia:iq" = "beta[4]",
                              "(phi)_dyslexia" = "gamma[2]",
                              "(phi)_iq" = "gamma[3]"))
ggplot(interpolation) +
    geom_point(aes(x = grid, y = value)) +
    geom_line(aes(x = grid, y = value), col = "grey") +
    geom_hline(aes(yintercept = qnorm(0.975)), col = "grey", lty = 3) +
    geom_hline(aes(yintercept = qnorm(0.025)), col = "grey", lty = 3) +
    geom_vline(data = intervals, aes(xintercept = low), col = "grey", lty = 2) +
    geom_vline(data = intervals, aes(xintercept = upp), col = "grey", lty = 2) +
    facet_grid(~ parameter, scale = "free_x", labeller = "label_parsed") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 7)) +
    labs(x = "parameter value", y = "statistic")

## ---- cache = TRUE-------------------------------------------------------
## Confidence intervals based on the location-adjusted Wald statistic
rs_cor_ml_cis <- waldi_confint(rs_beta_ml, level = 0.95, adjust = TRUE, parallel = FALSE)
rs_cor_br_cis <- waldi_confint(rs_beta_br, level = 0.95, adjust = TRUE, parallel = FALSE)
## Studentized bootstrap intervals
set.seed(123)
quantiles_ml <- dyslexia_bootstrap(rs_beta_ml, R = 500, ncores = 1)$quantiles
quantiles_br <- dyslexia_bootstrap(rs_beta_br, R = 500, ncores = 1)$quantiles
rs_ml_stud_cis <- waldi_confint(rs_beta_ml, adjust = FALSE, parallel = FALSE,
                                quantiles = quantiles_ml$zstat[, c("0.025", "0.975")])
rs_br_stud_cis <- waldi_confint(rs_beta_br, adjust = FALSE, parallel = FALSE,
                                quantiles = quantiles_br$zstat[, c("0.025", "0.975")])
rs_cor_ml_stud_cis <- waldi_confint(rs_beta_ml, adjust = TRUE, parallel = FALSE,
                      quantiles = quantiles_ml$zstat_cor[, c("0.025", "0.975")])
rs_cor_br_stud_cis <- waldi_confint(rs_beta_br, adjust = TRUE, parallel = FALSE,
                      quantiles = quantiles_br$zstat_cor[, c("0.025", "0.975")])

round(rbind(cbind(rs_cor_ml_cis, rs_cor_br_cis),
            cbind(rs_cor_ml_stud_cis, rs_cor_br_stud_cis)), 3)
rs_coverage %>% filter(statistic %in% c("ml_cor", "br_cor")) %>%
    select(statistic, parameter, coverage.90, coverage.95, coverage.99)
rs_coverage %>% filter(statistic %in% c("ml_cor_stud", "br_cor_stud")) %>%
    select(statistic, parameter, coverage.90, coverage.95, coverage.99)

## ---- cache = TRUE-------------------------------------------------------
## Intervals based on the location-adjusted Wald statistic
system.time({
    waldi_confint(rs_beta_ml, adjust = TRUE, parallel = FALSE, length = 5)
})
simu_fun <- get_simulate_function(rs_beta_ml)
generate_dyslexia <- function(data, mle) {
    simu_fun(mle)
}
stat <- function(data, psi) {
    temp <- ReadingSkills
    temp$accuracy <- data
    temp_fit <- try(update(rs_beta_ml, data = temp))
    if (inherits(temp_fit, "try-error")) {
        rep(NA, 7)
    }
    else {
        waldi(temp_fit, null = psi, adjust = TRUE)
    }
}
## Studentized bootstrap intervals
system.time({
    stats <- boot(ReadingSkills$accuracy, statistic = stat,
                  R = 500, sim = "parametric", ran.gen = generate_dyslexia,
                  mle = coef(rs_beta_ml), psi = coef(rs_beta_ml), ncpus = 1)$t
    quant <- t(apply(stats, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
    waldi_confint(rs_beta_br, adjust = TRUE, parallel = FALSE,
                  quantiles = quant)
})

## ---- warning = FALSE----------------------------------------------------
source(paste0(code_path, "/", "logodds_functions.R"))
## Distribution of the statistic against normal
settings <- expand.grid(m = c(8, 16, 32), theta0 = c(-2, -1, 0))
plot_data <- NULL
for (j in seq.int(nrow(settings))) {
    setting <- settings[j, ]
    z <- seq(-3, 3, length = 100)
    dat <- t(sapply(z, dist_function, n = setting$m, theta0 = setting$theta0))
    dd <- stack(as.data.frame(dat))
    dd$z <- z
    names(dd) <- c("prob", "method", "z")
    dd$theta0 <- setting$theta0
    dd$m <- setting$m
    plot_data <- rbind(plot_data, dd)
}
plot_data$theta0 <- paste0("theta[0] == ", plot_data$theta0)
plot_data$theta0 <- factor(plot_data$theta0, levels = unique(plot_data$theta0),
                           ordered = TRUE)
plot_data$m <- paste0("n == ", plot_data$m)
plot_data$m <- factor(plot_data$m, levels = unique(plot_data$m), ordered = TRUE)
plot_data$method <- factor(plot_data$method, levels = c("ml", "a_ml", "br", "a_br"),
                           ordered = TRUE)
plot_data$method <- recode(plot_data$method,
                           "ml" = "italic(t)",
                           "a_ml" = "italic(t)^{'*'}",
                           "br" = "tilde(italic(t))",
                           "a_br" = "tilde(italic(t))^{'*'}")
ggplot(plot_data) +
    geom_abline(aes(intercept = 0, slope = 0), col = "grey") +
    geom_line(aes(z, qnorm(prob) - z), alpha = 0.5) +
    facet_grid(method ~ theta0 + m, label = label_parsed) +
    theme_minimal() +
    labs(y = expression(paste(Phi^list(-1),(italic(G)(italic(z))) - italic(z))),
         x = expression(italic(z))) +
    theme(text=element_text(size = 11))

## ------------------------------------------------------------------------
probs <- seq(1e-08, 1 - 1e-08, length = 500)
df <- ddply(data.frame(m = c(8, 16, 32, 64, 128, 256)), ~ m, function(x) {
    m <- x$m
    cis <- compute_cis(m, level = 0.95)
    cc <- lapply(probs, function(pp) cover_ci_prop(n = m, p = pp, level = 0.95, cis = cis))
    do.call("rbind", cc)
})
df$m <- factor(paste("m ==", df$m), levels = paste("m ==", sort(unique(df$m))),
               ordered = TRUE)
df$method <- factor(df$method, levels = c("wald", "ac", "a_br", "ml", "a_ml", "br"),
                    ordered = TRUE)
df$method <- recode(df$method,
                    "wald" = "Wald",
                    "ml" = "italic(t)[trans]",
                    "a_ml" = "italic(t)^*[trans]",
                    "br" = "tilde(italic(t))[trans]",
                    "a_br" = "tilde(italic(t))^{'*'}",
                    "ac" = "Agresti-Coull")
## coverage
ggplot(df %>% filter(method %in% c("Wald", "Agresti-Coull", "tilde(italic(t))^{'*'}"))) +
    geom_hline(aes(yintercept = 0.95), col = "grey") +
    geom_line(aes(x = p, y = coverage)) +
    facet_grid(m ~ method, label = label_parsed) +
    coord_cartesian(ylim = c(0.7, 1)) +
    theme_minimal()
## expected length
ggplot(df %>% filter(method %in% c("Wald", "Agresti-Coull", "tilde(italic(t))^{'*'}"))) +
    geom_line(aes(x = p, y = length, col = method), alpha = 0.5, size = 0.8) +
    facet_wrap(~ m, label = label_parsed, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = c("#328900", "#C54E6D", "#0080C5"),
                        name = "",
                        labels = c("Wald",
                                   "Agresti-Coull",
                                   bquote(tilde(italic(t))^{'*'}))) +
    theme_minimal() +
    theme(legend.position = "top")

## ------------------------------------------------------------------------
sapply(28:32, t_ml, n =  32, theta0 = 0)
sapply(28:32, t_adjusted_ml, n =  32, theta0 = 0)
sapply(28:32, t_br, n =  32, theta0 = 0)
sapply(28:32, t_adjusted_br, n =  32, theta0 = 0)

## ------------------------------------------------------------------------
## The clotting data set
clotting <- data.frame(
  conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
  u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
  lot = factor(c(rep(1, 9), rep(2, 9))))
## The maximum likelihood fit of the gamma regression model
clotting_ml <- glm(conc ~ log(u)*lot, data = clotting, family = Gamma(link = "log"))
## Maximum likelihood estimates and Wald statistics using maximum likelihood estimator
## of the dispersion parameter
dispersion_ml <- MASS::gamma.dispersion(clotting_ml)
clotting_summary_ml <- summary(clotting_ml, dispersion = dispersion_ml)
clotting_ml_estimates <- coef(clotting_summary_ml)[, c("Estimate", "z value")]
## Reduced-bias estimates and Wald statistics
clotting_summary_rb <- summary(update(clotting_ml, method = "brglmFit"))
## Maximum likelihood estimates and Wald statistics using the moment-based estimator
## of the dispersion parameter
clotting_summary_mom <-  summary(clotting_ml)
dispersion_mom <- clotting_summary_mom$dispersion
clotting_mom_estimates <- coef(clotting_summary_mom)[, c("Estimate", "t value")]
## Location-adjusted Wald statistic
clotting_waldi <- waldi(clotting_ml, null = 0, adjust = TRUE)
round(cbind(c(clotting_ml_estimates[, 1], dispersion_ml, dispersion_mom),
            c(clotting_ml_estimates[, 2], NA, NA),
            c(clotting_mom_estimates[, 2], NA, NA),
            c(clotting_waldi, NA, NA)), 3)

## ------------------------------------------------------------------------
load(paste(results_path, "clotting_simulation.rda", sep = "/"))

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

ggplot(typeI %>% filter(parameter != 1, is.element(statistic, 
                                              c("ml", "ml_cor", "mom",
                                                "rb", "rb_cor")))) +
    geom_point(aes(parameter, typeI, pch = statistic), alpha = 0.7) +
    geom_hline(aes(yintercept = level), col = "grey", lty = 2) +
    facet_grid(test ~ level_chr, labeller = label_parsed, scales = "free") +
    scale_x_continuous(name = element_blank(),
                       breaks = c(2, 3, 4),
                       limits = c(1.8, 4.2),
                       labels = c(
                           expression(beta[2]),
                           expression(beta[3]),
                           expression(beta[4]))) +
    scale_y_continuous(name = expression(paste("Empirical rejection probability (",
                                               symbol('\045'), ")")),
                       labels = function (x) {
                           if (length(x) == 0)
                               return(character())
                           x <- round_any(x, scales:::precision(x)/100)
                           scales:::comma(x * 100)
                       }) +
  scale_shape_manual(name = " ", values = c(16, 17, 15, 3, 4), 
                     labels = c(expression(italic(t)), expression(italic(t)^list("*")), 
            expression(paste("Wald with ", tilde(phi))), expression(tilde(italic(t))), 
            expression(tilde(italic(t))^list("*")))) +
    theme_bw() +
    theme(legend.position = "top",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_blank())

## ----fig_boot, echo = TRUE, fig.height = 4, fig.width = 7, fig.align = "center", message = FALSE, warning = FALSE----
ggplot(typeI %>% filter(parameter != 1, is.element(statistic, 
                                              c("ml", "ml_cor", "mom",
                                                "sc_ml_cor_conv")))) +
  geom_point(aes(parameter, typeI, pch = statistic), alpha = 0.7) +
  geom_hline(aes(yintercept = level), col = "grey", lty = 2) +
  facet_grid(test ~ level_chr, labeller = label_parsed, scales = "free") +
  scale_x_continuous(name = element_blank(),
                     breaks = c(2, 3, 4),
                     limits = c(1.8, 4.2),
                     labels = c(expression(beta[2]),
                       expression(beta[3]),
                       expression(beta[4]))) +
  scale_y_continuous(name = expression(paste("Empirical rejection probability 
                                             (", symbol('\045'), ")")),
                     labels = function (x) {
                       if (length(x) == 0)
                         return(character())
                       x <- round_any(x, scales:::precision(x)/100)
                       scales:::comma(x * 100)
                     }) +
  theme_bw() +
  scale_shape_manual(name = " ", values = c(16, 17, 15, 18), 
              labels = c(expression(italic(t)), expression(italic(t)^list("*")), 
            expression(paste("Wald with ", tilde(phi))),  
            expression(italic(t)^list("**")))) +
  theme(legend.position = "top",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank())

## ------------------------------------------------------------------------
data("babies", package = "cond")
## clogit understands only 0-1 so expand
babies_expand <- ddply(babies, ~ lull + day, function(z) {
    data.frame(y = rep(c(0, 1), c(z$r2, z$r1)))
})
## Maximum likelihood fit
babies_ml <- glm(formula = y ~ day + lull - 1,
                 family = binomial, data = babies_expand)
babies_rb <- update(babies_ml, method = "brglmFit")
## Maximum conditional likelihood fit
babies_cond <- clogit(y ~ strata(day) + lull, data = babies_expand)
ml <- coef(summary(babies_ml))["lullyes", ]
rb <- coef(summary(babies_rb))["lullyes", ]
mcl <- coef(summary(babies_cond))["lullyes", ]
r <- lrtest(update(babies_ml, . ~ . - lull),
            babies_ml)
rc <- summary(babies_cond)$logtest[1]
scorec <- summary(babies_cond)$sctest[1]
out1 <- c(
    ml = unname(ml["Estimate"]),
    rb = unname(rb["Estimate"]),
    mcl = unname(mcl["coef"]),
    wald_ml = unname(ml["z value"]),
    wald_mcl = unname(mcl["z"]),
    wald_rb = unname(rb["z value"]),
    r = unname(sign(ml["Estimate"]) * sqrt(r$Chisq[2])),
    rc = unname(sign(mcl["coef"]) * sqrt(rc)),
    wald_ml_adjusted = unname(waldi(babies_ml, which = 19)),
    wald_rb_adjusted = unname(waldi(babies_rb, which = 19)))
out2 <- c(
    ml_se = unname(ml["Std. Error"]),
    rb_se = unname(rb["Std. Error"]),
    mcl_se = unname(mcl["se(coef)"]),
    ml_p = ml["Pr(>|z|)"],
    mcl_p = mcl["Pr(>|z|)"],
    rb_p = rb["Pr(>|z|)"],
    r_p = 2 * pnorm(-abs(out1["r"])),
    rc_p = 2 * pnorm(-abs(out1["rc"])),
    cor_ml_p = 2 * pnorm(-abs(out1["wald_ml_adjusted"])),
    cor_rb_p = 2 * pnorm(-abs(out1["wald_rb_adjusted"])))
round(matrix(c(out1, out2), ncol = 10, byrow = TRUE,
             dimnames = list(NULL,
                             c("mle", "rb", "mcle", "wald_ml", "wald_mlc",
                               "wald_rb", "r", "rc", "wald_ml_adjusted",
                               "wald_rb_adjusted"))), 4)

## ---- warning = FALSE----------------------------------------------------
load(paste(results_path, "babies_simulation.rda", sep = "/"))

## The bootstrap p-value for the babies data is
set.seed(123)
babies_bootstrap(babies_ml, R = 1000)$conv
## Compute pvalues from the various statistics account for the existence of bootstrap
## p-values
pval <- ddply(res %>% filter(!infinite & !is.na(value) & type != "summary"),
              ~ name,
              function(data) {
    if (all(data$type == "bootstrap_statistic")) {
        data.frame(sample = pnorm(data$value),
                   test = gsub("boot_prep_|boot_conv_", "", data$name))
    }
    else {
        p2 <- 2 * pnorm(-abs(data$value))
        pl <- pnorm(data$value)
        pr <- 1 - pl
        data.frame(sample = c(p2, pl, pr),
                   test = rep(c("2sided", "left", "right"), each = length(p2))) }
})
## Get rid of left right 2sided from statistic names
pval <- pval %>% mutate(name = gsub("_left|_right|_2sided", "", name))
pval <- pval %>%
    filter(!(name %in% c("scorec", "boot_prep")) & test != "right") %>%
    mutate(test = dplyr::recode(test,
                                "2sided" = "gamma != 0",
                                "left" = "gamma < 0",
                                "right" = "gamma > 0"),
           name = factor(name,
                         levels = c("mle", "rbe", "r", "cond", "scorec", "rc",
                                    "boot_conv", "cor", "cor_rb"),
                         ordered = TRUE)) %>%
    mutate(name = factor(name,
                         levels = c("mle", "r", "boot_conv", "rbe",
                                    "cond", "scorec", "rc",
                                    "cor", "cor_rb"),
                         ordered = TRUE)) %>%
    mutate(statistic = dplyr::recode(name,
                                     "mle" = "italic(t)",
                                     "rbe" = "italic(tilde(t))",
                                     "r" = "italic(r)",
                                     "cond" = "italic(t)[c]",
                                     "scorec" = "italic(s)[c]",
                                     "rc" = "italic(r)[c]",
                                     "cor" = "italic(t)^'*'",
                                     "cor_rb" = "tilde(italic(t))^'*'",
                                     "boot_conv" = "italic(boot)"))
## Bin sample
breaks <- (0:20)/20
pval <- pval %>%
    group_by(statistic, test) %>%
    mutate(sample = cut(sample, breaks = breaks, include.lowest = TRUE)) %>%
    group_by(statistic, test, sample)
ggplot(pval) +
    geom_hline(aes(yintercept = 1)) +
    geom_bar(aes(x = sample, y = ..count../2500), fill = "darkgray", alpha = 0.5) +
    facet_grid(test ~ statistic, labeller = label_parsed) +
    theme_bw() +
    scale_x_discrete(breaks = c("[0,0.05]", "(0.25,0.3]", "(0.5,0.55]",
                                "(0.75,0.8]", "(0.95,1]"),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
    theme(legend.position = "top",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(x = expression(paste(italic(p), "-value ")), y = "Density")

## ------------------------------------------------------------------------
source(paste0(code_path, "/", "overlay2_nifti.R"))
load(paste(results_path, "brains_case_study.rda", sep = "/"))

## Check how many times LR failed, excluding trivial voxels, and compute probability of
## infinite estimates
fits_mat %>% filter(statistic == "r" & voxel != 1) %>% group_by(parameter) %>%
    summarize(failed = 100 * sum((value == -Inf) * count) / sum(count),
              infinite = 100 * sum(infinite * count) / sum(count))
## detections
fits_mat %>%
    group_by(parameter, statistic) %>%
    filter(statistic %in% c("z_br", "corz_br")) %>%
    summarize(detections = mean(value < -1 | value > 1) * 100)
## Empirical lesion counts
lesion_counts <- colSums(lesions)
lesion_counts[lesion_counts == 0] <- NA
nifti_counts <- nifti(img = array(lesion_counts, dim(white_matter)))
lumin <- c(45, 100)
cols_counts <-heat_hcl(n = max(lesion_counts, na.rm = TRUE),
                       h = c(265, 200),
                       c = c(80, 50),
                       l = lumin,
                       power = c(0.7, 2))
overlay2.nifti(white_matter, y = nifti_counts, z = 32,
               plot.type = "single", plane = "sagittal",
               col.y = cols_counts, title = "lesion counts per voxel",
               col.main = "white")
## Significance maps
param <- "DD"
low <- 1
upp <- 4
lumin <- c(25, 120)
cols <- c(heat_hcl(n = 32,
                   h = c(265, 200),
                   c = c(80, 50),
                   l = lumin,
                   power = c(0.7, 2)),
          rev(heat_hcl(n = 32,
                       h = c(10, 40),
                       c = c(80, 50),
                       l = lumin,
                       power = c(0.4, 1.3))))
for (stat in c("z_br", "corz_br")) {
    zz <- (fits_mat %>% filter(statistic == stat & parameter == param))
    zz <- zz$value[array_indices]
    ## Threshold as in Ge et al (2014, AOAS)
    low_ind <- abs(zz) < low
    low_ind[is.na(low_ind)] <- FALSE
    zz[low_ind] <- NA
    upp_ind <- abs(zz) >= upp
    upp_ind[is.na(upp_ind)] <- FALSE
    zz[upp_ind] <- sign(zz[upp_ind]) * upp
    nifti_z <- nifti(img = array(zz, dim(white_matter)))
    nifti_z[1,1,1] <- -upp
    nifti_z[1,1,2] <- upp
    main <- switch(stat,
                   z_br = expression(tilde(italic(t))),
                   corz_br = expression(tilde(italic(t))^'*'))
    overlay2.nifti(white_matter, y = nifti_z, z = 32, plot.type = "single",
                   plane = "sagittal", col.y = cols, title = main,
                   col.main = "white")
}
### Plot z_br vs corz_br per parameter
v1 <- fits_mat %>%
    filter(statistic == "z_br", parameter == param) %>%
    select(z_br_value = value, voxel, parameter)
v2 <- fits_mat %>%
    filter(statistic == "corz_br", parameter == param) %>%
    select(corz_br_value = value, voxel, parameter)
v <- join(v1, v2, by = c("voxel", "parameter"))
ggplot(v) +
    geom_point(aes(x = z_br_value, y = corz_br_value), alpha = 0.1, size = 0.5) +
    geom_abline(aes(intercept = 0, slope = 1), col = "grey") +
    coord_cartesian(xlim = c(-2.4, 2.4), ylim = c(-2.4, 2.4)) +
    theme_minimal() +
    labs(x = expression(tilde(italic(t))), y = expression(tilde(italic(t))^'*'))

## ------------------------------------------------------------------------
numerical_time <- system.time(
    numerical <- waldi(babies_ml, numerical = TRUE, which = 19)
)
analytic_time <- system.time(
    analytic <- waldi(babies_ml, numerical = FALSE, which = 19)
)
(numerical_time/analytic_time)["elapsed"]

## ------------------------------------------------------------------------
load(paste(results_path, "babies_times.rda", sep = "/"))

times_ana <- sapply(out_ana, function(x) rowMeans(x)["elapsed"])
times_num <- sapply(out_num, function(x) rowMeans(x)["elapsed"])
cores <- seq_along(out_num)
res <- data.frame(times = c(times_ana, times_num), 
		   cores = c(cores, cores), 
		   derivatives = rep(c("analytical", "numerical"), each = length(cores)))

ggplot(res) +  
   geom_line(aes(x = cores, y = times, lty = derivatives)) +
   geom_point(aes(x = cores, y = times)) +
   theme_minimal() +
   labs(x = "cores", y = "elapsed time (sec)")

## ----metareg, echo = TRUE, message = FALSE, warning = FALSE--------------
load(paste(results_path, "meta_analysis_simulation.rda", sep = "/"))

## plot K = 10, 20
fig1 <- ggplot(cov_df_K  %>% filter(test == "2sided")) +
  geom_line(aes(psi, cov, group = statistic, col = statistic), size = 0.5, alpha = 0.8) +
  geom_hline(aes(yintercept = 95), col = "grey") +
  facet_wrap( ~ Klab, labeller = label_parsed, scales = "free") +
  labs(y = "Coverage probability (%)", x = expression(psi)) +
  lims(y = c(85, 100)) +
  scale_x_continuous(name = expression(psi), breaks = seq(0, 0.1, length = 6), 
                     labels = c("0.00", "0.02", "0.04", "0.06", "0.08", "0.10")) +
  scale_colour_manual(name = "", values = c("#328900", "#0080C5", "#C54E6D", "purple", 
        "orange", 4), labels = c(expression(italic(t)), expression(italic(t)^list("*")), 
        expression(tilde(italic(t))), 
        expression(tilde(italic(t))^list("*")), "DL", "ZL")) +
  theme_bw() +
  theme(legend.position = "top", panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), 
        strip.background = element_blank()) +
  guides(colour = guide_legend(nrow = 1))

## plot psi = 0.03, 0.07
fig2 <- ggplot(cov_df_psi  %>% filter(test == "2sided")) +
  geom_line(aes(log(K), cov, group = statistic, col = statistic), size = 0.5, 
            alpha = 0.8) +
  geom_hline(aes(yintercept = 95), col = "grey") +
  facet_wrap( ~ psilab, labeller = label_parsed, scales = "free") +
  labs(y = "Coverage probability (%)", x = expression(italic(K))) +
  lims(y = c(80, 100)) +
  scale_x_continuous(name = expression(italic(K)), breaks = c(log(5), log(10), 
                    log(25), log(50), log(100), log(200)), labels = 
                      c("5", "10", "25", "50", "100", "200")) +
  scale_colour_manual(name = "", values = c("#328900", "#0080C5", "#C54E6D", 
                      "purple", "orange", 4), labels = c(expression(italic(t)), 
                      expression(italic(t)^list("*")), expression(tilde(italic(t))), 
                      expression(tilde(italic(t))^list("*")), "DL", "ZL")) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), 
        strip.background = element_blank())


## ----figs, echo = FALSE, fig.height = 4, fig.width = 7, fig.align = "center", warning = FALSE----
print(fig1, width = 10)
print(fig2)

