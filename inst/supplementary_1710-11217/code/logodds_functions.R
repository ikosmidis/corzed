Tders <- function(theta, n, theta0) {
    T0_expr <- expression((theta - theta0) * sqrt(n * exp(theta) / (1 + exp(theta))^2))
    a <- eval(deriv3(T0_expr, "theta"))
    list(T0 = a[1],
         T1 = unname(drop(attr(a, "gradient"))),
         T2 = unname(drop(attr(a, "hessian"))))
}

info <- function(theta, n, theta0) {
    prob <- plogis(theta)
    n * prob * (1 - prob)
}

bias <- function(theta, n, theta0) {
    prob <- plogis(theta)
    -0.5  * (1 - 2 * prob) / info(theta, n, theta0)
}

t_ml <- function(y, n, theta0) {
    theta <- log(y/(n - y))
    if (y == 0 | y == n) {
        0
    }
    else {
        (theta - theta0) * sqrt(info(theta, n, theta0))
    }
}

t_br <- function(y, n, theta0) {
    theta <- log((y + 0.5)/(n - y + 0.5))
    (theta - theta0) * sqrt(info(theta, n, theta0))
}

t_adjusted_ml <- function(y, n, theta0) {
    theta <- log(y/(n - y))
    tders <- Tders(theta, n, theta0)
    b <- bias(theta, n, theta0)
    i_inv <- 1/info(theta, n, theta0)
    if (y == 0 | y == n) {
        0
    }
    else {
        with(tders, T0 - T1 * b - 0.5 * T2 * i_inv)
    }
}

t_bias <- function(theta, n, theta0) {
    tders <- Tders(theta, n, theta0)
    b <- bias(theta, n, theta0)
    i_inv <- 1/info(theta, n, theta0)
    with(tders, T1 * b + 0.5 * T2 * i_inv)
}

t_adjusted_br <- function(y, n, theta0) {
    theta <- log((y + 0.5)/(n - y + 0.5))
    tders <- Tders(theta, n, theta0)
    b <- bias(theta, n, theta0)
    i_inv <- 1/info(theta, n, theta0)
    with(tders, T0 - 0.5 * T2 * i_inv)
}

## pmfs
dist_function <- function(z, n, theta0, pvalue = FALSE) {
    y <- 0:n
    prob <- dbinom(y, n, exp(theta0)/(1 + exp(theta0)))
    ml <- sapply(y, t_ml, n = n, theta0 = theta0)
    br <- sapply(y, t_br, n = n, theta0 = theta0)
    a_ml <- sapply(y, t_adjusted_ml, n = n, theta0 = theta0)
    a_br <- sapply(y, t_adjusted_br, n = n, theta0 = theta0)
    if (pvalue) {
        out <- c(ml = sum(prob[2 * pnorm(-abs(ml)) <= z]),
                 br = sum(prob[2 * pnorm(-abs(br)) <= z]),
                 a_ml = sum(prob[2 * pnorm(-abs(a_ml)) <= z]),
                 a_br = sum(prob[2 * pnorm(-abs(a_br)) <= z]))
    }
    else {
        out <- c(ml = sum(prob[ml <= z]),
                 br = sum(prob[br <= z]),
                 a_ml = sum(prob[a_ml <= z]),
                 a_br = sum(prob[a_br <= z]))
        out[out == 0 | out == 1] <- NA
    }
    out
}

## coverage
cover <- function(n, theta0, level = 0.95) {
    y <- 0:n
    prob <- dbinom(y, n, exp(theta0)/(1 + exp(theta0)))
    ml <- sapply(y, t_ml, n = n, theta0 = theta0)
    br <- sapply(y, t_br, n = n, theta0 = theta0)
    a_ml <- sapply(y, t_adjusted_ml, n = n, theta0 = theta0)
    a_br <- sapply(y, t_adjusted_br, n = n, theta0 = theta0)
    quant <- qnorm(1 - (1 - level)/2)
    ml <- (ml <= quant) & (ml >= -quant)
    br <- (br <= quant) & (br >= -quant)
    a_ml <- (a_ml <= quant) & (a_ml >= -quant)
    a_br <- (a_br <= quant) & (a_br >= -quant)
    c(ml = sum(prob[ml]),
      br = sum(prob[br]),
      a_ml = sum(prob[a_ml]),
      a_br = sum(prob[a_br]))
}

ci <- function(y, n, level = 0.95, statistic = t_ml, interval = c(-10, 10)) {
    quant <- qnorm(1 - (1 - level)/2)
    low <- function(theta) {
        (statistic(y, n, theta) - quant)^2
    }
    upp <- function(theta) {
        (statistic(y, n, theta) + quant)^2
    }
    c(low = optimize(low, interval)$minimum,
      upper = optimize(upp, interval)$minimum)
}


ci_ml_prop <- function(y, n, level = 0.95) {
    quant <- qnorm(1 - (1 - level)/2)
    p <- y/n;
    c(low = p - quant * sqrt(p * (1- p) / n), upper = p + quant * sqrt(p * (1- p) / n))
}


ci_ac_prop <- function(y, n, level = 0.95) {
    quant <- qnorm(1 - (1 - level)/2)
    p <- (y + 2)/(n + 4);
    c(low = p - quant * sqrt(p * (1- p) / (n + 4)), upper = p + quant * sqrt(p * (1- p) / (n + 4)))
}

compute_cis <- function(n, level = 0.95) {
    y <- 0:n
    ci_wald <- data.frame(t(sapply(y, ci_ml_prop , n = n, level = level)),
                          method = "wald")
    ## Get CI by plogis of the endpoints
    ci_ml <- data.frame(t(plogis(sapply(y, ci, n = n, level = level, statistic = t_ml))),
                        method = "ml")
    ci_a_ml <- data.frame(t(plogis(sapply(y, ci, n = n, level = level, statistic = t_adjusted_ml))),
                          method = "a_ml")
    ci_ml[y == 0, c("low", "upper")] <-
        ci_ml[y == n, c("low", "upper")] <-
        ci_a_ml[y == 0, c("low", "upper")] <-
        ci_a_ml[y == n, c("low", "upper")] <- c(0, 1)
    ci_br <- data.frame(t(plogis(sapply(y, ci, n = n, level = level, statistic = t_br))),
                        method = "br")
    ci_a_br <- data.frame(t(plogis(sapply(y, ci, n = n, level = level, statistic = t_adjusted_br))),
                          method = "a_br")
    ci_ac_prop <- data.frame(t(sapply(y, ci_ac_prop , n = n, level = level)),
                             method = "ac")
    cis <- rbind(ci_wald, ci_ml, ci_a_ml, ci_br, ci_a_br, ci_ac_prop)
    cis$y <- y
    cis$n <- n
    cis
}

cover_ci_prop <- function(n, p, level = 0.95, cis = NULL) {
    if (is.null(cis)) {
        cis <- compute_cis(n, level = level)
    }
    cis$prob <- dbinom(0:n, n, p)
    cis$cover <- cis[, "low"] < p & cis[, "upper"] > p
    ddply(cis, ~ method, function(x) {
        with(x,
             c(coverage = sum(prob * cover),
               length = sum(prob * (x["upper"] - x["low"])),
               p = p))})
}


