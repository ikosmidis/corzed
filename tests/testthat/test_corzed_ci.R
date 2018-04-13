library("plyr")
library("foreach")
library("doMC")
registerDoMC(4)

data("babies", package = "waldi")

babies_ml0 <- glm(formula = y ~ day + lull - 1,
                  family = binomial, data = babies)


ci95_default <- confint.default(babies_ml0, level = 0.95, parallel = TRUE)
ci95_waldi <- waldi_confint(babies_ml0, level = 0.95, adjust = FALSE, parallel = TRUE)
ci99_default <- confint.default(babies_ml0, level = 0.99, parallel = TRUE)
ci99_waldi <- waldi_confint(babies_ml0, level = 0.99, adjust = FALSE, parallel = TRUE)

tolerance <- 1e-04

test_that("waldi_ci returns correct CIs when adjust = FALSE", {
    expect_equal(ci95_default, ci95_waldi, tolerance = tolerance)
    expect_equal(ci99_default, ci99_waldi, tolerance = tolerance)
})


babies_ml1 <- glm(formula = y ~ day + lull + I(as.integer(lull) + 2) + I(as.integer(lull) + 3) - 1,
                  family = binomial, data = babies)

ci95_waldi1 <- waldi_confint(babies_ml1, level = 0.95, adjust = TRUE, which = c(16, 18, 20), parallel = TRUE)
ci95_waldi2 <- waldi_confint(babies_ml1, level = 0.95, adjust = TRUE, which = names(coef(babies_ml1))[c(16, 18, 20)], parallel = TRUE)

test_that("which works with both names and indexes in waldi_confint", {
    expect_equal(ci95_waldi1, ci95_waldi1, tolerance = tolerance)
})

ci95_waldi2 <- waldi_confint(babies_ml1, level = 0.95, adjust = TRUE, which = c(20, 21))
test_that("waldi_confint returns NA for unidentified parameters", {
    expect_true(all(is.na(ci95_waldi2)))
})
