library("plyr")
library("foreach")
library("doMC")
registerDoMC(4)

data("babies", package = "corzed")

babies_ml0 <- glm(formula = y ~ day + lull - 1,
                  family = binomial, data = babies)


ci95_default <- confint.default(babies_ml0, level = 0.95, parallel = TRUE)
ci95_corzed <- corzed_confint(babies_ml0, level = 0.95, adjust = FALSE, parallel = TRUE)
ci99_default <- confint.default(babies_ml0, level = 0.99, parallel = TRUE)
ci99_corzed <- corzed_confint(babies_ml0, level = 0.99, adjust = FALSE, parallel = TRUE)

tolerance <- 1e-04

test_that("corzed_ci returns correct CIs when adjust = FALSE", {
    expect_equal(ci95_default, ci95_corzed, tolerance = tolerance)
    expect_equal(ci99_default, ci99_corzed, tolerance = tolerance)
})


babies_ml1 <- glm(formula = y ~ day + lull + I(as.integer(lull) + 2) + I(as.integer(lull) + 3) - 1,
                  family = binomial, data = babies)

ci95_corzed1 <- corzed_confint(babies_ml1, level = 0.95, adjust = TRUE, which = c(16, 18, 20), parallel = TRUE)
ci95_corzed2 <- corzed_confint(babies_ml1, level = 0.95, adjust = TRUE, which = names(coef(babies_ml1))[c(16, 18, 20)], parallel = TRUE)

test_that("which works with both names and indexes in corzed_confint", {
    expect_equal(ci95_corzed1, ci95_corzed1, tolerance = tolerance)
})

ci95_corzed2 <- corzed_confint(babies_ml1, level = 0.95, adjust = TRUE, which = c(20, 21))
test_that("corzed_confint returns NA for unidentified parameters", {
    expect_true(all(is.na(ci95_corzed2)))
})
