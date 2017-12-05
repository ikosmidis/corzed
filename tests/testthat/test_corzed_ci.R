library("plyr")

data("babies", package = "cond")

## clogit understands only 0-1 so expand
babies_expand <- ddply(babies, ~ lull + day, function(z) {
    data.frame(y = rep(c(0, 1), c(z$r2, z$r1)))
})

babies_ml0 <- glm(formula = y ~ day + lull - 1,
                  family = binomial, data = babies_expand)


ci95_default <- confint.default(babies_ml0, level = 0.95)
ci95_corzed <- corzed_confint(babies_ml0, level = 0.95, adjust = FALSE)
ci99_default <- confint.default(babies_ml0, level = 0.99)
ci99_corzed <- corzed_confint(babies_ml0, level = 0.99, adjust = FALSE)

tolerance <- 1e-04

test_that("corzed_ci returns correct CIs when adjust = FALSE", {
    expect_equal(ci95_default, ci95_corzed, tolerance = tolerance)
    expect_equal(ci99_default, ci99_corzed, tolerance = tolerance)
})


babies_ml1 <- glm(formula = y ~ day + lull + I(as.integer(lull) + 2) + I(as.integer(lull) + 3) - 1,
                  family = binomial, data = babies_expand)

ci95_corzed1 <- corzed_confint(babies_ml1, level = 0.95, adjust = TRUE, which = c(1, 4, 5, 6, 11, 16, 18, 20))
ci95_corzed2 <- corzed_confint(babies_ml1, level = 0.95, adjust = TRUE, which = names(coef(babies_ml1))[c(1, 4, 5, 6, 11, 16, 18, 20)])

test_that("which works with both names and indexes in corzed_confint", {
    expect_equal(ci95_corzed1, ci95_corzed1, tolerance = tolerance)
})

ci95_corzed2 <- corzed_confint(babies_ml1, level = 0.95, adjust = TRUE, which = c(20, 21))
test_that("corzed_confint returns NA for unidentified parameters", {
    expect_true(all(is.na(ci95_corzed2)))
})
