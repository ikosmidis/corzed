library("plyr")
library("foreach")
library("doMC")
registerDoMC(4)

data("babies", package = "corzed")

babies_ml0 <- glm(formula = y ~ day + lull - 1,
                  family = binomial, data = babies)
babies_ml1 <- glm(formula = y ~ day + lull + I(as.integer(lull) + 2) + I(as.integer(lull) + 3) - 1,
                  family = binomial, data = babies)

c0 <- corzed(babies_ml0, parallel = TRUE)
c1 <- corzed(babies_ml1, parallel = TRUE)

tolerance <- 1e-06
test_that("corzed.glm handles unidentifiable parameters correctly", {
    expect_equal(c1[names(c0)], c0, tolerance = tolerance)
})


c2 <- corzed(babies_ml1, which = c("I(as.integer(lull) + 3)", "day9", "lullyes"), parallel = TRUE)
test_that("corzed.glm handles unidentifiable parameters in which correctly", {
    expect_equal(c1[names(c2)], c2, tolerance = tolerance)
})


c3 <- corzed(babies_ml1, which = c(1, 4, 5, 6, 11, 16, 18, 20), parallel = TRUE)
test_that("corzed.glm handles unidentifiable parameters in which correctly", {
    expect_equal(c1[names(c3)], c3, tolerance = tolerance)
})
