library("plyr")

data("babies", package = "cond")

## clogit understands only 0-1 so expand
babies_expand <- ddply(babies, ~ lull + day, function(z) {
    data.frame(y = rep(c(0, 1), c(z$r2, z$r1)))
})

babies_ml0 <- glm(formula = y ~ day + lull - 1,
                  family = binomial, data = babies_expand)
babies_ml1 <- glm(formula = y ~ day + lull + I(as.integer(lull) + 2) + I(as.integer(lull) + 3) - 1,
                  family = binomial, data = babies_expand)

c0 <- corzed(babies_ml0)
c1 <- corzed(babies_ml1)

tolerance <- 1e-06
test_that("corzed.glm handles unidentifiable parameters correctly", {
    expect_equal(c1[names(c0)], c0, tolerance = tolerance)
})


c2 <- corzed(babies_ml1, which = c("I(as.integer(lull) + 3)", "day9", "lullyes"))
test_that("corzed.glm handles unidentifiable parameters in which correctly", {
    expect_equal(c1[names(c2)], c2, tolerance = tolerance)
})


c3 <- corzed(babies_ml1, which = c(1, 4, 5, 6, 11, 16, 18, 20))
test_that("corzed.glm handles unidentifiable parameters in which correctly", {
    expect_equal(c1[names(c3)], c3, tolerance = tolerance)
})
