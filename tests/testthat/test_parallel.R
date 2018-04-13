library("foreach")
library("doMC")
registerDoMC(4)

data("babies", package = "waldi")

babies_ml0 <- glm(formula = y ~ day + lull - 1,
                  family = binomial, data = babies)

c0 <- waldi(babies_ml0, parallel = FALSE)
c1 <- waldi(babies_ml0, parallel = TRUE)

tolerance <- 1e-06
test_that("use of parallel returns the same statistics as serial evaluation", {
    expect_equal(c1, c0, tolerance = tolerance)
})


ci0 <- waldi_confint(babies_ml0, parallel = FALSE, which = c(2, 5, 12, 19))
ci1 <- waldi_confint(babies_ml0, parallel = TRUE, which = c(2, 5, 12, 19))

test_that("use of parallel returns the same confidence intervals as serial evaluation", {
    expect_equal(ci1, ci0, tolerance = tolerance)
})
