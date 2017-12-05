clotting <- data.frame(
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    lot = factor(c(rep(1, 9), rep(2, 9))))

modML <- glm(conc ~ log(u)*lot, data = clotting, family = Gamma(link="log"))

c_num <- corzed(modML, numeric = FALSE)
c_ana <- corzed(modML, numeric = TRUE)

tolerance <- 1e-06
test_that("analytic implementation gives the same result as numerical one", {
    expect_equal(c_num, c_ana, tolerance = tolerance)
})

modML1 <- glm(conc ~ log(u)*lot + I(2 * log(u)), data = clotting, family = Gamma(link="log"))
c1_num <- corzed(modML1, numeric = FALSE)
c1_ana <- corzed(modML1, numeric = TRUE)

test_that("analytic implementation gives the same result as numerical one with aliasing", {
    expect_equal(c1_num, c1_ana, tolerance = tolerance)
})

