### must be run on a machine with at least 28 available cores
library("doMC")
library("plyr")
data("babies", package = "cond")

babies_expand <- ddply(babies, ~ lull + day, function(z) {
  data.frame(y = rep(c(0, 1), c(z$r2, z$r1)))
})

babies_ml <- glm(formula = y ~ day + lull - 1,
                 family = binomial, data = babies_expand)

out_num <- NULL
out_ana <- NULL

for (j in 1:28) {
  registerDoMC(j)
  out_num[[j]] <- replicate(10, system.time(waldi(babies_ml, parallel = TRUE, numerical = TRUE)))
  out_ana[[j]] <- replicate(10, system.time(waldi(babies_ml, parallel = TRUE, numerical = FALSE))) 
  cat(j, "\n")
}

save(out_num, out_ana, file = "babies_times.rda")