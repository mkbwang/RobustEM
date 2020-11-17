library(testthat)
library(RobustEM)

sim_info <- simMultGauss(n = 120, d = 2, cluster = 6, out_perc = 0.06, out_mag = 5)
result <- robustEM(sim_info[["simdata"]], cluster = 6)

real_mu <- sim_info$mus
real_mu <- real_mu[order(real_mu[, 1]), ]

output_mu <- result$mu
output_mu <- output_mu[order(output_mu[, 1]), ]


test_that("Mean of each cluster matches", {
  expect_equal(sum(abs(real_mu-output_mu)), 0,
               tolerance = nrow(real_mu) * ncol(real_mu))
})
