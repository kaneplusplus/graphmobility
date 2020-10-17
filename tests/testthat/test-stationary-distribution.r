context("Stationary distribution.")

P <- matrix(c(0, 1, 0, 1/3, 0, 2/3, 1/3, 1/3, 1/3), nrow = 3, byrow = TRUE)

expect_equal(stationary_distr(P), 
             matrix(c(1/4, 3/8, 3/8), ncol = 1))


# It would be nice to add a larger example.
