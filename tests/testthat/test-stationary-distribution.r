context("Stationary distribution.")

# From http://www.stat.yale.edu/~pollard/Courses/251.spring2013/Handouts/Chang-MarkovChains.pdf
P <- t(matrix(c(0, 1, 0, 1/3, 0, 2/3, 1/3, 1/3, 1/3), nrow = 3, byrow = TRUE))

expect_equal(stationary_distr(P), 
             matrix(c(1/4, 3/8, 3/8), ncol = 1))


