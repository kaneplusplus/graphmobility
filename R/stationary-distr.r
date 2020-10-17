#' Find the Stationary Distribution of a Markov Transition Matrix
#' @param P a Markov transition matrix.
#' @param .eigen the function to calculate the eigen vectors with. By 
#' default, when the matrix is small the full eigenvalue decomposition is 
#' peformed. Whe the matrix is large, only the leading 20 eigen vectors are 
#' calculated.
#' @importFrom RSpectra eigs
#' @importFrom purrr partial
stationary_distr <- function(P, 
  .eigen = ifelse(nrow(P) > 50, 
                  partial(eigs, k = 20),
                  eigen)) {

  v <- .eigen(t(P), symmetric = FALSE)$vectors
  vn <- colSums(v)

  v <- sweep(v, 2, vn, FUN = "/")
  im <- apply(v, 2, function(x) sum(abs(Im(x))))
  v <- v[, abs(im < 1e-7), drop = FALSE]
  if (ncol(v) == 0) {
    warning("No real solution found. Returning zeros")
    matrix(0, nrow = nrow(v), ncol = 1)
  } else {
    v <- Re(v)
    good_vals <- apply(v, 2, function(x) all(x > 0))
    if (!any(good_vals)) {
      warning("No non-zero solution found. Returning zeros")
      matrix(0, nrow = nrow(v), ncol = 1)
    } else {
      v[, which(good_vals)[1], drop = FALSE]
    }
  }
}

