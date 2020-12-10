#' Find the Stationary Distribution of a Markov Transition Matrix
#' @details Note that this 
#' @param P a Markov transition matrix.
#' @param .eigen the function to calculate the eigen vectors with. By 
#' default, when the matrix is small the full eigenvalue decomposition is 
#' peformed. Whe the matrix is large, only the leading 20 eigen vectors are 
#' calculated.
#' @param tol the numerical tolerance. Default 1e-12.
#' @importFrom RSpectra eigs
#' @importFrom Matrix rowSums diag crossprod tcrossprod sparseMatrix
#' @importFrom purrr partial
#' @importFrom stats optimize optim
#' @export
stationary_distr <- function(P, make_hollow = FALSE) {

  if (make_hollow) {
    I <- sparseMatrix(i = seq_len(nrow(P)),
                      j = seq_len(ncol(P)),
                      x = Matrix(diag(P)))
    P <- P - I
  }
  if (length(find_ptm_sinks(P)) > 0) {
    stop("Probability transition matrix has sinks.")
  }

  components(
    graph_from_adjacency_matrix(P * 100, 
                                mode = "directed",
                                add.colnames = TRUE))
  as.numeric(
    Matrix::solve(crossprod(P), 
                  tcrossprod(P,  rep(1/nrow(P), nrow(P)))))
}

#' @importFrom Matrix rowSums
#' @export
ptm_sink_rows <- function(P) {
  which(Matrix::rowSums(P) == 1)
}

#' Mobility Graph to Probability Transition Matrix
#' @param mg the mobility graph.
#' @importFrom foreach foreach %dopar% getDoParName registerDoSEQ
#' @importFrom itertools isplitRows
#' @importFrom Matrix rowSums diag
#' @export
mg_to_ptm <- function(mg) {

  if (is.null(getDoParName())) {
    registerDoSEQ()
  }

  if (nrow(mg) > 100) {
    num_workers <- getDoParWorkers()
  } else{
    num_workers <- 1
  }

  ret <- 
    foreach(mgs = isplitRows(mg, chunks = num_workers), 
            .combine = rbind) %dopar% {
    ss <- Matrix::rowSums(mgs)
    foreach (i = seq_along(ss), .combine = rbind) %do% {  
      if (ss[i] == 0) {
        mgs[i,, drop = FALSE]
      } else {
        mgs[i, ] <- mgs[i, , drop = FALSE] / ss[i]
        mgs[i,, drop = FALSE]
      }
    }
  }
  #browser()
  #ones <- which(Matrix::diag(ret) == Matrix::rowSums(ret))
  ret
}
