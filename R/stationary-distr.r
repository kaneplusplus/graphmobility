#' Find the Stationary Distribution of a Markov Transition Matrix
#' @param P a Markov transition matrix. Rows correspond to the "from" 
#' vertex, columns correspond to the "to" vertex.
#' @importFrom Matrix sparseMatrix
#' @importFrom CVXR Variable Minimize Problem solve
#' @export
#' @importFrom purrr partial
#' @importFrom igraph is_directed
#' @importFrom checkmate assert check_class
stationary_distr <- function(P) {

  assert(
    check_class(P, "Matrix"),
    check_class(P, "matrix"))

  P <- t(P)

  I <- sparseMatrix(i = seq_len(nrow(P)),
                      j = seq_len(ncol(P)),
                      x = rep(1, ncol(P))) 
  B <- sparseMatrix(i = nrow(P) + 1, j = 1, x = c(rep(nrow(P), 0), 1))
  A <- P - I
  A <- rbind(
    A, 
    sparseMatrix(
      i = rep(1, nrow(P)), 
      j = seq_len(nrow(P)), 
      x = rep(1, nrow(P))))

  pi <- Variable(ncol(A))
  Ad <- as.matrix(A)
  Bd <- as.matrix(B)
  obj <- Minimize( sum((Ad %*% pi - Bd)^2) )
  prob <- Problem(obj, constraints = list(pi >= 0))
  ret <- solve(prob)
  ret$getValue(pi)
}

#' @importFrom checkmate assert check_class
prune_edges <- function(g, tresh = NULL) {

  assert(
    check_class(g, "igraph"),
    is_directed(g),
    combine = "and"
  )

  
}

#' @importFrom checkmate assert check_class
num_clusters <- function(cg) {
  assert(check_class(cg, "clusGap"))
  tab <- cg$Tab
  max_gap <- which.max(diff(tab[,"gap"])) + 1
  mt <- tab[max_gap, , drop = FALSE]
  tab <- tab[-max_gap, , drop = FALSE]
  uppers <- tab[,"gap"] + 2 * (tab[,"SE.sim"] + mt[,"SE.sim"])
  if (any(uppers >= mt[,"gap"])) {
    warning("Gap separation is less than 2 sd.",
            "Verify the number of clusters manually.")
  }
  max_gap
}



#' @importFrom Matrix rowSums
ptm_sink_rows <- function(P) {
  which(Matrix::rowSums(P) == 1)
}

#' Mobility Graph to Probability Transition Matrix
#' @param am the adjacency matrix.
#' @importFrom foreach foreach %dopar% getDoParName registerDoSEQ 
#' getDoParWorkers %do%
#' @importFrom itertools isplitRows
#' @importFrom Matrix rowSums diag
#' @aliases am_to_ptm
#' @export
mg_to_ptm <- function(am) {

  if (is.null(getDoParName())) {
    registerDoSEQ()
  }

  if (nrow(am) > 100) {
    num_workers <- getDoParWorkers()
  } else{
    num_workers <- 1
  }

  ret <- 
    foreach(ams = isplitRows(am, chunks = num_workers), 
            .combine = rbind) %dopar% {
    ss <- Matrix::rowSums(ams)
    foreach (i = seq_along(ss), .combine = rbind) %do% {  
      if (ss[i] == 0) {
        ams[i,, drop = FALSE]
      } else {
        ams[i, ] <- ams[i, , drop = FALSE] / ss[i]
        ams[i,, drop = FALSE]
      }
    }
  }
  dimnames(ret) <- dimnames(am)
  ret
}

#' @export
am_to_ptm <- mg_to_ptm

