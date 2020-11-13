
#' @importFrom Matrix spMatrix diag 
make_diag <- function(mg) {
  spMatrix(nrow = nrow(mg), 
           ncol = ncol(mg), 
           i = seq_len(nrow(mg)),
           j = seq_len(ncol(mg)),
           x = Matrix::diag(mg))
}

#' Mobility Proportion 
#' @param mg the mobility graph.
#' @export
mobility_proportion <- function(mg) {
  sum(mg - make_diag(mg)) / sum(mg)
}

#' Stationary Proportion 
#' @param mg the mobility graph.
#' @importFrom Matrix diag
#' @export
stationary_proportion <- function(mg) {
  1 - mobility_proportion(mg)
}

#' Mobility Asymmetry
#' @param mg the mobility graph.
#' @param prop should the mobility asymmetry be reported as a proportion? 
#' Default TRUE.
#' @importFrom Matrix t
#' @export
mobility_asymmetry <- function(mg, prop = TRUE) {
  hg <- mg - make_diag(mg)
  ret <- sum(abs(hg - Matrix::t(hg))) / 2
  if (prop) {
    ret <- ret / sum(abs(mg))
  } 
  ret
}


