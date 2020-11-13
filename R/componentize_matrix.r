
#' Break up an Adjacency Matrix by Components
#' @param x an adjacency matrix
#' @param comp_name the name of the output component column. 
#' Default "component".
#' @param adj_mat_name the name of the output adjacency matrix. Default
#' "adj_mat".
#' @param size_name the name of the column denoting the component size. Default 
#' "size".
#' @importFrom foreach foreach %dopar% getDoParName registerDoSEQ
#' @importFrom igraph components graph_from_adjacency_matrix
#' @importFrom tibble tibble
#' @importFrom itertools isplitVector
#' @importFrom dplyr bind_rows
#' @export
componentize_matrix <- 
  function(x, comp_name = "component", adj_mat_name = "adj_mat",
           size_name= "size") {

  if (is.null(mg)) {
    stop("No rownames or columns.")
  } else {
    rne <- dimnames(x)
    if (is.null(rne[[1]]) || is.null(rne[[2]]) || any(rne[[1]] != rne[[2]])) {
      stop("Dimension names not properly specified.")
    }
  }

  if (is.null(getDoParName())) {
    registerDoSEQ()
  }
  num_workers <- getDoParWorkers()

  singleton_b <- Matrix::diag(x) == Matrix::rowSums(x) & 
    Matrix::diag(x) == Matrix::colSums(x)
  singleton_inds <- which(singleton_b)
  
  ret_single <- 
    foreach(it = isplitVector(singleton_inds, chunks = getDoParWorkers()), 
            .combine = bind_rows) %dopar% {
      foreach(ind = it, .combine = bind_rows) %do% {
        tibble({{comp_name}} := ind, 
               {{adj_mat_name}} := list(x[ind, ind, drop = FALSE]),
               {{size_name}} := 1)
      }
  }
  ret_single[[comp_name]] <- seq_len(nrow(ret_single))

  x <- x[-singleton_inds, -singleton_inds, drop = FALSE]
  comps <- components(
    graph_from_adjacency_matrix(x, mode = "max", 
                                diag = TRUE, 
                                add.colnames = TRUE)
  )

  ret <- tibble({{comp_name}} := seq_len(comps$no))
  ret[[adj_mat_name]] <- 
    foreach(
      it = isplitVector(ret[[comp_name]], chunks = getDoParWorkers()),
      .combine = c) %dopar% {

      foreach(i = it) %do% {
        comp_inds <- which(i == comps$membership)
        x[comp_inds, comp_inds, drop = FALSE]
      }
    }
  ret[[size_name]] <-     
    foreach(
      it = isplitVector(ret[[comp_name]], chunks = getDoParWorkers()),
      .combine = c) %dopar% {

      foreach(i = it, .combine = c) %do% {
        nrow(ret[[adj_mat_name]][[i]])  
      }
    }
  rbind(ret[order(ret$size, decreasing = TRUE),], ret_single)
}
