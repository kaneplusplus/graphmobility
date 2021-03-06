% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/componentize_matrix.r
\name{componentize_matrix}
\alias{componentize_matrix}
\title{Break up an Adjacency Matrix by Components}
\usage{
componentize_matrix(
  x,
  comp_name = "component",
  adj_mat_name = "adj_mat",
  size_name = "size"
)
}
\arguments{
\item{x}{an adjacency matrix}

\item{comp_name}{the name of the output component column.
Default "component".}

\item{adj_mat_name}{the name of the output adjacency matrix. Default
"adj_mat".}

\item{size_name}{the name of the column denoting the component size. Default
"size".}
}
\description{
Break up an Adjacency Matrix by Components
}
