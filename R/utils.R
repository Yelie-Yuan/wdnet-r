##
## wdnet: Weighted directed network
## Copyright (C) 2024  Yelie Yuan, Tiandong Wang, Jun Yan and Panpan Zhang
## Jun Yan <jun.yan@uconn.edu>
##
## This file is part of the R package wdnet.
##
## The R package wdnet is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package wdnet is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

#' @importFrom igraph graph_from_adjacency_matrix as_edgelist E
NULL

#' Converts an adjacency matrix to edgelist and edgeweight using the
#' \code{igraph} package.
#'
#' @param adj Adjacency matrix of a network.
#' @param directed Logical, whether the network is directed. This value is
#'   passed to \code{igraph::graph_from_adjacency_matrix()}.
#' @param weighted Logical, whether the network is weighted.
#'
#' @return A list of edgelist, edgeweight and directed.
#'
#' @keywords internal
#' 
adj_to_edgelist <- function(adj, directed = TRUE, weighted = TRUE) {
  if (!directed) {
    if (! isSymmetric(adj)) {
      directed <- TRUE
      cat('Returned network is directed because "adj" is asymmetric.\n\n')
    }
  }
  if (!weighted) {
    weighted <- NULL
  }
  g <- igraph::graph_from_adjacency_matrix(adj,
    mode = ifelse(directed, "directed", "undirected"),
    weighted = weighted,
    diag = TRUE
  )
  edgelist <- igraph::as_edgelist(g)
  edgeweight <- igraph::E(g)$weight
  list(
    "edgelist" = edgelist,
    "edgeweight" = edgeweight,
    "directed" = directed
  )
}

#' Convert edgelist and edgeweight to adjacency matrix.
#'
#' @param edgelist A two column matrix representing edges.
#' @param edgeweight A vector representing the weight of edges. If \code{NULL},
#'   all edges are considered to have a weight of 1.
#' @param directed Logical, whether the network is directed.
#'
#' @return Returns an adjacency matrix.
#'
#' @keywords internal
#' 
edgelist_to_adj <- function(edgelist, edgeweight, directed = TRUE) {
  nnode <- max(edgelist)
  adj <- matrix(0, nrow = nnode, ncol = nnode)
  if (missing(edgeweight)) {
    edgeweight <- rep(1, nrow(edgelist))
  }
  if (length(edgelist) == 2) {
    edgelist <- matrix(edgelist, ncol = 2)
    edgeweight <- c(edgeweight, 0)
  }
  adj <- fill_weight_cpp(adj, edgelist - 1, edgeweight)
  if (!directed) {
    adj <- adj + t(adj)
    diag(adj) <- diag(adj) / 2
  }
  return(adj)
}
