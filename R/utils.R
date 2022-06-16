##
## wdnet: Weighted directed network
## Copyright (C) 2022  Yelie Yuan, Tiandong Wang, Jun Yan and Panpan Zhang
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

#' Convert adjacency matrix to edgelist and edgeweight.
#'
#' @param adj Adjacency matrix of a network.
#' @param directed Logical, whether the network is directed.
#'
#' @return A list of edgelist and edgeweight.
#'
adj_to_edges <- function(adj, directed = TRUE) {
  if (! directed) {
    stopifnot('"adj" must be symmetric if the network is undirected.' = 
                isSymmetric(adj))
  }
  mode <- ifelse(directed, "directed", "undirected")
  g <- igraph::graph_from_adjacency_matrix(adj, mode = mode, 
                                           weighted = TRUE, diag = TRUE)
  edgelist <- igraph::as_edgelist(g)
  edgeweight <- igraph::E(g)$weight
  return(list("edgelist" = edgelist, 
              "edgeweight" = edgeweight))
}

#' Convert edgelist and edgeweight to adjacency matrix.
#'
#' @param edgelist A two column matrix represents edges. If \code{NULL},
#'   \code{edgelist} and \code{edgeweight} will be extracted from the adjacency
#'   matrix \code{adj}.
#' @param edgeweight A vector represents the weight of edges. If \code{NULL},
#'   all the edges are considered have weight 1.
#' @param directed Logical, whether the network is directed.
#'
#' @return An adjacency matrix.
#'
edges_to_adj <- function(edgelist, edgeweight = NULL, directed = TRUE) {
  nnode <- max(edgelist)
  adj <- matrix(0, nrow = nnode, ncol = nnode)
  if (is.null(edgeweight)) {
    edgeweight <- rep(1, nrow(edgelist))
  }
  adj <- fillWeight_cpp(adj, edgelist - 1, edgeweight)
  if (! directed) {
    adj <- adj + t(adj)
    diag(adj) <- diag(adj) / 2
  }
  return(adj)
}