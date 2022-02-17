##
## wdnet: Weighted directed network
## Copyright (C) 2022  Yelie Yuan, Panpan Zhang and Jun Yan
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

#' @importFrom stats weighted.mean
#' @importFrom wdm wdm
#' @importFrom igraph graph_from_adjacency_matrix as_edgelist
NULL

## Directed assortativity coefficient

#' Compute the assortativity coefficient of a weighted and directed network.
#'
#' @param adj is an adjacency matrix of a weighted and directed network.
#' @param type which type of assortativity coefficient to compute: "out-in" (default), 
#' "in-in", "out-out" or "in-out"?
#'
#' @return a scalar of assortativity coefficient
#'
#' @references
#' \itemize{
#' \item Foster, J.G., Foster, D.V., Grassberger, P. and Paczuski, M. (2010). Edge direction 
#' and the structure of networks. \emph{Proceedings of the National Academy of Sciences of the
#' United States}, 107(24), 10815--10820.
#' \item Yuan, Y. Zhang, P. and Yan, J. (2020+). Assortativity coefficients for 
#' weighted and directed networks
#' }
#'
#' @note 
#' When the adjacency matrix is binary (i.e., directed but unweighted networks), \code{dw_assort}
#' returns the assortativity coefficient proposed in Foster et al. (2010).
#'
#' @examples
#' ## Generate a network according to the Erd\"{o}s-Renyi model of order 20
#' ## and parameter p = 0.3
#' edge_ER <- rbinom(400,1,0.3)
#' weight_ER <- sapply(edge_ER, function(x) x*sample(3,1))
#' adj_ER <- matrix(weight_ER,20,20)
#' system.time(myassort <- dw_assort(adj_ER, type = "out-in"))
#' myassort
#' 
#' @export

dw_assort <- function(adj, type = c("out-in", "in-in", "out-out", "in-out")) {
  stopifnot(dim(adj)[1] == dim(adj)[2])
  ## determine the location of edges in the network
  in_str <- colSums(adj)
  out_str <- rowSums(adj)
  vert_from <- unlist(apply(adj, 2, function(x){which(x > 0)}))
  number_to <- apply(adj, 2, function(x){length(which(x > 0) == TRUE)})
  temp_to <- cbind(seq(1:dim(adj)[1]),number_to)
  vert_to <- rep(temp_to[,1],temp_to[,2])
  weight <- adj[which(adj > 0)]
  type  <- match.arg(type)
  .type <- unlist(strsplit(type, "-"))
  x <- switch(.type[1], "out" = out_str, "in" = in_str)[vert_from]
  y <- switch(.type[2], "out" = out_str, "in" = in_str)[vert_to]
  weighted.cor <- function(x, y, w) {
    mean_x <- stats::weighted.mean(x, w)
    mean_y <- stats::weighted.mean(y, w)
    var_x <- sum((x - mean_x)^2 * w)
    var_y <- sum((y - mean_y)^2 * w)
    return(sum(w * (x - mean_x) * (y - mean_y)) / 
             sqrt(var_x * var_y))
  }
  return(weighted.cor(x, y, weight))
}

#' Return four directed, weighted assortativity coefficient with a given
#' network.
#'
#' @param edgelist A two column matrix represents edges.
#' @param directed Logical. Whether the edges will be considered as directed. If
#'   FALSE, the input network will be considered as undirected.
#' @param edgeweight A vector represents the weight of edges.
#'
#' @return Assortativity coefficient for undirected network, or four directed
#'   assortativity coefficients for directed network.
#' @export
#'
#' @examples
#' net <- rpanet(nstep = 10^3)
#' result <- edge_assort(net$edgelist, directed = TRUE)
edge_assort <- function(edgelist, edgeweight = NA, directed = TRUE) {
  if (! directed) {
    edgelist <- rbind(edgelist, edgelist[, c(2, 1)])
    edgeweight <- c(edgeweight, edgeweight)
  }
  temp <- range(c(edgelist))
  stopifnot("Node index should start from 1." = temp[1] == 1)
  nnode <- temp[2]
  sourceNode <- edgelist[, 1]
  targetNode <- edgelist[, 2]
  if (is.na(edgeweight[1])) {
    temp <- nodeStrength_cpp(start_node = sourceNode, 
                             end_node = targetNode, 
                             nnode = nnode, 
                             weight = 1,
                             weighted = FALSE)
    edgeweight <- rep(1, length(sourceNode))
  } else {
    temp <- nodeStrength_cpp(start_node = sourceNode, 
                             end_node = targetNode, 
                             nnode = nnode, 
                             weight = edgeweight,
                             weighted = TRUE)
  }
  outs <- temp$outstrength
  ins <- temp$instrength
  sourceOut <- outs[sourceNode]
  targetIn <- ins[targetNode]
  if (! directed) {
    return(wdm::wdm(x = sourceOut, y = targetIn, 
                    weights = edgeweight, method = 'pearson'))
  }
  sourceIn <- ins[sourceNode]
  targetOut <- outs[targetNode]
  return(list('out-out' = wdm::wdm(x = sourceOut, y = targetOut,
                                   weights = edgeweight, method = 'pearson'), 
              'out-in' = wdm::wdm(x = sourceOut, y = targetIn, 
                                  weights = edgeweight, method = 'pearson'), 
              'in-out' = wdm::wdm(x = sourceIn, y = targetOut,
                                  weights = edgeweight, method = 'pearson'),
              'in-in' = wdm::wdm(x = sourceIn, y = targetIn,
                                 weights = edgeweight, method = 'pearson')))
}

#' Assortativity coefficients between features of a weighted and directed
#' network.
#'
#' @param edgelist A two column matrix represents the directed edges. If
#'   \code{edgelist} and \code{edgeweight} are \code{NA}, the adjacency matrix
#'   \code{adj} will be used.
#' @param edgeweight Vector, represents the weight of edges.
#' @param adj The adjacency matrix of a weighted directed network. If \code{NA},
#'   \code{edgelist} and \code{edgeweight} will be used.
#' @param feature1 Vector, represents feature values of node 1, node 2, etc,.
#'   Number of nodes \code{= length(feature1) = length(feautre2)}. If \code{NA},
#'   out-strength will be used.
#' @param feature2 Vector, represents feature values of node 1, node 2, etc,. If
#'   \code{NA}, in-strength will be used.
#'
#' @return Directed weighted assortativity coefficients between source nodes'
#'   \code{feature1} (or \code{feature2}) and target nodes' \code{feature2}(or
#'   \code{feature1}).
#' @export
#'
#' @examples
#' adj <- matrix(rbinom(400, 1, 0.2) * sample(1:3, 400, replace = TRUE), 20, 20)
#' feature1 <- runif(20)
#' feature2 <- abs(rnorm(20))
#' ret <- dw_feature_assort(adj = adj, feature1 = feature1, feature2 = feature2)
#' 
dw_feature_assort <- function(adj = NA, edgelist = NA, edgeweight = NA, 
                              feature1 = NA, feature2 = NA) {
  if (is.na(edgelist)[1] & is.na(edgeweight)[1]) {
    temp_adj <- (adj > 0) * 1
    g <- igraph::graph_from_adjacency_matrix(temp_adj)
    rm(temp_adj)
    edgelist <- igraph::as_edgelist(g)
    adj <- t(adj)
    edgeweight <- adj[adj > 0]
  }
  nNodes <- max(edgelist)
  if (is.na(feature1)[1] | is.na(feature2)[1]) {
    temp <- nodeStrength_cpp(start_node = edgelist[, 1], 
                             end_node = edgelist[, 2], 
                             weight = edgeweight, 
                             nnode = nNodes, 
                             weighted = TRUE)
  }
  if (is.na(feature1)[1]) {
    feature1 <- temp$outstrength
  }
  if (is.na(feature2)[1]) {
    feature2 <- temp$instrength
  }
  stopifnot(length(feature1) == nNodes & length(feature2) == nNodes)
  ret <- list()
  ret[[paste('feature1', 'feature1', sep = '-')]] <- wdm::wdm(feature1[edgelist[, 1]],
                                                              feature1[edgelist[, 2]],
                                                              weights = edgeweight, 
                                                              method = 'pearson')
  ret[[paste('feature1', 'feature2', sep = '-')]] <- wdm::wdm(feature1[edgelist[, 1]],
                                                              feature2[edgelist[, 2]],
                                                              weights = edgeweight, 
                                                              method = 'pearson')
  ret[[paste('feature2', 'feature1', sep = '-')]] <- wdm::wdm(feature2[edgelist[, 1]],
                                                              feature1[edgelist[, 2]],
                                                              weights = edgeweight, 
                                                              method = 'pearson')
  ret[[paste('feature2', 'feature2', sep = '-')]] <- wdm::wdm(feature2[edgelist[, 1]],
                                                              feature2[edgelist[, 2]],
                                                              weights = edgeweight, 
                                                              method = 'pearson')
  return(ret)
}
