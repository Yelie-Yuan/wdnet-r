##
## wdnet: Weighted directed network
## Copyright (C) 2020  Panpan Zhang and Jun Yan
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

#' @importFrom stats runif rpois
#' @importFrom igraph graph_from_edgelist degree
NULL

#' Parameter settings for function rpanet.
#'
#' @param alpha The probability of adding an edge from the new node to an
#'   existing node.
#' @param beta The probability of adding an edge between two existing nodes.
#' @param gamma The probability of adding an edge from an existing node to a
#'   new node. 1 - alpha - beta - gamma then represents the probability of
#'   adding an edge between two newly added node.
#' @param xi The probability of adding an edge between two new nodes.
#' @param delta_out is a tuning parameter related to growth rate. Probability of
#'   choosing an existing node as the source node of the newly added edge is
#'   proportional to nodes outstrength + outdelta.
#' @param delta_in is a tuning parameter related to growth rate. Probability of
#'   choosing an existing node as the target node of the newly added edge is
#'   proportional to nodes instrength + indelta.
#' @param mdist Distribution function for number of newly added edges per step.
#' @param mpar Additional parameters passed on to mdist.
#' @param wdist Dsitribution function for edge weights.
#' @param wpar Additional parameters passed on to wdist.
#' @param ... Additional arguments
#' 
#' @return List of parameters.
#' @export


panet.control  <- function(alpha = 0.5, beta = 0, gamma = 0.5, xi = 0, 
                           delta_out = 0.1, delta_in = 0.1, 
                           mdist = stats::rpois, 
                           mpar = list(lambda = 1), 
                           wdist = stats::runif, 
                           wpar = list(min = 1, max = 1), ...) {
  ## set default value here
  ## how to set m as poisson(lambda) + 1 and m > 0.
  list(alpha = alpha, beta = beta, gamma = gamma, xi = xi,
       delta_out = delta_out, delta_in = delta_in, 
       mdist = mdist, mpar = mpar, 
       wdist = wdist, wpar = wpar, 
       batsize = round(mean(do.call(mdist, c(100, mpar)))), ...)
}



#' Generate a growing network with preferential attachment.
#'
#' @param edgelist A two column matrix represents the starting graph.
#' @param edgeweight A vector represents the weight of each edges in edgelist.
#' @param nsteps Number of steps when generate a network.
#' @param directed Logical. Whether to generate a directed graph.
#' @param control A list of parameters to be used when generate network.
#' @param ... Additional arguments
#'
#' @return A three column matrix represents the generated network with given
#'   parameters. The first two columns represent the source and target node,
#'   respectively, for directed network. For undirected network, the first two
#'   columns represents the nodes at either end of edges. The third column
#'   represents the weight of edges.
#' @export
#'
#' @examples
#' net <- rpanet(nsteps = 100, directed = FALSE,
#'         control = panet.control(alpha = 0.4, beta = 0, gamma = 0.6))
#' net <- rpanet(edgelist = matrix(c(1:8), ncol = 2), nsteps = 100, 
#'       control = panet.control(mdist = stats::rbinom,
#'       mpar = list(size = 5, prob = 0.2),
#'       wdist = stats::runif, wpar = list(min = 1, max = 10)))


rpanet <- function(nsteps = 10^3, edgelist = matrix(c(1, 2), ncol = 2), 
                   edgeweight = NA, directed = TRUE, 
                   control = panet.control(), ...) {
  stopifnot(nsteps > 0)
  if (control$alpha + control$beta + control$gamma + control$xi> 1) {
    stop("Alpha + beta + bamma + xi must be less or equal to 1.")
  }
  if (is.na(edgeweight[1])) edgeweight[1:dim(edgelist)[1]] <- 1
  
  if (! directed) stopifnot(control$delta_in == control$delta_out)
  
  nnode <- tnode <- max(c(edgelist))
  m <- do.call(control$mdist, c(nsteps, control$mpar)) + 1
  temp_m <- sum(m)
  w <- do.call(control$wdist, c(temp_m, control$wpar))

  g <- igraph::graph_from_edgelist(edgelist, directed = directed)
  outstrength <- rep(0, 2 * temp_m + nnode) + control$delta_out
  instrength <- rep(0, 2 *temp_m + nnode) + control$delta_in
  outstrength[1:nnode] <- outstrength[1:nnode] + igraph::degree(g, mode = 'out')
  instrength[1:nnode] <- instrength[1:nnode] + igraph::degree(g, mode = 'in')
  
  control_cpp <- c(control$alpha, control$beta, control$gamma)
  ret <- rpanet_cpp(nsteps,
                    control_cpp, directed, m, w,
                    outstrength, instrength,
                    nnode, tnode)
  ret2 <- list()
  ret2$edgelist <- rbind(edgelist, cbind(ret$startnode, ret$endnode))
  ret2$edgeweight <- c(edgeweight, w)
  ret2$outstrength <- ret$outstrength - control$delta_out
  ret2$instrength <- ret$instrength - control$delta_in
  return(ret2)
}
