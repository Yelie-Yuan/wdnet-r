##
## wdnet: Weighted directed network
## Copyright (C) 2021  Panpan Zhang and Jun Yan
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
NULL

#' Parameter settings for function rpanet.
#'
#' @param alpha The probability of adding an edge from the new node to an
#'   existing node.
#' @param beta The probability of adding an edge between existing nodes.
#' @param gamma The probability of adding an edge from an existing node to a new
#'   node.
#' @param xi The probability of adding an edge between two new nodes.
#' @param rho The probability of introducing a new node with a self looped edge.
#' @param delta_out A tuning parameter related to growth rate. Probability of
#'   choosing an existing node as the source node of the newly added edge is
#'   proportional to nodes outstrength + outdelta.
#' @param delta_in A tuning parameter related to growth rate. Probability of
#'   choosing an existing node as the target node of the newly added edge is
#'   proportional to nodes instrength + indelta.
#' @param mdist Distribution function or a constant for number of newly added
#'   edges per step. The default value is 0.
#' @param mpar Additional parameters passed on to mdist.
#' @param m_c A constant add to mdist. The number of newly added edges per step
#'   then follows mdist(mpar) + m_c. The default value is 1.
#' @param wdist Dsitribution function or a constant for edge weights. The
#'   default value is 0.
#' @param wpar Additional parameters passed on to wdist.
#' @param w_c A constant add to wdist. The number of newly added edges per step
#'   then follows wdist(wpar) + w_c. The default value is 1.
#'
#' @return List of parameters.
#' @export


panet.control  <- function(alpha = 0.5, beta = 0.2, gamma = 0.3, xi = 0, rho = 0,
                           delta_out = 0.1, delta_in = 0.1, 
                           mdist = 0, 
                           mpar = list(), m_c = 1, 
                           wdist = 0, 
                           wpar = list(), w_c = 1) {
  ## set default value here
  list(alpha = alpha, beta = beta, gamma = gamma, xi = xi, rho = rho,
       delta_out = delta_out, delta_in = delta_in, 
       mdist = mdist, mpar = mpar, m_c = m_c,
       wdist = wdist, wpar = wpar, w_c = w_c)
}



#' Generate a growing network with preferential attachment.
#'
#' @param edgelist A two column matrix represents the seed graph.
#' @param edgeweight A vector represents the weight of edges of the seed graph.
#'   Its length equals the number of edges of the seed graph. If NA, all the
#'   edges of the seed graph have weight 1.
#' @param nsteps Number of steps when generate a network.
#' @param control A list of parameters to be used when generate network.
#'
#' @return A list with the following components: edgelist, edgeweight, out- and
#'   in-strength, number of edges per step (m), scenario of each new edge
#'   (1~alpha, 2~beta, 3~gamma, 4~xi, 5~rho).
#' @export
#'
#' @examples
#' net <- rpanet(nsteps = 100, 
#'         control = panet.control(alpha = 0.4, beta = 0, gamma = 0.6))
#' net <- rpanet(edgelist = matrix(c(1:8), ncol = 2), nsteps = 10^5,
#'       control = panet.control(mdist = stats::rpois,
#'       mpar = list(lambda = 1), m_c = 1,
#'       wdist = stats::runif, wpar = list(min = 1, max = 10), w_c = 0))


rpanet <- function(nsteps = 10^3, edgelist = matrix(c(1, 2), ncol = 2), 
                   edgeweight = NA,
                   control = panet.control()) {
  stopifnot("nsteps must be greater than 0." = nsteps > 0)
  stopifnot("alpha + beta + bamma + xi + rho must be less or equal to 1." = 
              control$alpha + control$beta + control$gamma + 
              control$xi + control$rho <= 1)
  if (is.na(edgeweight[1])) edgeweight[1:nrow(edgelist)] <- 1
  stopifnot(length(edgeweight) == nrow(edgelist))
  if (! is.numeric(control$mdist)) {
    m <- do.call(control$mdist, c(nsteps, control$mpar)) + control$m_c
  }
  else {
    m <- rep(control$mdist + control$m_c, nsteps)
  }
  stopifnot("Number of new edges per step must be positive integers." = m %% 1 == 0)
  stopifnot("Number of new edges per step must be positive integers." = m > 0)
  sum_m <- sum(m)
  if (! is.numeric(control$wdist)) {
    w <- do.call(control$wdist, c(sum_m, control$wpar)) + control$w_c
  }
  else {
    w <- rep(control$wdist + control$w_c, sum_m)
  }
  stopifnot("Edge weight must be greater than 0." = w > 0)
  
  exNodes <- max(c(edgelist))
  exEdges <- nrow(edgelist)
  exWeight <- sum(edgeweight)
  nEdges <- sum_m + exEdges
  
  scenario <- sample(1:5, size = sum_m, replace = TRUE,
                     prob = c(control$alpha, control$beta,
                              control$gamma, control$xi,
                              control$rho))
  scenario1 <- scenario == 1
  scenario4 <- scenario == 4
  
  noNewStart <- !((scenario > 3) | scenario1)
  noNewEnd <- scenario < 3
  totalNode <- endNode <- cumsum(c((scenario != 2) + scenario4)) + exNodes
  startNode <- totalNode - scenario4
  endNode[noNewEnd] <- 0
  startNode[noNewStart] <- 0
  nNodes <- totalNode[length(totalNode)]
  edgeweight <- c(edgeweight, w)
  weightIntv <- cumsum(c(0, edgeweight))
  
  temp_m <- cumsum(m[-nsteps])
  temp <- c(exWeight, weightIntv[temp_m + exEdges + 1])
  totalWeight <- rep(temp, m)
  temp <- c(exNodes, totalNode[temp_m])
  totalNode <- rep(temp, m)
  
  randOut <- runif(sum(noNewStart)) * (totalWeight + control$delta_out * totalNode)[noNewStart]
  randIn <- runif(sum(noNewEnd)) * (totalWeight + control$delta_in * totalNode)[noNewEnd]
  
  temp <- randOut < totalWeight[noNewStart]
  startEdge <- findInterval(randOut[temp], weightIntv, left.open = TRUE)
  if (! all(temp)) startNode[noNewStart][! temp] <- sapply(totalNode[noNewStart][! temp], 
                                                           function(i) sample(i, size = 1))
  startNode <- c(edgelist[, 1], startNode)
  temp <- which(startNode == 0)
  startNode <- nodes_cpp(startNode, startEdge, temp)
  
  temp <- randIn < totalWeight[noNewEnd]
  endEdge <- findInterval(randIn[temp], weightIntv, left.open = TRUE)
  if (! all(temp)) endNode[noNewEnd][! temp] <- sapply(totalNode[noNewEnd][!temp], 
                                                       function(i) sample(i, size = 1))
  endNode <- c(edgelist[, 2], endNode)
  temp <- which(endNode == 0)
  endNode <- nodes_cpp(endNode, endEdge, temp)
  
  outstrength <- strength_cpp(startNode, edgeweight, nNodes)
  instrength <- strength_cpp(endNode, edgeweight, nNodes)
  
  ret <- list(edgelist = cbind(startNode, endNode),
              edgeweight = edgeweight,
              outstrength = c(outstrength),
              instrength = c(instrength),
              scenario = scenario,
              m = m)
  return(ret)
}
