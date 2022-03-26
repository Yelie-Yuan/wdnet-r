##
## wdnet: Weighted directed network
## Copyright (C) 2022  Yelie Yuan, Panpan Zhang, and Jun Yan
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

#' Generate a directed preferential attachment network.
#'
#' @param edgelist A two column matrix represents the seed graph.
#' @param edgeweight A vector represents the weight of edges of the seed graph.
#'   Its length equals the number of edges of the seed graph. If NA, all the
#'   edges of the seed graph have weight 1.
#' @param nstep Number of steps when generating a network.
#' @param control A list of parameters to be used when generate network.
#'
#' @return A list with the following components: edgelist, edgeweight, out- and
#'   in-strength, number of edges per step (m), scenario of each new edge
#'   (1~alpha, 2~beta, 3~gamma, 4~xi, 5~rho). The edges in the seed graph 
#'   are denoted as scenario 0.
#' @export
#'
#' @examples
#' net <- rpanet_naive(nstep = 100, 
#'         control = panet.control(alpha = 0.4, beta = 0, gamma = 0.6))
#' net <- rpanet_naive(edgelist = matrix(c(1:8), ncol = 2), nstep = 10^3,
#'       control = panet.control(m_dist = stats::rpois,
#'       m_par = list(lambda = 1), m_const = 1,
#'       w_dist = stats::runif, w_par = list(min = 1, max = 10), w_const = 0))


rpanet_naive <- function(nstep = 10^3, edgelist = matrix(c(1, 2), ncol = 2), 
                          edgeweight = NA, control = panet.control()) {
  stopifnot("nstep must be greater than 0." = nstep > 0)
  stopifnot("alpha + beta + bamma + xi + rho must be less or equal to 1." = 
              control$alpha + control$beta + control$gamma + 
              control$xi + control$rho <= 1)
  temp <- c(edgelist)
  nnode <- max(temp)
  stopifnot("Nodes' index should be consecutive natural numbers start from 1." =
              sum(! duplicated(temp)) == nnode)
  if (is.na(edgeweight[1])) edgeweight[1:nrow(edgelist)] <- 1
  stopifnot(length(edgeweight) == nrow(edgelist))
  
  if (! is.numeric(control$m_dist)) {
    m <- do.call(control$m_dist, c(nstep, control$m_par)) + control$m_const
  }
  else m <- rep(control$m_dist + control$m_const, nstep)
  stopifnot("Number of new edges per step must be positive integers." = m %% 1 == 0)
  stopifnot("Number of new edges per step must be positive integers." = m > 0)
  sum_m <- sum(m)
  if (! is.numeric(control$w_dist)) {
    w <- do.call(control$w_dist, c(sum_m, control$w_par)) + control$w_const
  }
  else w <- rep(control$w_dist + control$w_const, sum_m)
  stopifnot("Edge weight must be greater than 0." = w > 0)
  
  strength <- nodeStrength_cpp(edgelist[, 1], edgelist[, 2], 
                               edgeweight, nnode, weighted = TRUE)
  outstrength <- instrength <- rep(0, 2 * sum_m + nnode)
  outstrength[1:nnode] <- strength$outstrength
  instrength[1:nnode] <- strength$instrength
  sumstrength <- sum(strength$outstrength)
  
  control_cpp <- c(control$alpha, control$beta, control$gamma, control$xi, 
                   control$delta_out, control$delta_in)
  
  ret <- rpanet_naive_cpp(nstep, control_cpp, m, w,
                           outstrength, instrength, sumstrength,
                           nnode)
  control$delta <- NULL
  list("edgelist" = rbind(edgelist, cbind(ret$start_node, ret$end_node) + 1), 
       "edgeweight" = c(edgeweight, w),
       "outstrength" = ret$outstrength, 
       "instrength" = ret$instrength, 
       "scenario" = c(rep(0, nrow(edgelist)), ret$scenario),
       "m" = m, 
       "control" = control)
}