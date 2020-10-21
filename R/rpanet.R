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

#' Parameter settings for function rpanet.
#'
#' @param alpha The probability of adding an edge from the new node to an
#'   existing node.
#' @param beta The probability of adding an edge between two existing nodes.
#' @param gamma The probability of adding an edge from an existing node to a
#'   new node. 1 - alpha - beta - gamma then represents the probability of
#'   adding an edge between two newly added node.
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
#' @param ...
#' 
#' @return List of parameters.
#' @export

panet.control  <- function(alpha = 0.5, beta = 0, gamma = 0.5,
                           delta_out = 0.1, delta_in = 0.1, 
                           mdist = stats::rpois, 
                           mpar = list(lambda = 1), 
                           wdist = stats::runif, 
                           wpar = list(min = 1, max = 1), ...) {
  ## set default value here
  ## how to set m as poisson(lambda) + 1 and m > 0.
  list(alpha = alpha, beta = beta, gamma = gamma, 
       delta_out = delta_out, delta_in = delta_in, 
       mdist = mdist, mpar = mpar, 
       wdist = wdist, wpar = wpar, 
       batsize = round(mean(do.call(mdist, c(100, mpar)))), ...)
}

#' Generate a growing network with preferential attachment.
#'
#' @param nsteps Number of steps when generate a network.
#' @param directed Logical. Whether to generate a directed graph.
#' @param control A list of parameters to be used when generate network.
#' @param ...
#'
#' @return A three column matrix represents the generated network with given
#'   parameters. The first two columns represent the source and target node,
#'   respectively, for directed network. For undirected network, the first two
#'   columns represents the nodes at either end of edges. The third column
#'   represents the weight of edges.
#' @export
#'
#' @examples
#' net <- rpanet(10^3, directed = FALSE,
#'         control = panet.control(alpha = 0.4, beta = 0, gamma = 0.2))
#' net <- rpanet(10^3, control = panet.control(mdist = stats::rbinom,
#'       mpar = list(size = 5, prob = 0.2),
#'       wdist = stats::runif, wpar = list(min = 1, max = 10)))

rpanet <- function(nsteps, directed = TRUE, 
                   control = panet.control(), # input
                   ...) {
  if (control$alpha + control$beta + control$gamma > 1) {
    stop("Alpha + beta + bamma must be less or equal to 1.")
  }
  ## set output container
  outmat <- matrix(NA, nrow = nsteps * control$batsize, ncol = 3) ## to be changed
  outmat <- matrix(NA, nrow = 2, ncol = 3)
  outmat[1, ] <- 1
  numrow <- 1
  numnode <- tnumnode <- 1
  outs <- ins <- 1
  if (! directed) outs <- ins <- outs + ins
  for (i in 1:nsteps) {
    m  <- do.call(control$mdist, c(1, control$mpar))
    if (m > 0) {
      w  <- do.call(control$wdist, c(m, control$wpar))
      v1 <- v2 <- rep(NA, m)
      for (j in 1:m) {
        u <- runif(1, min = 0, max = 1)
        if (u < control$alpha) {
          v1[j] <- numnode + 1
          v2[j] <- sample(tnumnode, 1, prob = ins + control$delta_in)
          numnode <- numnode + 1
        }
        else if (u < control$alpha + control$beta) {
          v1[j] <- sample(tnumnode, 1, prob = outs + control$delta_out)
          v2[j] <- sample(tnumnode, 1, prob = ins + control$delta_in)
        }
        else if (u < control$alpha + control$beta + control$gamma) {
          v1[j] <- sample(tnumnode, 1, prob = outs + control$delta_out)
          v2[j] <- numnode + 1
          numnode <- numnode + 1
        }
        else {
          v1[j] <- numnode + 1
          v2[j] <- numnode + 2
          numnode <- numnode + 2
        }
        numrow <- numrow + 1
      }
      outs[is.na(outs[1:numnode])] <- 0
      ins[is.na(ins[1:numnode])] <- 0
      for (j in 1:m) {
        outs[v1[j]] <- outs[v1[j]] + w[j]
        ins[v2[j]] <- ins[v2[j]] + w[j]
      }
      if (! directed) {
        for (j in 1:m) {
          outs[v2[j]] <- outs[v2[j]] + w[j]
          ins[v1[j]] <- ins[v1[j]] + w[j]
        }
      }
      ## check if outmat is full
      ## increase size if full
      if (numrow > dim(outmat)[1]) {
        outmat <- rbind(outmat, matrix(NA, nrow = nsteps * control$batsize, ncol = 3))
      }
      ## fill
      outmat[(numrow - m + 1):numrow, ] <- cbind(v1, v2, w)
      tnumnode <- numnode
    }
  }
  outmat[1:numrow, ]
}
