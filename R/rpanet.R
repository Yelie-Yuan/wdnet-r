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

#' @importFrom utils modifyList
#' @importFrom stats rgamma rpois
NULL

#' Generate a PA network with non-linear preference functions
#'
#' @param nstep Number of steps when generating a network.
#' @param seednetwork A list represents the seed network. If \code{NULL},
#'   \code{seednetwork} will have one edge from node 1 to node 2 with weight 1.
#'   It consists of the following components: a two column matrix
#'   \code{edgelist} represents the edges; a vector \code{edgeweight} represents
#'   the weight of edges; a integer vector \code{nodegroup} represents the group
#'   of each node. \code{nodegroup} is defined for directed networks, if
#'   \code{NULL}, all nodes from the seed graph are considered from group 1.
#' @param control A list of parameters that controls the PA generation process.
#'   The default value is \code{scenario.control() + edgeweight.control() +
#'   newedge.control() + preference.control() + reciprocal.control()}. By
#'   default, in each step, a new edge of weight 1 is added from a new node
#'   \code{A} to an existing node \code{B} (\code{alpha} scenario), where
#'   \code{B} is chosen with probability proportional to its in-strength + 1.
#' @param directed Logical, whether to generate directed networks. If
#'   \code{FALSE}, the edge directions are omitted.
#' @param method Which method to use when generating PA networks: "binary" or
#'   "naive". When \code{beta.loop = TRUE} is specified in
#'   \code{scenario.control}, \code{node.unique = FALSE}, \code{tnode.unique =
#'   FALSE}, \code{snode.unique = FALSE} are specified in \code{newedge.control}
#'   and \code{sparams = (1, 1, 0, 0, c), tparams = c(0, 0, 1, 1, c)} are use in
#'   \code{preference.control}, where \code{c} is a constant, the PA
#'   configuration is simple, the function will switch to a more efficient
#'   algorithm than \code{binary} and \code{naive}.
#'
#'
#' @return A list with the following components: edgelist, edgeweight, strength
#'   for undirected networks, out- and in-strength for directed networks,
#'   controlling parameters, node group (if applicable) and edge scenario
#'   (1~alpha, 2~beta, 3~gamma, 4~xi, 5~rho, 6~reciprocal). The scenario of
#'   edges from the seed network \code{seednetwork} are denoted as 0.
#'
#' @export
#'
#' @examples
#' # Control edge scenarios and edge weight through scenario.control()
#' # and edgeweight.control(), respectively, while keeping newedge.control(),
#' # preference.control() and reciprocal.control() as default.
#' set.seed(123)
#' control <- scenario.control(alpha = 0.5, beta = 0.5) +
#'     edgeweight.control(distribution = rgamma,
#'         dparams = list(shape = 5, scale = 0.2), shift = 0)
#' ret1 <- rpanet(nstep = 1e3, control = control)
#'
#' control <- control + reciprocal.control(group.prob = c(0.4, 0.6),
#'     recip.prob = matrix(runif(4), ncol = 2))
#' ret2 <- rpanet(nstep = 1e3, control = control)
#'
#' control <- control + newedge.control(distribution = rpois,
#'     dparams = list(lambda = 2), shift = 1)
#' ret3 <- rpanet(nstep = 1e3, seednetwork = ret2, control = control)
#' 
rpanet <- function(nstep = 10^3, seednetwork = NULL,
                   control = NULL,
                   directed = TRUE,
                   method = c("binary", "naive")) {
  stopifnot("nstep must be greater than 0." = nstep > 0)
  method <- match.arg(method)
  if (is.null(seednetwork)) {
    seednetwork <- list("edgelist" = matrix(1:2, ncol = 2),
               "edgeweight" = NULL,
               "nodegroup" = NULL)
  }
  temp <- c(seednetwork$edgelist)
  nnode <- max(temp)
  stopifnot("Nodes' index should be consecutive numbers starting from 1." =
              sum(! duplicated(temp)) == nnode)
  rm(temp)
  if (is.null(seednetwork$nodegroup)) {
    seednetwork$nodegroup <- rep(1, nnode)
  }
  else {
    seednetwork$nodegroup <- as.integer(seednetwork$nodegroup)
    stopifnot('"nodegroup" is not valid.' =
                all(seednetwork$nodegroup > 0) &
                length(seednetwork$nodegroup) == nnode &
                max(seednetwork$nodegroup) > length(control$reciprocal.control$group.prob))
  }
  nedge <- nrow(seednetwork$edgelist)
  if (is.null(seednetwork$edgeweight)) {
    seednetwork$edgeweight[1:nedge] <- 1
  }
  stopifnot(length(seednetwork$edgeweight) == nedge)
  
  default.control <- scenario.control() + edgeweight.control() +
    newedge.control() + preference.control() + reciprocal.control()
  if (! is.list(control)) {
    control <- list()
  }
  control <- modifyList(default.control, control, keep.null = TRUE)
  rm(default.control)
  
  if (is.function(control$newedge$distribution)) {
    m <- do.call(control$newedge$distribution, c(nstep, control$newedge$dparams)) + 
      control$newedge$shift
  } else {
    m <- rep(control$newedge$shift, nstep)
  }
  stopifnot("Number of new edges per step must be positive integers." =
              all(m %% 1 == 0))
  stopifnot("Number of new edges per step must be positive integers." =
              all(m > 0))
  sum_m <- sum(m)
  if (is.function(control$edgeweight$distribution)) {
    w <- do.call(control$edgeweight$distribution,
                 c(sum_m * 2, control$edgeweight$dparams)) +
      control$edgeweight$shift
  } else {
    w <- rep(control$edgeweight$shift, sum_m * 2)
  }
  stopifnot("Edgeweight must be greater than 0." = w > 0)
  
  simplecase <- FALSE
  if (! any(control$newedge$node.unique, control$newedge$snode.unique, 
            control$newedge$tnode.unique)) {
    if (control$scenario$beta.loop & is.null(control$reciprocal$group.prob) & 
        is.null(control$reciprocal$recip.prob)) {
      if (directed) {
        if (all(control$preference$sparams[1:2] == 1,
                control$preference$sparams[3:4] == 0,
                control$preference$tparams[1:2] == 0,
                control$preference$tparams[3:4] == 1)) {
          simplecase <- TRUE
        }
      }
      else {
        if(control$preference$params[1] == 1) {
          simplecase <- TRUE
        }
      }
    }
  }
  if (simplecase) {
    cat("PA generation setup is simple. Switch to a more efficient method. \n")
    ret <- rpanet_simple(nstep, seednetwork, control, directed, 
                         m, sum_m, w, nnode, nedge)
  }
  else {
    ret <- rpanet_general(nstep, seednetwork, control, directed, 
                          m, sum_m, w, nnode, nedge, method)
  }
  return(ret)
}
