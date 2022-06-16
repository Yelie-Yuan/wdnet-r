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

#' @importFrom utils modifyList
#' @importFrom stats rgamma rpois
NULL

#' Generate PA networks with non-linear preference functions
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
#'   newedge.control() + preference.control() + reciprocal.control()}. Under the
#'   default setup, in each step, a new edge of weight 1 is added from a new
#'   node \code{A} to an existing node \code{B} (\code{alpha} scenario), where
#'   \code{B} is chosen with probability proportional to its in-strength + 1.
#' @param directed Logical, whether to generate directed networks. If
#'   \code{FALSE}, the edge directions are omitted.
#' @param method Which method to use: \code{binary}, \code{naive},
#'   \code{edgesampler} or \code{nodelist}. For \code{nodelist} and
#'   \code{edgesampler} methods, the source preference function must be
#'   out-degree (out-strength) plus a nonnegative constant, the target
#'   preference function must be in-degree (in-strength) plus a nonnegative
#'   constant, \code{beta.loop} must be TRUE. Besides, \code{nodelist} method
#'   only works for unweighted networks, \code{edgeweight.control},
#'   \code{newedge.control}, \code{reciprocal.control} must set as default;
#'   \code{node.replace}, \code{snode.replace}, \code{tnode.replace} must be
#'   TRUE for \code{edgesampler} method.
#'
#'
#' @return A list with the following components: \code{edgelist},
#'   \code{edgeweight}, \code{strength} for undirected networks,
#'   \code{outstrength} and \code{instrength} for directed networks, number of
#'   new edges in each step \code{newedge}, controlling parameters
#'   \code{control}, node group \code{nodegroup} (if applicable) and edge
#'   scenario \code{scenario} (1~alpha, 2~beta, 3~gamma, 4~xi, 5~rho,
#'   6~reciprocal). The scenario of edges from \code{seednetwork} are denoted as
#'   0 unless specified in the seed network.
#'
#' @note The \code{nodelist} method implements the algorithm from Wan et al.
#'   (2017). The \code{edgesampler} first samples edges then find the
#'   source/target node of the sampled edge.
#'
#' @references \itemize{ \item Wan P, Wang T, Davis RA, Resnick SI (2017).
#'   Fitting the Linear Preferential Attachment Model. Electronic Journal of
#'   Statistics, 11(2), 3738–3780.}
#'
#' @export
#'
#' @examples
#' # Control edge scenario and edge weight through scenario.control()
#' # and edgeweight.control(), respectively, while keeping newedge.control(),
#' # preference.control() and reciprocal.control() as default.
#' set.seed(123)
#' control <- scenario.control(alpha = 0.5, beta = 0.5) +
#'     edgeweight.control(distribution = rgamma,
#'         dparams = list(shape = 5, scale = 0.2), shift = 0)
#' ret1 <- rpanet(nstep = 1e3, control = control)
#'
#' # In addition, set node groups and probability of creating reciprocal edges.
#' control <- control + reciprocal.control(group.prob = c(0.4, 0.6),
#'     recip.prob = matrix(runif(4), ncol = 2))
#' ret2 <- rpanet(nstep = 1e3, control = control)
#'
#' # Further, set the number of new edges in each step as Poisson(2) + 1 and use
#' # ret2 as a seed network.
#' control <- control + newedge.control(distribution = rpois,
#'     dparams = list(lambda = 2), shift = 1)
#' ret3 <- rpanet(nstep = 1e3, seednetwork = ret2, control = control)
#' 
rpanet <- function(nstep = 10^3, seednetwork = NULL,
                   control = NULL,
                   directed = TRUE,
                   method = c("binary", "naive", "edgesampler", "nodelist")) {
  method <- match.arg(method)
  stopifnot("nstep must be greater than 0." = nstep > 0)
  if (is.null(seednetwork)) {
    seednetwork <- list("edgelist" = matrix(1:2, ncol = 2),
                        "edgeweight" = NULL,
                        "nodegroup" = NULL)
  }
  nnode <- max(seednetwork$edgelist)
  nedge <- nrow(seednetwork$edgelist)
  if (is.null(seednetwork$edgeweight)) {
    seednetwork$edgeweight[1:nedge] <- 1
  }
  stopifnot(length(seednetwork$edgeweight) == nedge)
  if (is.null(seednetwork$nodegroup)) {
    seednetwork$nodegroup <- rep(1, nnode)
  }
  else {
    seednetwork$nodegroup <- as.integer(seednetwork$nodegroup)
    stopifnot('"nodegroup" of seednetwork is not valid.' =
                all(seednetwork$nodegroup > 0) &
                length(seednetwork$nodegroup) == nnode &
                max(seednetwork$nodegroup) <= length(control$reciprocal$group.prob))
  }

  control.default <- scenario.control() + edgeweight.control() +
    newedge.control() + preference.control() + reciprocal.control()
  if (! is.list(control)) {
    control <- structure(list(), class = "panet.control")
  }
  control <- control.default + control
  if (! control$newedge$node.replace) {
    control$scenario$beta.loop <- FALSE
    control$newedge$snode.replace <- control$newedge$tnode.replace <- FALSE
  }
  rm(control.default)

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
  temp <- 1
  if (! identical(control$reciprocal, reciprocal.control()$reciprocal)) {
    temp <- 2
  }
  if (is.function(control$edgeweight$distribution)) {
    w <- do.call(control$edgeweight$distribution,
                 c(sum_m * temp, control$edgeweight$dparams)) +
      control$edgeweight$shift
  } else {
    w <- rep(control$edgeweight$shift, sum_m * temp)
  }
  rm(temp)
  stopifnot("Edgeweight must be greater than 0." = w > 0)

  if (method == "nodelist" | method == "edgesampler") {
    if (directed) {
      stopifnot('Source preference must be out-degree plus a constant for "nodelist" and "edgesampler" methods.' = 
                  all(control$preference$sparams[1:2] == 1,
                      control$preference$sparams[3:4] == 0))
      stopifnot('Target preference must be in-degree plus a constant for "nodelist" and "edgesampler" methods.' = 
                  all(control$preference$tparams[1:2] == 0,
                      control$preference$tparams[3:4] == 1))
    }
    else {
      stopifnot('Preference must be degree plus a constant for "nodelist" and "edgesampler" methods.' = 
                  control$preference$params[1] == 1)
    }
    stopifnot('"reciprocal.control" must set as default for "nodelist" and "edgesampler" methods.' = 
                identical(control$reciprocal, reciprocal.control()$reciprocal))
    stopifnot('"beta.loop" must be TRUE for "nodelist" and "edgesampler" methods.' = 
                control$scenario$beta.loop)
  }
  if (method == "nodelist") {
    stopifnot('"edgeweight.control" must set as default for "nodelist" method.' = 
                identical(control$edgeweight, edgeweight.control()$edgeweight))
    stopifnot('Weight of existing edges must be 1 for "nodelist" method.' =
                all(seednetwork$edgeweight == 1))
    stopifnot('"newedge.control" must set as default for "nodelist" method.' = 
                identical(control$newedge, newedge.control()$newedge))
    return(rpanet_simple(nstep, seednetwork, control, directed,
                         m, sum_m, w, nnode, nedge, method))
  }
  if (method == "edgesampler") {
    stopifnot('"node.replace" must be TRUE for "edgesampler" method.' = 
                control$newedge$node.replace)
    stopifnot('"snode.replace" must be TRUE for "edgesampler" method.' = 
                control$newedge$snode.replace)
    stopifnot('"tnode.replace" must be TRUE for "edgesampler" method.' = 
                control$newedge$tnode.replace)
    return(rpanet_simple(nstep, seednetwork, control, directed,
                         m, sum_m, w, nnode, nedge, method))
  }
  return(rpanet_general(nstep, seednetwork, control, directed, 
                        m, sum_m, w, nnode, nedge, method))
}
