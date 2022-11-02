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
#' @importFrom RcppXPtrUtils checkXPtr
NULL

#' Generate PA networks.
#'
#' Generate preferential attachment (PA) networks with linear or non-linear
#' preference functions.
#'
#' @param nstep Number of steps when generating a network.
#' @param seednetwork A list represents the seed network. By default,
#'   \code{seednetwork} has one edge from node 1 to node 2 with weight 1. It
#'   consists of the following components: a two column matrix \code{edgelist}
#'   represents the edges; a vector \code{edgeweight} represents the weight of
#'   edges; an integer vector \code{nodegroup} represents the group of nodes.
#'   \code{nodegroup} is defined for directed networks, if \code{NULL}, all
#'   nodes from the seed graph are considered from group 1.
#' @param control A list of parameters that controls the PA network generation
#'   process. Defaults to an empty list, i.e., all the controlling parameters
#'   are set as default. For more details about available controlling
#'   parameters, see \code{rpacontrol.scenario}, \code{rpacontrol.newedge},
#'   \code{rpacontrol.edgeweight}, \code{rpacontrol.preference} and
#'   \code{rpacontrol.reciprocal}. Under the default setup, in each step, a new
#'   edge of weight 1 is added from a new node \code{A} to an existing node
#'   \code{B} (\code{alpha} scenario), where \code{B} is chosen with probability
#'   proportional to its in-strength + 1.
#' @param directed Logical, whether to generate directed networks. If
#'   \code{FALSE}, the edge directions are omitted.
#' @param method Which method to use: \code{binary}, \code{linear},
#'   \code{edgesampler} or \code{nodelist}. For \code{nodelist} and
#'   \code{edgesampler} methods, \code{beta.loop} must be \code{TRUE}; default
#'   preference functions must be used and \code{sparams = c(1, 1, 0, 0, a)},
#'   \code{tparams = c(0, 0, 1, 1, b)}, \code{param = c(1, c)}, where \code{a},
#'   \code{b} and \code{c} are non-negative constants; reciprocal edges and
#'   sampling without replacement are not considered, i.e., option
#'   \code{rpacontrol.reciprocal} must be set as default, \code{snode.replace},
#'   \code{tnode.replace} and \code{node.replace} must be \code{TRUE}. In
#'   addition, \code{nodelsit} method only works for unweighted networks and
#'   does not consider multiple edges, i.e., \code{rpacontrol.edgeweight} and
#'   \code{rpacontrol.newedge} must be set as default.
#'
#'
#' @return A list with the following components: \code{edgelist},
#'   \code{edgeweight}, \code{strength} for undirected networks,
#'   \code{outstrength} and \code{instrength} for directed networks, number of
#'   new edges in each step \code{newedge} (including reciprocal edges), control
#'   list \code{control}, node group \code{nodegroup} (if applicable) and edge
#'   scenario \code{scenario} (1~alpha, 2~beta, 3~gamma, 4~xi, 5~rho,
#'   6~reciprocal). The scenario of edges from \code{seednetwork} are denoted as
#'   0.
#'
#' @note The \code{nodelist} method implements the algorithm from Wan et al.
#'   (2017). The \code{edgesampler} first samples edges then find the
#'   source/target node of the sampled edge. If all the edges are of weight 1,
#'   the network can be considered as unweighted, node strength then equals node
#'   degree.
#'
#' @references \itemize{ \item Wan P, Wang T, Davis RA, Resnick SI (2017).
#'   Fitting the Linear Preferential Attachment Model. Electronic Journal of
#'   Statistics, 11(2), 3738–3780.}
#'
#' @export
#'
#' @examples
#' # Control edge scenario and edge weight through rpacontrol.scenario()
#' # and rpacontrol.edgeweight(), respectively, while keeping rpacontrol.newedge(),
#' # rpacontrol.preference() and rpacontrol.reciprocal() as default.
#' set.seed(123)
#' control <- rpacontrol.scenario(alpha = 0.5, beta = 0.5) +
#'     rpacontrol.edgeweight(distribution = rgamma,
#'         dparams = list(shape = 5, scale = 0.2), shift = 0)
#' ret1 <- rpanet(nstep = 1e3, control = control)
#'
#' # In addition, set node groups and probability of creating reciprocal edges.
#' control <- control + rpacontrol.reciprocal(group.prob = c(0.4, 0.6),
#'     recip.prob = matrix(runif(4), ncol = 2))
#' ret2 <- rpanet(nstep = 1e3, control = control)
#'
#' # Further, set the number of new edges in each step as Poisson(2) + 1 and use
#' # ret2 as a seed network.
#' control <- control + rpacontrol.newedge(distribution = rpois,
#'     dparams = list(lambda = 2), shift = 1)
#' ret3 <- rpanet(nstep = 1e3, seednetwork = ret2, control = control)
#' 
rpanet <- function(nstep = 10^3, seednetwork = list(
                    edgelist = matrix(c(1, 2), nrow = 1)),
                   control = list(),
                   directed = TRUE,
                   method = c("binary", "linear", "edgesampler", "nodelist")) {
  method <- match.arg(method)
  stopifnot("nstep must be greater than 0." = nstep > 0)
  nnode <- max(seednetwork$edgelist)
  stopifnot("Nodes must be consecutive integers starting from 1." = 
            min(seednetwork$edgelist) == 1 & 
            nnode == length(unique(c(seednetwork$edgelist))))
  stopifnot(ncol(seednetwork$edgelist) == 2)
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
                length(seednetwork$nodegroup) == nnode)
  }
  if (length(control$reciprocal$group.prob) > 0) {
    stopifnot('Length of "group.prob" in the control list in not valid.' = 
              max(seednetwork$nodegroup) <= length(control$reciprocal$group.prob))
  }
  
  control.default <- rpacontrol.scenario() + rpacontrol.edgeweight() +
    rpacontrol.newedge() + rpacontrol.reciprocal() + rpacontrol.preference()
  stopifnot(is.list(control))
  control <- structure(control, class = "rpacontrol")
  control <- control.default + control
  rm(control.default)
  if (control$preference$ftype == "customized") {
    if (directed) {
      RcppXPtrUtils::checkXPtr(ptr = control$preference$spref.pointer,
                               type = "double",
                               args = c("double", "double"))
      RcppXPtrUtils::checkXPtr(ptr = control$preference$tpref.pointer,
                               type = "double",
                               args = c("double", "double"))
    }
    else {
      RcppXPtrUtils::checkXPtr(ptr = control$preference$pref.pointer,
                               type = "double",
                               args = "double")
    }
  }
  
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
  sample.recip <- TRUE
  if (identical(control$reciprocal, rpacontrol.reciprocal()$reciprocal)) {
    sample.recip <- FALSE
  }
  if (is.function(control$edgeweight$distribution)) {
    w <- do.call(control$edgeweight$distribution,
                 c(sum_m * (1 + sample.recip), control$edgeweight$dparams)) +
      control$edgeweight$shift
  } else {
    w <- rep(control$edgeweight$shift, sum_m * (1 + sample.recip))
  }
  stopifnot("Edgeweight must be greater than 0." = w > 0)
  
  if ((! directed) & 
      ((! control$newedge$snode.replace) | (! control$newedge$tnode.replace))) {
    warning('"snode.replace" and "tnode.replace" are ignored for undirected networks.')
    control$newedge$snode.replace <- control$tnode.replace <- TRUE
  }
  if (directed & (! control$newedge$node.replace)) {
    warning('"node.replace" is ignored for directed networks.')
    control$newedge$node.replace <- TRUE
  }
  if (method == "nodelist" | method == "edgesampler") {
    stopifnot('"nodelist" and "edgesampler" methods require "default" preference functions.' = 
                control$preference$ftype == "default")
    if (directed) {
      stopifnot('Source preference must be out-degree plus a non-negative constant for "nodelist" and "edgesampler" methods.' = 
                  all(control$preference$sparams[1:2] == 1,
                      control$preference$sparams[3:4] == 0,
                      control$preference$sparams[5] >= 0))
      stopifnot('Target preference must be in-degree plus a non-negative constant for "nodelist" and "edgesampler" methods.' = 
                  all(control$preference$tparams[1:2] == 0,
                      control$preference$tparams[3:4] == 1,
                      control$preference$tparams[5] >= 0))
    }
    else {
      stopifnot('Preference must be degree plus a non-negative constant for "nodelist" and "edgesampler" methods.' = 
                  control$preference$params[1] == 1 & 
                    control$preference$params[2] >= 0)
    }
    stopifnot('"rpacontrol.reciprocal" must set as default for "nodelist" and "edgesampler" methods.' = 
                identical(control$reciprocal, rpacontrol.reciprocal()$reciprocal))
    stopifnot('"beta.loop" must be TRUE for "nodelist" and "edgesampler" methods.' = 
                control$scenario$beta.loop)
    if (method == "nodelist") {
      stopifnot('"rpacontrol.edgeweight" must set as default for "nodelist" method.' = 
                  identical(control$edgeweight, rpacontrol.edgeweight()$edgeweight))
      stopifnot('Weight of existing edges must be 1 for "nodelist" method.' =
                  all(seednetwork$edgeweight == 1))
      stopifnot('"rpacontrol.newedge" must set as default for "nodelist" method.' = 
                  identical(control$newedge, rpacontrol.newedge()$newedge))
    }
    if (method == "edgesampler") {
      if (directed) {
        stopifnot('"snode.replace" must be TRUE for "edgesampler" method.' = 
                    control$newedge$snode.replace)
        stopifnot('"tnode.replace" must be TRUE for "edgesampler" method.' = 
                    control$newedge$tnode.replace)
      }
      else {
        stopifnot('"node.replace" must be TRUE for "edgesampler" method.' = 
                    control$newedge$node.replace)
      }
    }
    return(rpanet_simple(nstep = nstep, seednetwork = seednetwork, 
                         control = control, directed = directed,
                         m = m, sum_m = sum_m, 
                         w = w, ex_node = nnode, 
                         ex_edge = nedge, method = method))
  }
  if ((! control$newedge$node.replace) & control$scenario$beta.loop) {
    control$scenario$beta.loop <- FALSE
    warning('"beta.loop" is set as FALSE since "node.replace" is FALSE.')
  }
  if (directed & 
      (! all(control$newedge$snode.replace, control$newedge$tnode.replace))) {
    if (all(control$preference$sparams[c(3, 5)] <= 0) |
        all(control$preference$tparams[c(1, 5)] <= 0) )
      stop("Source preference function and target preference function must be strictly positive when sampling source/target nodes without replacement.")
  }
  if ((! directed) & (! control$newedge$node.replace)) {
    if (control$preference$params[2] <= 0) {
      stop("Preference function must be strictly positive when sampling nodes without replacement.")
    }
  }
  return(rpanet_general(nstep = nstep, seednetwork = seednetwork, 
                        control = control, directed = directed, 
                        m = m, sum_m = sum_m, 
                        w = w, nnode = nnode, 
                        nedge = nedge, method = method, 
                        sample.recip = sample.recip))
}
