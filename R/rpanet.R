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
#' @param seednetwork A list represents the seed network. If \code{NULL},
#'   \code{seednetwork} will have one edge from node 1 to node 2 with weight 1.
#'   It consists of the following components: a two column matrix
#'   \code{edgelist} represents the edges; a vector \code{edgeweight} represents
#'   the weight of edges; an integer vector \code{nodegroup} represents the
#'   group of nodes. \code{nodegroup} is defined for directed networks, if
#'   \code{NULL}, all nodes from the seed graph are considered from group 1.
#' @param control A list of parameters that controls the PA generation process.
#'   The default value is \code{rpactl.scenario() + rpactl.edgeweight() +
#'   rpactl.newedge() + rpactl.preference() + rpactl.reciprocal()}. Under the
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
#'   only works for unweighted networks, \code{rpactl.edgeweight},
#'   \code{rpactl.newedge}, \code{rpactl.reciprocal} must set as default;
#'   \code{node.replace}, \code{snode.replace}, \code{tnode.replace} must be
#'   TRUE for \code{edgesampler} method.
#'
#'
#' @return A list with the following components: \code{edgelist},
#'   \code{edgeweight}, \code{strength} for undirected networks,
#'   \code{outstrength} and \code{instrength} for directed networks, number of
#'   new edges in each step \code{newedge} (reciprocal edges are not included),
#'   control list
#'   \code{control}, node group \code{nodegroup} (if applicable) and edge
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
#' # Control edge scenario and edge weight through rpactl.scenario()
#' # and rpactl.edgeweight(), respectively, while keeping rpactl.newedge(),
#' # rpactl.preference() and rpactl.reciprocal() as default.
#' set.seed(123)
#' control <- rpactl.scenario(alpha = 0.5, beta = 0.5) +
#'     rpactl.edgeweight(distribution = rgamma,
#'         dparams = list(shape = 5, scale = 0.2), shift = 0)
#' ret1 <- rpanet(nstep = 1e3, control = control)
#'
#' # In addition, set node groups and probability of creating reciprocal edges.
#' control <- control + rpactl.reciprocal(group.prob = c(0.4, 0.6),
#'     recip.prob = matrix(runif(4), ncol = 2))
#' ret2 <- rpanet(nstep = 1e3, control = control)
#'
#' # Further, set the number of new edges in each step as Poisson(2) + 1 and use
#' # ret2 as a seed network.
#' control <- control + rpactl.newedge(distribution = rpois,
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
  stopifnot("Nodes must be consecutive integers starting from 1." = 
            min(seednetwork$edgelist) == 1 & 
            nnode == length(unique(c(seednetwork$edgelist))))
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
  
  control.default <- rpactl.scenario() + rpactl.edgeweight() +
    rpactl.newedge() + rpactl.reciprocal() + rpactl.preference()
  if (! is.list(control)) {
    control <- structure(list(), class = "rpactl")
  }
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
  if (identical(control$reciprocal, rpactl.reciprocal()$reciprocal)) {
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
    stopifnot('"rpactl.reciprocal" must set as default for "nodelist" and "edgesampler" methods.' = 
                identical(control$reciprocal, rpactl.reciprocal()$reciprocal))
    stopifnot('"beta.loop" must be TRUE for "nodelist" and "edgesampler" methods.' = 
                control$scenario$beta.loop)
    if (method == "nodelist") {
      stopifnot('"rpactl.edgeweight" must set as default for "nodelist" method.' = 
                  identical(control$edgeweight, rpactl.edgeweight()$edgeweight))
      stopifnot('Weight of existing edges must be 1 for "nodelist" method.' =
                  all(seednetwork$edgeweight == 1))
      stopifnot('"rpactl.newedge" must set as default for "nodelist" method.' = 
                  identical(control$newedge, rpactl.newedge()$newedge))
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
