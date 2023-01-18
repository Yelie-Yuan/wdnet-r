##
## wdnet: Weighted directed network
## Copyright (C) 2023  Yelie Yuan, Tiandong Wang, Jun Yan and Panpan Zhang
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
#' @param initial.network A list represents the seed network. By default,
#'   \code{initial.network} has one edge from node 1 to node 2 with weight 1. It
#'   consists of the following components: a two column matrix \code{edgelist}
#'   represents the edges; a vector \code{edgeweight} represents the weight of
#'   edges; an integer vector \code{nodegroup} represents the group of nodes.
#'   \code{nodegroup} is defined for directed networks, if \code{NULL}, all
#'   nodes from the seed graph are considered from group 1.
#' @param control A list of parameters that controls the PA network generation
#'   process. Defaults to an empty list, i.e., all the controlling parameters
#'   are set as default. For more details about available controlling
#'   parameters, see \code{rpa_control_scenario}, \code{rpa_control_newedge},
#'   \code{rpa_control_edgeweight}, \code{rpa_control_preference} and
#'   \code{rpa_control_reciprocal}. Under the default setup, in each step, a new
#'   edge of weight 1 is added from a new node \code{A} to an existing node
#'   \code{B} (\code{alpha} scenario), where \code{B} is chosen with probability
#'   proportional to its in-strength + 1.
#' @param directed Logical, whether to generate directed networks. If
#'   \code{FALSE}, the edge directions are omitted.
#' @param method Which method to use: \code{binary}, \code{linear}, \code{bagx}
#'   or \code{bag}. For \code{bag} and \code{bagx} methods, \code{beta.loop}
#'   must be \code{TRUE}; default preference functions must be used and
#'   \code{sparams = c(1, 1, 0, 0, a)}, \code{tparams = c(0, 0, 1, 1, b)},
#'   \code{param = c(1, c)}, where \code{a}, \code{b} and \code{c} are
#'   non-negative constants; reciprocal edges and sampling without replacement
#'   are not considered, i.e., option \code{rpa_control_reciprocal} must be set
#'   as default, \code{snode.replace}, \code{tnode.replace} and
#'   \code{node.replace} must be \code{TRUE}. In addition, \code{nodelsit}
#'   method only works for unweighted networks and does not consider multiple
#'   edges, i.e., \code{rpa_control_edgeweight} and \code{rpa_control_newedge}
#'   must be set as default.
#'
#'
#' @return A list with the following components: \code{edgelist};
#'   \code{edgeweight}; number of new edges in each step \code{newedge}
#'   (reciprocal edges are not included); \code{node.attribute}, including node
#'   strengths, preference scores and node group (if applicable); control list
#'   \code{control}; edge scenario \code{scenario} (1~alpha, 2~beta, 3~gamma,
#'   4~xi, 5~rho, 6~reciprocal). The edges from \code{initial.network} are
#'   denoted as scenario 0.
#'
#' @note The \code{bianry} method implements binary search algorithm;
#'   \code{linear} represents linear search algorithm; \code{bag} method
#'   implements the algorithm from Wan et al. (2017); \code{bagx} puts all the
#'   edges into a bag, then samples edges and find the source/target node of the
#'   sampled edge.
#'
#' @references \itemize{ \item Wan P, Wang T, Davis RA, Resnick SI (2017).
#'   Fitting the Linear Preferential Attachment Model. Electronic Journal of
#'   Statistics, 11(2), 3738–3780.}
#'
#' @export
#'
#' @examples
#' # Control edge scenario and edge weight through rpa_control_scenario()
#' # and rpa_control_edgeweight(), respectively, while keeping rpa_control_newedge(),
#' # rpa_control_preference() and rpa_control_reciprocal() as default.
#' set.seed(123)
#' control <- rpa_control_scenario(alpha = 0.5, beta = 0.5) +
#'     rpa_control_edgeweight(distribution = rgamma,
#'         dparams = list(shape = 5, scale = 0.2), shift = 0)
#' ret1 <- rpanet(nstep = 1e3, control = control)
#'
#' # In addition, set node groups and probability of creating reciprocal edges.
#' control <- control + rpa_control_reciprocal(group.prob = c(0.4, 0.6),
#'     recip.prob = matrix(runif(4), ncol = 2))
#' ret2 <- rpanet(nstep = 1e3, control = control)
#'
#' # Further, set the number of new edges in each step as Poisson(2) + 1 and use
#' # ret2 as a seed network.
#' control <- control + rpa_control_newedge(distribution = rpois,
#'     dparams = list(lambda = 2), shift = 1)
#' ret3 <- rpanet(nstep = 1e3, initial.network = ret2, control = control)
#' 
rpanet <- function(nstep = 10^3, initial.network = list(
                    edgelist = matrix(c(1, 2), nrow = 1)),
                   control = list(),
                   directed = TRUE,
                   method = c("binary", "linear", "bagx", "bag")) {
  method <- match.arg(method)
  stopifnot("nstep must be greater than 0." = nstep > 0)
  nnode <- max(initial.network$edgelist)
  stopifnot("Nodes must be consecutive integers starting from 1." = 
            min(initial.network$edgelist) == 1 & 
            nnode == length(unique(c(initial.network$edgelist))))
  stopifnot(ncol(initial.network$edgelist) == 2)
  nedge <- nrow(initial.network$edgelist)
  if (is.null(initial.network$edgeweight)) {
    initial.network$edgeweight[1:nedge] <- 1
  }
  stopifnot(length(initial.network$edgeweight) == nedge)
  if (is.null(initial.network$nodegroup)) {
    initial.network$nodegroup <- rep(1, nnode)
  }
  else {
    initial.network$nodegroup <- as.integer(initial.network$nodegroup)
    stopifnot('"nodegroup" of initial.network is not valid.' =
                all(initial.network$nodegroup > 0) &
                length(initial.network$nodegroup) == nnode)
  }
  if (length(control$reciprocal$group.prob) > 0) {
    stopifnot('Length of "group.prob" in the control list in not valid.' = 
              max(initial.network$nodegroup) <= length(control$reciprocal$group.prob))
  }
  
  control.default <- rpa_control_scenario() + rpa_control_edgeweight() +
    rpa_control_newedge() + rpa_control_reciprocal() + rpa_control_preference()
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
  if (identical(control$reciprocal, rpa_control_reciprocal()$reciprocal)) {
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
  if (method == "bag" | method == "bagx") {
    stopifnot('"bag" and "bagx" methods require "default" preference functions.' = 
                control$preference$ftype == "default")
    if (directed) {
      stopifnot('Source preference must be out-degree plus a non-negative constant for "bag" and "bagx" methods.' = 
                  all(control$preference$sparams[1:2] == 1,
                      control$preference$sparams[3:4] == 0,
                      control$preference$sparams[5] >= 0))
      stopifnot('Target preference must be in-degree plus a non-negative constant for "bag" and "bagx" methods.' = 
                  all(control$preference$tparams[1:2] == 0,
                      control$preference$tparams[3:4] == 1,
                      control$preference$tparams[5] >= 0))
    }
    else {
      stopifnot('Preference must be degree plus a non-negative constant for "bag" and "bagx" methods.' = 
                  control$preference$params[1] == 1 & 
                    control$preference$params[2] >= 0)
    }
    stopifnot('"rpa_control_reciprocal" must set as default for "bag" and "bagx" methods.' = 
                identical(control$reciprocal, rpa_control_reciprocal()$reciprocal))
    stopifnot('"beta.loop" must be TRUE for "bag" and "bagx" methods.' = 
                control$scenario$beta.loop)
    if (method == "bag") {
      stopifnot('"rpa_control_edgeweight" must set as default for "bag" method.' = 
                  identical(control$edgeweight, rpa_control_edgeweight()$edgeweight))
      stopifnot('Weight of existing edges must be 1 for "bag" method.' =
                  all(initial.network$edgeweight == 1))
      stopifnot('"rpa_control_newedge" must set as default for "bag" method.' = 
                  identical(control$newedge, rpa_control_newedge()$newedge))
    }
    if (method == "bagx") {
      if (directed) {
        stopifnot('"snode.replace" must be TRUE for "bagx" method.' = 
                    control$newedge$snode.replace)
        stopifnot('"tnode.replace" must be TRUE for "bagx" method.' = 
                    control$newedge$tnode.replace)
      }
      else {
        stopifnot('"node.replace" must be TRUE for "bagx" method.' = 
                    control$newedge$node.replace)
      }
    }
    return(rpanet_simple(nstep = nstep, initial.network = initial.network, 
                         control = control, directed = directed,
                         m = m, sum_m = sum_m, 
                         w = w, ex_node = nnode, 
                         ex_edge = nedge, method = method))
  }
  if ((! control$newedge$node.replace) & control$scenario$beta.loop) {
    control$scenario$beta.loop <- FALSE
    warning('"beta.loop" is set as FALSE since "node.replace" is FALSE.')
  }
  return(rpanet_general(nstep = nstep, initial.network = initial.network, 
                        control = control, directed = directed, 
                        m = m, sum_m = sum_m, 
                        w = w, nnode = nnode, 
                        nedge = nedge, method = method, 
                        sample.recip = sample.recip))
}
