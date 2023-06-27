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

#' @importFrom stats rgamma rpois
NULL

#' Generate PA networks.
#'
#' Generate preferential attachment (PA) networks with linear or non-linear
#' preference functions.
#'
#' @param nstep Number of steps.
#' @param initial.network A \code{wdnet} object or a list representing the
#'   initial network. By default, \code{initial.network} has one directed edge
#'   from node 1 to node 2 with weight 1. It can contain the following components:
#'   a two-column matrix \code{edgelist} representing the edges; a vector
#'   \code{edgeweight} representing the weight of edges; a logical argument
#'   \code{directed} indicating whether the initial network is directed. If
#'   \code{edgeweight} is not specified, all edges from the initial network are
#'   assumed to have weight 1. In addition, an integer vector
#'   \code{nodegroup} can be added to the list for specifing node groups;
#'   \code{nodegroup} is defined for directed networks, if \code{NULL}, all
#'   nodes from the seed network are assumed to be in group 1.
#' @param control An \code{rpacontrol} object controlling the PA network
#'   generation process. If not specified, all the control parameters will be
#'   set to default. For more details, see \code{rpa_control_scenario()},
#'   \code{rpa_control_newedge()}, \code{rpa_control_edgeweight()},
#'   \code{rpa_control_preference} and \code{rpa_control_reciprocal()}. Under
#'   the default setup, at each step, a new edge of weight 1 is added from a new
#'   node \code{A} to an existing node \code{B} (\code{alpha} scenario), where
#'   \code{B} is chosen with probability proportional to its in-strength + 1.
#' @param method Which method to use: \code{binary}, \code{linear}, \code{bagx}
#'   or \code{bag}. For \code{bag} and \code{bagx} methods, \code{beta.loop}
#'   must be \code{TRUE}, default preference functions must be used, and
#'   \code{sparams} should be set to \code{c(1, 1, 0, 0, a)}, \code{tparams} to
#'   \code{c(0, 0, 1, 1, b)}, and \code{param} to \code{c(1, c)}, where
#'   \code{a}, \code{b}, and \code{c} are non-negative constants; furthermore,
#'   reciprocal edges and sampling without replacement are not considered, i.e.,
#'   option \code{rpa_control_reciprocal()} must be set as default,
#'   \code{snode.replace}, \code{tnode.replace} and \code{node.replace} must be
#'   \code{TRUE}. In addition, \code{bag} method only works for unweighted
#'   networks and does not consider multiple edges, i.e.,
#'   \code{rpa_control_edgeweight()} and \code{rpa_control_newedge()} must be
#'   set as default.
#'
#'
#' @return Returns a \code{wdnet} object that includes the following components:
#' \itemize{
#' \item \code{directed}: Logical, whether the network is directed.
#' \item \code{weighted}: Logical, whether the network is weighted.
#' \item \code{edgelist}: A two-column matrix representing the edges.
#' \item \code{edge.attr}: A data frame including edge weights and edge
#' scenarios (0: from initial network; 1: \code{alpha}; 2: \code{beta};
#' 3: \code{gamma}; 4: \code{xi}; 5; \code{rho}; 6: reciprocal edge).
#' \item \code{node.attr}: A data frame including node out- and
#' in-strength, node source and target preference scores (for directed
#' networks), node strength and preference scores (for undirected
#' networks), and node group (if applicable).
#' \item \code{newedge}: The number of new edges at each step, including
#' reciprocal edges.
#' \item \code{control}: An \code{rpacontrol} object that is used to
#' generate the network.
#' }
#'
#' @note The \code{binary} method implements binary search algorithm;
#'   \code{linear} represents linear search algorithm; \code{bag} method
#'   implements the algorithm from Wan et al. (2017); \code{bagx} puts all the
#'   edges into a bag, then samples edges and find the source/target node of the
#'   sampled edge.
#'
#' @references \itemize{\item Wan P, Wang T, Davis RA, Resnick SI (2017).
#'   Fitting the Linear Preferential Attachment Model. Electronic Journal of
#'   Statistics, 11(2), 3738–3780.}
#'
#' @export
#'
#' @examples
#' # Control edge scenario and edge weight through rpa_control_scenario()
#' # and rpa_control_edgeweight(), respectively,
#' # while keeping rpa_control_newedge(),
#' # rpa_control_preference() and rpa_control_reciprocal() as default.
#' set.seed(123)
#' control <- rpa_control_scenario(alpha = 0.5, beta = 0.5) +
#'   rpa_control_edgeweight(
#'     sampler = function(n) rgamma(n, shape = 5, scale = 0.2)
#'   )
#' ret1 <- rpanet(nstep = 1e3, control = control)
#'
#' # In addition, set node groups and probability of creating reciprocal edges.
#' control <- control + rpa_control_reciprocal(
#'   group.prob = c(0.4, 0.6),
#'   recip.prob = matrix(runif(4), ncol = 2)
#' )
#' ret2 <- rpanet(nstep = 1e3, control = control)
#'
#' # Further, set the number of new edges in each step as Poisson(2) + 1 and use
#' # ret2 as a seed network.
#' control <- control + rpa_control_newedge(
#'   sampler = function(n) rpois(n, lambda = 2) + 1
#' )
#' ret3 <- rpanet(nstep = 1e3, initial.network = ret2, control = control)
#'
rpanet <- function(
    nstep, initial.network = list(
      edgelist = matrix(c(1, 2), nrow = 1),
      edgeweight = 1,
      directed = TRUE
    ),
    control,
    method = c("binary", "linear", "bagx", "bag")) {
  method <- match.arg(method)
  stopifnot('"nstep" must be greater than 0.' = nstep > 0)
  if (is.null(initial.network$edgelist)) {
    cat("Assume the initial network has only one edge between nodes 1 and 2.\n")
    initial.network$edgelist <- matrix(c(1, 2), nrow = 1)
  }
  if (is.null(initial.network$directed)) {
    cat("Assume the initial network is directed.\n")
    initial.network$directed <- TRUE
  }
  initial.network <- create_wdnet(
    netwk = initial.network,
    edgelist = initial.network$edgelist,
    edgeweight = initial.network$edgeweight,
    nodegroup = initial.network$nodegroup,
    directed = initial.network$directed
  )
  if (missing(control) || is.null(control)) {
    control <- rpa_control_default()
  } else {
    is_rpacontrol(control)
    control <- rpa_control_default() + control
  }

  nnode <- max(initial.network$edgelist)
  nedge <- nrow(initial.network$edgelist)
  if (is.null(initial.network$node.attr$group)) {
    initial.network$node.attr$group <- rep(1, nnode)
  } else {
    initial.network$node.attr$group <- as.integer(initial.network$node.attr$group)
    stopifnot(
      '"nodegroup" of initial network is not valid.' =
        all(initial.network$node.attr$group > 0) &
          length(initial.network$node.attr$group) == nnode
    )
  }
  if (length(control$reciprocal$group.prob) > 0) {
    stopifnot(
      'Length of "group.prob" is not valid.' =
        max(initial.network$node.attr$group) <= length(control$reciprocal$group.prob)
    )
  }
  if (control$preference$ftype == "customized") {
    control$preference <- compile_pref_func(
      control$preference,
      directed = initial.network$directed
    )
  }

  if (is.null(control$newedge$sampler)) {
    m <- rep(1L, nstep)
  } else {
    tryCatch(do.call(control$newedge$sampler, list(5)), error = function(e) {
      message('Invalid "sampler" for rpa_control_newedge().')
      stop(e)
    })
    m <- do.call(control$newedge$sampler, list(nstep))
  }
  stopifnot(
    'Invalid "sampler" for rpa_control_newedge().' =
      all(length(m) == nstep, m %% 1 == 0, m > 0)
  )

  sum_m <- sum(m)
  sample.recip <- TRUE
  if (is.null(control$reciprocal$recip.prob)) {
    sample.recip <- FALSE
  }

  if (is.null(control$edgeweight$sampler)) {
    w <- rep(1L, sum_m * (1 + sample.recip))
  } else {
    tryCatch(do.call(control$edgeweight$sampler, list(5)), error = function(e) {
      message('Invalid "sampler" for rpa_control_edgeweight().')
      stop(e)
    })
    w <- do.call(control$edgeweight$sampler, list(sum_m * (1 + sample.recip)))
  }
  stopifnot(
    'Invalid "sampler" for rpa_control_edgeweight().' =
      all(length(w) == sum_m * (1 + sample.recip), w > 0)
  )

  if ((!initial.network$directed) &&
    ((!control$newedge$snode.replace) || (!control$newedge$tnode.replace))) {
    cat('"snode.replace" and "tnode.replace" are ignored for undirected networks.')
    control$newedge$snode.replace <- control$tnode.replace <- TRUE
  }
  if (initial.network$directed && (!control$newedge$node.replace)) {
    cat('"node.replace" is ignored for directed networks.')
    control$newedge$node.replace <- TRUE
  }
  if (method == "bag" || method == "bagx") {
    stopifnot(
      '"bag" and "bagx" methods require "default" preference functions.' =
        control$preference$ftype == "default"
    )
    if (initial.network$directed) {
      stopifnot(
        'Source preference must be out-degree plus a non-negative constant for "bag" and "bagx" methods.' =
          all(
            control$preference$sparams[1:2] == 1,
            control$preference$sparams[3:4] == 0,
            control$preference$sparams[5] >= 0
          )
      )
      stopifnot(
        'Target preference must be in-degree plus a non-negative constant for "bag" and "bagx" methods.' =
          all(
            control$preference$tparams[1:2] == 0,
            control$preference$tparams[3:4] == 1,
            control$preference$tparams[5] >= 0
          )
      )
    } else {
      stopifnot(
        'Preference must be degree plus a non-negative constant for "bag" and "bagx" methods.' =
          control$preference$params[1] == 1 &
            control$preference$params[2] >= 0
      )
    }
    stopifnot(
      '"rpa_control_reciprocal" must set as default for "bag" and "bagx" methods.' =
        identical(control$reciprocal, rpa_control_reciprocal()$reciprocal)
    )
    stopifnot(
      '"beta.loop" must be TRUE for "bag" and "bagx" methods.' =
        control$scenario$beta.loop
    )
    if (method == "bag") {
      stopifnot(
        '"rpa_control_edgeweight" must set as default for "bag" method.' =
          identical(control$edgeweight, rpa_control_edgeweight()$edgeweight)
      )
      stopifnot(
        'Weight of existing edges must be 1 for "bag" method.' =
          all(initial.network$edge.attr$weight == 1)
      )
      stopifnot(
        '"rpa_control_newedge" must set as default for "bag" method.' =
          identical(control$newedge, rpa_control_newedge()$newedge)
      )
    }
    if (method == "bagx") {
      if (initial.network$directed) {
        stopifnot(
          '"snode.replace" must be TRUE for "bagx" method.' =
            control$newedge$snode.replace
        )
        stopifnot(
          '"tnode.replace" must be TRUE for "bagx" method.' =
            control$newedge$tnode.replace
        )
      } else {
        stopifnot(
          '"node.replace" must be TRUE for "bagx" method.' =
            control$newedge$node.replace
        )
      }
    }
    return(rpanet_simple(
      nstep = nstep, initial.network = initial.network,
      control = control,
      m = m, sum_m = sum_m,
      w = w, ex_node = nnode,
      ex_edge = nedge, method = method
    ))
  }
  if ((!control$newedge$node.replace) && control$scenario$beta.loop) {
    control$scenario$beta.loop <- FALSE
    cat('"beta.loop" is set to FALSE since "node.replace" is FALSE.')
  }
  return(rpanet_general(
    nstep = nstep, initial.network = initial.network,
    control = control,
    m = m, sum_m = sum_m,
    w = w, nnode = nnode,
    nedge = nedge, method = method,
    sample.recip = sample.recip
  ))
}
