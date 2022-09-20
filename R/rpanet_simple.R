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

#' @importFrom stats runif
NULL

#' Generate a PA network with linear preference functions.
#'
#' Source preference function must be out-degree (out-strength) plus a
#' nonnegative constant; target preference function must be in-degree
#' (in-strength) plus a nonnegative constant.
#'
#' @param nstep Number of steps when generating a network.
#' @param seednetwork A list represents the seed network. If \code{NULL},
#'   \code{seednetwork} will have one edge from node 1 to node 2 with weight 1.
#'   It consists of the following components: a two column matrix
#'   \code{edgelist} represents the edges; a vector \code{edgeweight} represents
#'   the weight of edges; a integer vector \code{nodegroup} represents the group
#'   of each node. \code{nodegroup} is defined for directed networks, if
#'   \code{NULL}, all nodes from the seed graph are considered from group 1.
#' @param control A list of parameters to be used when generate network.
#' @param directed Logical, whether to generate directed networks. When FALSE,
#'   the edge directions are omitted.
#' @param m Integer vector, number of new edges in each step.
#' @param sum_m Integer, summation of \code{m}.
#' @param w Vector, weight of new edges.
#' @param ex_node Integer, number of nodes in \code{seednetwork}.
#' @param ex_edge Integer, number of edges in \code{seednetwork}.
#' @param method Which method to use, \code{nodelist} or \code{edgesampler}.
#'
#' @return A list with the following components: edgelist, edgeweight, out- and
#'   in-strength, number of edges per step (m), scenario of each new edge
#'   (1~alpha, 2~beta, 3~gamma, 4~xi, 5~rho). The edges in the seed graph are
#'   denoted as scenario 0.
#'   
#' @keywords internal
#' 
rpanet_simple <- function(nstep, seednetwork, control, directed,
                          m, sum_m, w, ex_node, ex_edge, method) {
  delta <- control$preference$params[2]
  delta_out <- control$preference$sparams[5]
  delta_in <- control$preference$tparams[5]
  ex_weight <- sum(seednetwork$edgeweight)
  
  edgeweight <- c(seednetwork$edgeweight, w)
  scenario <- sample(1:5, size = sum_m, replace = TRUE,
                     prob = c(control$scenario$alpha, control$scenario$beta,
                              control$scenario$gamma, control$scenario$xi,
                              control$scenario$rho))
  if (! directed) {
    delta_out <- delta_in <- delta / 2
  }
  if (method == "nodelist") {
    # stopifnot(all(edgeweight == 1) & all(m == 1))
    snode <- c(seednetwork$edgelist[, 1], rep(0, sum_m))
    tnode <- c(seednetwork$edgelist[, 2], rep(0, sum_m))
    ret <- rpanet_nodelist_cpp(snode, tnode,
                               scenario,
                               ex_node, ex_edge,
                               delta_out, delta_in,
                               directed)
    snode <- ret$snode
    tnode <- ret$tnode
    nnode <- ret$nnode
  }
  else {
    scenario1 <- scenario == 1
    scenario4 <- scenario == 4
    
    no_new_start <- !((scenario > 3) | scenario1)
    no_new_end <- scenario < 3
    total_node <- tnode <- cumsum(c((scenario != 2) + scenario4)) +
      ex_node
    snode <- total_node - scenario4
    tnode[no_new_end] <- 0
    snode[no_new_start] <- 0
    nnode <- total_node[length(total_node)]
    
    weight_intv <- cumsum(c(0, edgeweight))
    temp_m <- cumsum(m[-nstep])
    temp <- c(ex_weight, weight_intv[temp_m + ex_edge + 1])
    total_weight <- rep(temp, m)
    temp <- c(ex_node, total_node[temp_m])
    rm(temp_m)
    total_node <- rep(temp, m)
    rm(temp)
    
    rand_out <- runif(sum(no_new_start)) *
      (total_weight + delta_out * total_node)[no_new_start]
    rand_in <- runif(sum(no_new_end)) *
      (total_weight + delta_in * total_node)[no_new_end]
    temp_out <- rand_out <= total_weight[no_new_start]
    if (! all(temp_out)) {
      snode[no_new_start][! temp_out] <- sample_node_cpp(
        total_node[no_new_start][! temp_out])
    }
    temp_in <- rand_in <= total_weight[no_new_end]
    if (! all(temp_in)) {
      tnode[no_new_end][! temp_in] <- sample_node_cpp(
        total_node[no_new_end][!temp_in])
    }
    
    snode <- c(seednetwork$edgelist[, 1], snode)
    tnode <- c(seednetwork$edgelist[, 2], tnode)
    start_edge <- findInterval(
      rand_out[temp_out], weight_intv, left.open = TRUE)
    end_edge <- findInterval(rand_in[temp_in], weight_intv, left.open = TRUE)
    if (directed) {
      snode <- find_node_cpp(snode, start_edge)
      tnode <- find_node_cpp(tnode, end_edge)
    }
    else {
      ret <- find_node_undirected_cpp(
        snode, tnode, start_edge, end_edge)
      snode <- ret$node1
      tnode <- ret$node2
    }
  }
  edgelist <- cbind(snode, tnode)
  strength <- node_strength_cpp(snode, tnode,
                               edgeweight, nnode, weighted = TRUE)
  colnames(edgelist) <- NULL
  ret <- list("edgelist" = edgelist,
              "edgeweight" = edgeweight,
              "scenario" = c(rep(0, ex_edge), scenario),
              "newedge" = m,
              "control" = control,
              "seednetwork" = seednetwork[c("edgelist", "edgeweight")], 
              "directed" = directed)
  if (directed) {
    ret$outstrength <- c(strength$outstrength)
    ret$instrength <- c(strength$instrength)
    ret$control$preference$params <- NULL
  }
  else {
    ret$strength <- c(strength$outstrength) + c(strength$instrength)
    ret$control$newedge$snode.replace <- ret$control$newedge$tnode.replace <- NULL
    ret$control$preference$sparams <- ret$control$preference$tparams <- NULL
  }
  return(ret)
}
