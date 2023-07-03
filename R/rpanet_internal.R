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
#' @importFrom stats runif
NULL

#' Internal functions for generating PA networks
#'
#' These functions generate a PA network with linear (\code{rpanet_simple}) or
#' non-linear (\code{rpanet_general}) preference functions
#'
#' @param nstep Number of steps when generating a network.
#' @param initial.network A \code{wdnet} object or a list that represents the
#'   initial network. By default, \code{initial.network} has one directed edge from node 1
#'   to node 2 with weight 1. It may have the following components: a two-column
#'   matrix \code{edgelist} representing the edges; a vector \code{edgeweight}
#'   representing the weight of edges; a logical argument \code{directed} indicating
#'   whether the initial network is directed;
#'   an integer vector \code{nodegroup}
#'   representing the group of nodes. \code{nodegroup} is defined for directed
#'   networks, if \code{NULL}, all nodes from the seed network are considered
#'   from group 1.
#' @param control A list of parameters that controls the PA generation process.
#'   The default value is \code{wdnet:::rpa_control_default()}. By default, in
#'   each step, a new edge of weight 1 is added from a new node \code{A} to an
#'   existing node \code{B} (\code{alpha} scenario), where $\code{B} is chosen
#'   with probability proportional to its in-strength + 1.
#' @param m Integer vector, number of new edges in each step.
#' @param sum_m Integer, summation of \code{m}.
#' @param w Vector, weight of new edges.
#' @param nnode Integer, number of nodes in \code{initial.network}.
#' @param nedge Integer, number of edges in \code{initial.network}.
#' @param method Which method to use when generating PA networks: "binary" or
#'   "linear".
#' @param sample.recip Whether reciprocal edges will be added.
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
#' @rdname rpanet.internal
#' @keywords internal
#' 
rpanet_general <- function(
    nstep, initial.network, control,
    m, sum_m, w,
    nnode, nedge, method, sample.recip) {
  edgeweight <- c(initial.network$edge.attr$weight, w)
  node_vec_length <- sum_m * 2 + max(nedge, nnode)
  node_vec1 <- integer(node_vec_length)
  node_vec2 <- integer(node_vec_length)
  scenario <- integer(node_vec_length)
  node_vec1[1:nedge] <- initial.network$edgelist[, 1]
  node_vec2[1:nedge] <- initial.network$edgelist[, 2]
  scenario[1:nedge] <- 0

  control$preference$ftype.temp <- ifelse(
    control$preference$ftype == "default",
    1, 2
  )
  if (initial.network$directed) {
    outs <- double(node_vec_length)
    ins <- double(node_vec_length)
    outs[1:nnode] <- initial.network$node.attr$outs
    ins[1:nnode] <- initial.network$node.attr$ins
    spref <- double(node_vec_length)
    tpref <- double(node_vec_length)
    if (!sample.recip) {
      control$reciprocal$group.prob <- 1
      control$reciprocal$recip.prob <- matrix(0)
    }
    nodegroup <- integer(node_vec_length)
    nodegroup[1:nnode] <- initial.network$node.attr$group
    if (method == "binary") {
      ret_c <- rpanet_binary_directed(
        nstep,
        m,
        nnode,
        nedge,
        node_vec1,
        node_vec2,
        outs,
        ins,
        edgeweight,
        scenario,
        sample.recip,
        nodegroup,
        spref,
        tpref,
        control
      )
    } else {
      ret_c <- rpanet_linear_directed_cpp(
        nstep,
        m,
        nnode,
        nedge,
        node_vec1,
        node_vec2,
        outs,
        ins,
        edgeweight,
        scenario,
        sample.recip,
        nodegroup,
        spref,
        tpref,
        control
      )
    }
  } else {
    sample.recip <- FALSE
    s <- double(node_vec_length)
    s[1:nnode] <- initial.network$node.attr$s
    pref <- double(node_vec_length)
    if (method == "binary") {
      ret_c <- rpanet_binary_undirected_cpp(
        nstep,
        m,
        nnode,
        nedge,
        node_vec1,
        node_vec2,
        s,
        edgeweight,
        scenario,
        pref,
        control
      )
    } else {
      ret_c <- rpanet_linear_undirected_cpp(
        nstep,
        m,
        nnode,
        nedge,
        node_vec1,
        node_vec2,
        s,
        edgeweight,
        scenario,
        pref,
        control
      )
    }
  }
  control$preference$ftype.temp <- NULL
  control$preference$spref.pointer <- NULL
  control$preference$tpref.pointer <- NULL
  control$preference$pref.pointer <- NULL
  nnode <- ret_c$nnode
  nedge <- ret_c$nedge
  ret <- structure(
    list(
      "edgelist" = cbind(
        ret_c$node_vec1[1:nedge],
        ret_c$node_vec2[1:nedge]
      ),
      "newedge" = ret_c$m,
      "control" = control,
      "directed" = initial.network$directed,
      "edge.attr" = data.frame(
        "weight" = edgeweight[1:nedge],
        "scenario" = ret_c$scenario[1:nedge]
      )
    ),
    class = "wdnet"
  )
  ret$weighted <- any(ret$edge.attr$weight != 1)
  if (initial.network$directed) {
    ret$node.attr <- data.frame(
      "outs" = ret_c$outs[1:nnode],
      "ins" = ret_c$ins[1:nnode],
      "spref" = ret_c$spref[1:nnode],
      "tpref" = ret_c$tpref[1:nnode]
    )
  } else {
    ret$node.attr <- data.frame(
      "s" = ret_c$s[1:nnode],
      "pref" = ret_c$pref[1:nnode]
    )
  }
  if (sample.recip) {
    ret$node.attr$group <- ret_c$nodegroup[1:nnode]
  } else {
    ret$control$reciprocal$group.prob <- NULL
    ret$control$reciprocal$recip.prob <- NULL
  }
  is_wdnet(ret)
  return(ret)
}


#' @rdname rpanet.internal
#' @keywords internal
#'
rpanet_simple <- function(
    nstep, initial.network, control,
    m, sum_m, w, ex_node, ex_edge, method) {
  delta <- control$preference$params[2]
  delta_out <- control$preference$sparams[5]
  delta_in <- control$preference$tparams[5]
  ex_weight <- sum(initial.network$edge.attr$weight)

  edgeweight <- c(initial.network$edge.attr$weight, w)
  scenario <- sample(1:5,
    size = sum_m, replace = TRUE,
    prob = c(
      control$scenario$alpha, control$scenario$beta,
      control$scenario$gamma, control$scenario$xi,
      control$scenario$rho
    )
  )
  if (!initial.network$directed) {
    delta_out <- delta_in <- delta / 2
  }
  if (method == "bag") {
    snode <- c(initial.network$edgelist[, 1], rep(0, sum_m))
    tnode <- c(initial.network$edgelist[, 2], rep(0, sum_m))
    ret <- rpanet_bag_cpp(
      snode, tnode,
      scenario,
      ex_node, ex_edge,
      delta_out, delta_in,
      initial.network$directed
    )
    snode <- ret$snode
    tnode <- ret$tnode
  } else {
    scenario1 <- scenario == 1
    scenario4 <- scenario == 4

    no_new_start <- !((scenario > 3) | scenario1)
    no_new_end <- scenario < 3
    total_node <- tnode <- cumsum(c((scenario != 2) + scenario4)) +
      ex_node
    snode <- total_node - scenario4
    tnode[no_new_end] <- 0
    snode[no_new_start] <- 0

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
    if (!all(temp_out)) {
      snode[no_new_start][!temp_out] <- sample_node_cpp(
        total_node[no_new_start][!temp_out]
      )
    }
    temp_in <- rand_in <= total_weight[no_new_end]
    if (!all(temp_in)) {
      tnode[no_new_end][!temp_in] <- sample_node_cpp(
        total_node[no_new_end][!temp_in]
      )
    }

    snode <- c(initial.network$edgelist[, 1], snode)
    tnode <- c(initial.network$edgelist[, 2], tnode)
    start_edge <- findInterval(
      rand_out[temp_out], weight_intv,
      left.open = TRUE
    )
    end_edge <- findInterval(rand_in[temp_in], weight_intv, left.open = TRUE)
    if (initial.network$directed) {
      snode <- find_node_cpp(snode, start_edge)
      tnode <- find_node_cpp(tnode, end_edge)
    } else {
      ret <- find_node_undirected_cpp(
        snode, tnode, start_edge, end_edge
      )
      snode <- ret$node1
      tnode <- ret$node2
    }
  }
  edgelist <- cbind(snode, tnode)
  ret <- create_wdnet(
    edgelist = edgelist,
    edgeweight = edgeweight,
    newedge = m,
    control = control,
    directed = initial.network$directed,
    weighted = any(edgeweight != 1)
  )
  ret$edge.attr$scenario <- c(rep(0, ex_edge), scenario)
  if (initial.network$directed) {
    ret$node.attr$spref <- ret$node.attr$outs +
      control$preference$sparams[5]
    ret$node.attr$tpref <- ret$node.attr$ins +
      control$preference$tparams[5]
  } else {
    ret$node.attr$pref <- ret$node.attr$s +
      control$preference$params[2]
  }
  # is_wdnet(ret)
  return(ret)
}
