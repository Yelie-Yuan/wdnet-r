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
NULL

#' Generate a PA network with non-linear preference functions
#'
#' @param nstep Number of steps when generating a network.
#' @param seednetwork A list represents the seed network. If \code{NULL}, \code{seednetwork} will
#'   have one edge from node 1 to node 2 with weight 1. It consists of the
#'   following components: a two column matrix \code{edgelist} represents the
#'   edges; a vector \code{edgeweight} represents the weight of edges; a integer
#'   vector \code{nodegroup} represents the group of each node. \code{nodegroup}
#'   is defined for directed networks, if \code{NULL}, all nodes from the seed
#'   graph are considered from group 1.
#' @param control A list of parameters that controls the PA generation process.
#'   The default value is \code{scenario.control() + edgeweight.control() +
#'   newedge.control() + preference.control() + reciprocal.control()}. By
#'   default, in each step, a new edge of weight 1 is added from a new node
#'   \code{A} to an existing node \code{B} (\code{alpha} scenario), where
#'   $\code{B} is chosen with probability proportional to its in-strength + 1.
#' @param directed Logical, whether to generate directed networks. If
#'   \code{FALSE}, the edge directions are omitted.
#' @param m Integer vector, number of new edges in each step.
#' @param sum_m Integer, summation of \code{m}.
#' @param w Vector, weight of new edges.
#' @param nnode Integer, number of nodes in \code{seednetwork}.
#' @param nedge Integer, number of edges in \code{seednetwork}.
#' @param method Which method to use when generating PA networks: "binary" or
#'   "naive".
#'
#' @return A list with the following components: edgelist, edgeweight, strength
#'   for undirected networks, out- and in-strength for directed networks,
#'   controlling parameters, node group (if applicable) and edge scenario
#'   (1~alpha, 2~beta, 3~gamma, 4~xi, 5~rho, 6~reciprocal). The scenario of
#'   edges from the seed network \code{seednetwork} are denoted as 0.
#'   
rpanet_general <- function(nstep, seednetwork, control, directed,
                           m, sum_m, w,
                           nnode, nedge, method) {
  
  if (! all(control$newedge$snode.replace, control$newedge$tnode.replace, 
          control$newedge$node.replace)) {
    if (all(control$preference$sparams[c(3, 5)] <= 0) |
        all(control$preference$tparams[c(1, 5)] <= 0) )
      stop("Source preference function and target preference 
            function must be strictly positive when unique 
            source/target nodes is required.")
  }
  
  if ((! control$newedge$node.replace) & (! directed)) {
    if (control$preference$param[2] <= 0) {
      stop("Preference function must be strictly positive when unique 
            nodes is required.")
    }
  }
  
  edgeweight <- c(seednetwork$edgeweight, w)
  node_vec_length <- (sum_m + nedge) * 2
  node_vec1 <- node_vec2 <-  scenario <- integer(node_vec_length)
  node_vec1[1:nedge] <- seednetwork$edgelist[, 1] - 1
  node_vec2[1:nedge] <- seednetwork$edgelist[, 2] - 1
  scenario[1:nedge] <- 0
  seed_strength <- nodeStrength_cpp(seednetwork$edgelist[, 1], seednetwork$edgelist[, 2], seednetwork$edgeweight,
                                    nnode, weighted = TRUE)
  if (directed) {
    outstrength <- instrength <- double(node_vec_length)
    outstrength[1:nnode] <- seed_strength$outstrength
    instrength[1:nnode] <- seed_strength$instrength
    source_pref <- target_pref <- double(node_vec_length)
    if (is.null(control$reciprocal$group.prob)) {
      sample_recip <- FALSE
      control$reciprocal$group.prob <- 1
      control$reciprocal$recip.prob <- matrix(0)
    }
    else {
      sample_recip <- TRUE
    }
    nodegroup <- c(seednetwork$nodegroup - 1, integer(node_vec_length))
    if (method == "binary") {
      ret_c <- .C("rpanet_binary_directed_cpp",
                  as.integer(nstep),
                  m = as.integer(m),
                  nnode = as.integer(nnode),
                  nedge = as.integer(nedge),
                  node_vec1 = as.integer(node_vec1),
                  node_vec2 = as.integer(node_vec2),
                  outstrength = as.double(outstrength),
                  instrength = as.double(instrength),
                  as.double(edgeweight),
                  scenario = as.integer(scenario),
                  as.double(control$scenario$alpha),
                  as.double(control$scenario$beta),
                  as.double(control$scenario$gamma),
                  as.double(control$scenario$xi),
                  as.integer(control$scenario$beta.loop),
                  as.integer(control$scenario$source.first),
                  as.integer(! control$newedge$node.replace),
                  as.integer(! control$newedge$snode.replace),
                  as.integer(! control$newedge$tnode.replace),
                  as.double(control$preference$sparams),
                  as.double(control$preference$tparams),
                  as.integer(sample_recip),
                  as.double(control$reciprocal$group.prob),
                  as.double(t(control$reciprocal$recip.prob)),
                  nodegroup = as.integer(nodegroup),
                  as.integer(length(control$reciprocal$group.prob)),
                  source_pref = as.double(source_pref),
                  target_pref = as.double(target_pref),
                  PACKAGE = "wdnet")
    } 
    else {
      ret_c <- .C("rpanet_naive_directed_cpp",
                  as.integer(nstep),
                  m = as.integer(m),
                  nnode = as.integer(nnode),
                  nedge = as.integer(nedge),
                  node_vec1 = as.integer(node_vec1),
                  node_vec2 = as.integer(node_vec2),
                  outstrength = as.double(outstrength),
                  instrength = as.double(instrength),
                  as.double(edgeweight),
                  scenario = as.integer(scenario),
                  as.double(control$scenario$alpha),
                  as.double(control$scenario$beta),
                  as.double(control$scenario$gamma),
                  as.double(control$scenario$xi),
                  as.integer(control$scenario$beta.loop),
                  as.integer(control$scenario$source.first),
                  as.integer(! control$newedge$node.replace),
                  as.integer(! control$newedge$snode.replace),
                  as.integer(! control$newedge$tnode.replace),
                  as.double(control$preference$sparams),
                  as.double(control$preference$tparams),
                  as.integer(sample_recip),
                  as.double(control$reciprocal$group.prob),
                  as.double(t(control$reciprocal$recip.prob)),
                  nodegroup = as.integer(nodegroup),
                  as.integer(length(control$reciprocal$group.prob)),
                  source_pref = as.double(source_pref),
                  target_pref = as.double(target_pref),
                  PACKAGE = "wdnet")
    }
    control$preference$params <- NULL
  }
  else {
    sample_recip <- FALSE
    strength <- double(node_vec_length)
    strength[1:nnode] <- seed_strength$outstrength +
      seed_strength$instrength
    pref <- double(node_vec_length)
    if (method == "binary") {
      ret_c <- .C("rpanet_binary_undirected_cpp",
                  as.integer(nstep),
                  m = as.integer(m),
                  nnode = as.integer(nnode),
                  nedge = as.integer(nedge),
                  node_vec1 = as.integer(node_vec1),
                  node_vec2 = as.integer(node_vec2),
                  strength = as.double(strength),
                  as.double(edgeweight),
                  scenario = as.integer(scenario),
                  as.double(control$scenario$alpha),
                  as.double(control$scenario$beta),
                  as.double(control$scenario$gamma),
                  as.double(control$scenario$xi),
                  as.integer(control$scenario$beta.loop),
                  as.integer(! control$newedge$node.replace),
                  as.double(control$param),
                  pref = as.double(pref),
                  PACKAGE = "wdnet")
    }
    else {
      ret_c <- .C("rpanet_naive_undirected_cpp",
                  as.integer(nstep),
                  m = as.integer(m),
                  nnode = as.integer(nnode),
                  nedge = as.integer(nedge),
                  node_vec1 = as.integer(node_vec1),
                  node_vec2 = as.integer(node_vec2),
                  strength = as.double(strength),
                  as.double(edgeweight),
                  scenario = as.integer(scenario),
                  as.double(control$scenario$alpha),
                  as.double(control$scenario$beta),
                  as.double(control$scenario$gamma),
                  as.double(control$scenario$xi),
                  as.integer(control$scenario$beta.loop),
                  as.integer(! control$newedge$node.replace),
                  as.double(control$param),
                  pref = as.double(pref),
                  PACKAGE = "wdnet")
    }
    control$preference$sparams <- control$preference$tparams <- NULL
  }
  nnode <- ret_c$nnode
  nedge <- ret_c$nedge
  node_vec1 <- ret_c$node_vec1[1:nedge] + 1
  node_vec2 <- ret_c$node_vec2[1:nedge] + 1
  scenario <- ret_c$scenario[1:nedge]
  # scenario <- dplyr::recode(scenario, "0" = "NA", "1" = "alpha",
  #    "2" = "beta",
  #    "3" = "gamma", "4" = "xi", "5" = "rho", "6" = "reciprocal")
  m <- ret_c$m
  edgeweight <- edgeweight[1:nedge]
  edgelist <- cbind(node_vec1, node_vec2)
  colnames(edgelist) <- NULL
  ret <- list("edgelist" = edgelist,
              "edgeweight" = edgeweight,
              "scenario" = scenario, 
              "newedge" = m,
              "control" = control,
              "seednetwork" = seednetwork[c("edgelist", "edgeweight", "nodegroup")])
  if (directed) {
    ret$outstrength <- ret_c$outstrength[1:nnode]
    ret$instrength <- ret_c$instrength[1:nnode]
    # ret$source_pref <- ret_c$source_pref[1:nnode]
    # ret$target_pref <- ret_c$target_pref[1:nnode]
    # if (ret$control$newedge$node.replace) {
    #   ret$control$newedge$snode.replace <- ret$control$newedge$tnode.replace <- NULL
    # }
  }
  else {
    ret$strength <- ret_c$strength[1:nnode]
    # ret$pref <- ret_c$pref[1:nnode]
    ret$control$newedge$snode.replace <- ret$control$newedge$tnode.replace <- NULL
  }
  if (sample_recip) {
    ret$nodegroup <- ret_c$nodegroup[1:nnode] + 1
  }
  else {
    ret$control$reciprocal$group.prob <- ret$control$reciprocal$recip.prob <- NULL
    ret$seednetwork$nodegroup <- NULL
  }
  return(ret)
}
