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
NULL

#' Generate a PA network with non-linear preference functions
#'
#' @param nstep Number of steps when generating a network.
#' @param seednetwork A list represents the seed network. If \code{NULL},
#'   \code{seednetwork} will have one edge from node 1 to node 2 with weight 1.
#'   It consists of the following components: a two column matrix
#'   \code{edgelist} represents the edges; a vector \code{edgeweight} represents
#'   the weight of edges; a integer vector \code{nodegroup} represents the group
#'   of nodes. \code{nodegroup} is defined for directed networks, if
#'   \code{NULL}, all nodes from the seed graph are considered from group 1.
#' @param control A list of parameters that controls the PA generation process.
#'   The default value is \code{rpactl.scenario() + rpactl.edgeweight() +
#'   rpactl.newedge() + rpactl.preference() + rpactl.reciprocal()}. By
#'   default, in each step, a new edge of weight 1 is added from a new node
#'   \code{A} to an existing node \code{B} (\code{alpha} scenario), where
#'   $\code{B} is chosen with probability proportional to its in-strength + 1.
#' @param directed Logical, whether to generate directed networks. If
#'   \code{FALSE}, the edge directions are ignored.
#' @param m Integer vector, number of new edges in each step.
#' @param sum_m Integer, summation of \code{m}.
#' @param w Vector, weight of new edges.
#' @param nnode Integer, number of nodes in \code{seednetwork}.
#' @param nedge Integer, number of edges in \code{seednetwork}.
#' @param method Which method to use when generating PA networks: "binary" or
#'   "naive".
#' @param sample.recip Whether reciprocal edges will be added.
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
#' @keywords internal
#' 
rpanet_general <- function(nstep, seednetwork, control, directed,
                           m, sum_m, w,
                           nnode, nedge, method, sample.recip) {  
  edgeweight <- c(seednetwork$edgeweight, w)
  node_vec_length <- (sum_m + nedge) * 2
  node_vec1 <- integer(node_vec_length)
  node_vec2 <- integer(node_vec_length)
  scenario <- integer(node_vec_length)
  node_vec1[1:nedge] <- seednetwork$edgelist[, 1] - 1
  node_vec2[1:nedge] <- seednetwork$edgelist[, 2] - 1
  scenario[1:nedge] <- 0
  seed_strength <- node_strength_cpp(seednetwork$edgelist[, 1], 
                                     seednetwork$edgelist[, 2], 
                                     seednetwork$edgeweight,
                                     nnode, weighted = TRUE)
  control$preference$ftype.temp <- ifelse(control$preference$ftype == "default", 
                                          yes = 1, no = 2)
  if (directed) {
    outstrength <- double(node_vec_length)
    instrength <- double(node_vec_length)
    outstrength[1:nnode] <- seed_strength$outstrength
    instrength[1:nnode] <- seed_strength$instrength
    source_pref <- double(node_vec_length)
    target_pref <- double(node_vec_length)
    if (! sample.recip) {
      control$reciprocal$group.prob <- 1
      control$reciprocal$recip.prob <- matrix(0)
    }
    nodegroup <- integer(node_vec_length)
    nodegroup[1:nnode] <- seednetwork$nodegroup - 1
    # nodegroup <- c(seednetwork$nodegroup - 1, integer(node_vec_length))
    if (method == "binary") {
      ret_c <- rpanet_binary_directed(nstep,
                                      m,
                                      nnode,
                                      nedge,
                                      node_vec1,
                                      node_vec2,
                                      outstrength,
                                      instrength,
                                      edgeweight,
                                      scenario,
                                      sample.recip,
                                      nodegroup,
                                      source_pref,
                                      target_pref,
                                      control)
    } 
    else {
      ret_c <- rpanet_naive_directed_cpp(nstep,
                                         m,
                                         nnode,
                                         nedge,
                                         node_vec1,
                                         node_vec2,
                                         outstrength,
                                         instrength,
                                         edgeweight,
                                         scenario,
                                         sample.recip,
                                         nodegroup,
                                         source_pref,
                                         target_pref,
                                         control)
    }
  }
  else {
    sample.recip <- FALSE
    strength <- double(node_vec_length)
    strength[1:nnode] <- seed_strength$outstrength +
      seed_strength$instrength
    pref <- double(node_vec_length)
    if (method == "binary") {
      ret_c <- rpanet_binary_undirected_cpp(nstep,
                                            m,
                                            nnode,
                                            nedge,
                                            node_vec1,
                                            node_vec2,
                                            strength,
                                            edgeweight,
                                            scenario,
                                            pref,
                                            control)
    }
    else {
      ret_c <- rpanet_naive_undirected_cpp(nstep,
                                           m,
                                           nnode,
                                           nedge,
                                           node_vec1,
                                           node_vec2,
                                           strength,
                                           edgeweight,
                                           scenario,
                                           pref,
                                           control)
    }
  }
  control$preference$ftype.temp <- NULL
  control$preference$spref.pointer <- NULL
  control$preference$tpref.pointer <- NULL
  control$preference$pref.pointer <- NULL
  if (control$preference$ftype == "customized") {
    control$preference$sparams <- NULL
    control$preference$tparams <- NULL
    control$preference$params <- NULL
  }
  else {
    control$preference$spref <- NULL
    control$preference$tpref <- NULL
    control$preference$pref <- NULL
  }
  if (directed) {
    control$newedge$node.replace <- NULL
    control$preference$params <- NULL
    control$preference$pref <- NULL
  }
  else {
    control$newedge$snode.replace <- NULL
    control$newedge$tnode.replace <- NULL
    control$preference$sparams <- NULL
    control$preference$tparams <- NULL
    control$preference$spref <- NULL
    control$preference$tpref <- NULL
  }
  nnode <- ret_c$nnode
  nedge <- ret_c$nedge
  ret <- list("edgelist" = cbind(ret_c$node_vec1[1:nedge] + 1, 
                                 ret_c$node_vec2[1:nedge] + 1),
              "edgeweight" = edgeweight[1:nedge],
              "scenario" = ret_c$scenario[1:nedge], 
              "newedge" = ret_c$m,
              "control" = control,
              "seednetwork" = seednetwork[c("edgelist", "edgeweight", "nodegroup")], 
              "directed" = directed)
  colnames(ret$edgelist) <- NULL
  if (directed) {
    ret$outstrength <- ret_c$outstrength[1:nnode]
    ret$instrength <- ret_c$instrength[1:nnode]
    ret$spref <- ret_c$source_pref[1:nnode]
    ret$tpref <- ret_c$target_pref[1:nnode]
  }
  else {
    ret$strength <- ret_c$strength[1:nnode]
    ret$pref <- ret_c$pref[1:nnode]
  }
  if (sample.recip) {
    ret$nodegroup <- ret_c$nodegroup[1:nnode] + 1
  }
  else {
    ret$control$reciprocal$group.prob <- NULL
    ret$control$reciprocal$recip.prob <- NULL
    ret$seednetwork$nodegroup <- NULL
  }
  return(ret)
}
