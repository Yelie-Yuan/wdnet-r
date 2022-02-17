##
## wdnet: Weighted directed network
## Copyright (C) 2022  Yelie Yuan, Panpan Zhang, and Jun Yan
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

#' Parameter settings for function rpanet_general.
#'
#' @param alpha The probability of adding an edge from the new node to an
#'   existing node.
#' @param beta The probability of adding an edge between existing nodes.
#' @param gamma The probability of adding an edge from an existing node to a new
#'   node.
#' @param xi The probability of adding an edge between two new nodes.
#' @param rho The probability of introducing a new node with a self looped edge.
#' @param wdist Dsitribution function or a constant for edge weights.
#'   Default value is 0.
#' @param wpar Additional parameters passed on to wdist.
#' @param wconst A constant add to wdist. Weight of new edges follow
#'   distribution wdist(wpar) + wconst. Default value is 1.
#' @param param Parameters of the preference function for undirected networks.
#'   Probability of choosing an exising node is proportional to
#'   \code{strength^param[1] + param[2]}.
#' @param source_param Parameters of the preference function for directed
#'   networks. Probability of choosing an exising node as the source node
#'   is proportional to \code{source_param[1] *
#'   out-strength^source_param[2] + source_param[3] *
#'   in-strength^source_param[4] + source_param[5]}.
#' @param target_param Parameters of the preference function for directed
#'   networks. Probability of choosing an exising node as the source node
#'   is proportional to \code{target_param[1] *
#'   in-strength^target_param[2] + target_param[3] *
#'   out-strength^target_param[4] + target_param[5]}.
#'
#' @return List of parameters.
#' @export

general.control  <- function(alpha = 0.5, beta = 0.5, gamma = 0,
    xi = 0, rho = 0,
    wdist = 0, wpar = list(), wconst = 1,
    param = c(1, 1),
    source_param = c(1, 1, 0, 0, 1),
    target_param = c(1, 1, 0, 0, 1)) {
  ## set default value here
  list(alpha = alpha, beta = beta, gamma = gamma, xi = xi, rho = rho,
       wdist = wdist, wpar = wpar, wconst = wconst,
       param = param,
       source_param = source_param,
       target_param = target_param)
}

#' Generate a preferential attachment network using augument tree method.
#'
#' @param edgelist A two column matrix represents the seed graph.
#' @param edgeweight A vector represents the weight of edges of the seed graph.
#'   Its length equals the number of edges of the seed graph. If NA, all the
#'   edges of the seed graph have weight 1.
#' @param nstep Number of steps when generating a network.
#' @param control A list of parameters that controls the process.
#' @param directed Logical, whether to generate directed networks. When 
#' \code{FALSE}, the edge directions are omitted.
#'
#' @return A list with the following components: edgelist, edgeweight, strength
#'   for undirected networks, out- and in-strength for directed networks,
#'   scenario of each new edge
#'   (1~alpha, 2~beta, 3~gamma, 4~xi, 5~rho). The edges in the seed graph
#'   are denoted as scenario 0.
#' @export
#'
#' @examples
#' net <- rpanet_general(nstep = 100,
#'   control = general.control(alpha = 0.4, beta = 0, gamma = 0.6))
#' net <- rpanet_general(edgelist = matrix(c(1:8), ncol = 2), nstep = 10^5,
#'   control = general.control(wdist = stats::runif,
#'   wpar = list(min = 1, max = 10), wconst = 0))

rpanet_general <- function(nstep = 10^3, edgelist = matrix(c(1, 2), ncol = 2),
                           edgeweight = NA,
                           control = general.control(),
                           directed = TRUE) {
    stopifnot("nstep must be greater than 0." = nstep > 0)
    stopifnot("alpha + beta + bamma + xi + rho must less or equal to 1." =
              control$rho >= 0)
    temp <- c(edgelist)
    nnode <- max(temp)
    stopifnot("Nodes' index should be consecutive numbers start from 1." =
              sum(! duplicated(temp)) == nnode)
    nedge <- nrow(edgelist)
    if (is.na(edgeweight[1])) edgeweight[1:nedge] <- 1
    stopifnot(length(edgeweight) == nedge)
    if (! is.numeric(control$wdist)) {
        w <- do.call(control$wdist, c(nstep, control$wpar)) + control$wconst
    } else {
        w <- rep(control$wdist + control$wconst, nstep)
    }
    stopifnot("Edgeweight must be greater than 0." = w > 0)

    node_vec1 <- node_vec2 <- scenario <- integer(nstep + nedge)
    node_vec1[1:nedge] <- edgelist[, 1] - 1
    node_vec2[1:nedge] <- edgelist[, 2] - 1
    scenario[1:nedge] <- 0

    strength <- nodeStrength_cpp(edgelist[, 1], edgelist[, 2], edgeweight,
        nnode, weighted = TRUE)
    edgeweight <- c(edgeweight, w)
    if (directed) {
        ret_c <- .C("rpanet_directed_general_cpp",
            as.integer(nstep),
            nnode = as.integer(nnode),
            nedge = as.integer(nedge),
            node_vec1 = as.integer(node_vec1),
            node_vec2 = as.integer(node_vec2),
            as.double(strength$outstrength),
            as.double(strength$instrength),
            as.double(edgeweight),
            scenario = as.integer(scenario),
            as.double(control$alpha),
            as.double(control$beta),
            as.double(control$gamma),
            as.double(control$xi),
            as.double(control$source_param),
            as.double(control$target_param),
            PACKAGE = "wdnet")
        control$param <- NULL
    }
    else {
        ret_c <- .C("rpanet_undirected_general_cpp",
            as.integer(nstep),
            nnode = as.integer(nnode),
            nedge = as.integer(nedge),
            node_vec1 = as.integer(node_vec1),
            node_vec2 = as.integer(node_vec2),
            as.double(strength$outstrength + strength$instrength),
            as.double(edgeweight),
            scenario = as.integer(scenario),
            as.double(control$alpha),
            as.double(control$beta),
            as.double(control$gamma),
            as.double(control$xi),
            as.double(control$param),
            PACKAGE = "wdnet")
        control$source_param <- control$target_param <- NULL
    }
    ret_c$node_vec1 <- ret_c$node_vec1 + 1
    ret_c$node_vec2 <- ret_c$node_vec2 + 1
    edgelist <- cbind(ret_c$node_vec1, ret_c$node_vec2)
    colnames(edgelist) <- NULL
    ret <- list("edgelist" = edgelist,
        "edgeweight" = edgeweight,
        "scenario" = ret_c$scenario,
        "control" = control)
    strength <- nodeStrength_cpp(ret_c$node_vec1, ret_c$node_vec2,
        edgeweight, ret_c$nnode, TRUE)
    if (directed) {
        ret$outstrength <- c(strength$outstrength)
        ret$instrength <- c(strength$instrength)
    }
    else {
        ret$strengh <- c(strength$outstrength) + c(strength$instrength)
    }
    return(ret)
}