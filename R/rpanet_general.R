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
#' @param alpha The probability of adding an edge from a new node to an
#'   existing node.
#' @param beta The probability of adding an edge between existing nodes.
#' @param gamma The probability of adding an edge from an existing node to a new
#'   node.
#' @param xi The probability of adding an edge between two new nodes.
#' @param rho The probability of introducing a new node with a self looped edge.
#' @param beta_loop Logical, wheter self-loops are allowed in beta scenario.
#' @param w_dist Dsitribution function or a constant for edge weights. Default
#'   value is 0.
#' @param w_par Parameters passed on to w_dist.
#' @param w_const A constant add to w_dist. Weight of new edges follow
#'   distribution w_dist(w_par) + w_const. Default value is 1.
#' @param m_dist Distribution function or a constant for number of new edges per
#'   step. The default value is 0.
#' @param m_par Parameters passed on to m_dist.
#' @param m_const A constant add to m_dist. The number of newly added edges per
#'   step then follows m_dist(m_par) + m_const. The default value is 1.
#' @param m_unique Logical, whether the sampled nodes in the same step should be
#'   different. Defined for undirected and directed networks. For directed
#'   networks, when \code{m_unique} is \code{TRUE}, all the sampled source and
#'   target nodes in the same step are different.
#' @param m_source_unique Logical, whether the sampled source nodes in the same
#'   step should be different. Defined for directed networks.
#' @param m_target_unique Logical, whether the sampled target nodes in the same
#'   step should be different. Defined for directed networks.
#' @param param Parameters of the preference function for undirected networks.
#'   Probability of choosing an exising node is proportional to
#'   \code{strength^param[1] + param[2]}.
#' @param source_param Parameters of the preference function for directed
#'   networks. Probability of choosing an exising node as the source node is
#'   proportional to \code{source_param[1] * out-strength^source_param[2] +
#'   source_param[3] * in-strength^source_param[4] + source_param[5]}.
#' @param target_param Parameters of the preference function for directed
#'   networks. Probability of choosing an exising node as the source node is
#'   proportional to \code{target_param[1] * out-strength^target_param[2] +
#'   target_param[3] * in-strength^target_param[4] + target_param[5]}.
#' @param group_dist The distribution of node groups, its length must equal to
#'   the number of rows of \code{recip_matrix}. If \code{NA}, all groups will
#'   have the same probability. \code{group_dist} and \code{recip_matrix} are
#'   defined for directed networks.
#' @param recip_matrix A square matrix giving the probability of adding a
#'   reciprocal edge after a new edge is introduced. Element \code{p_{ij}}
#'   represents the probability of adding a reciprocal edge from node \code{A},
#'   which belongs to group \code{i}, to node \code{B}, which belongs to group
#'   \code{j}, immediately after a directed edge from \code{B} to \code{A} is
#'   introduced. Note that reciprocal edges are not considered when rho-scenario
#'   edges are introduced.
#'
#' @return List of parameters.
#' @export

general.control  <- function(alpha = 0.5, beta = 0.5, gamma = 0,
    xi = 0, rho = 0, beta_loop = TRUE,
    w_dist = 0, w_par = list(), w_const = 1,
    m_dist = 0, m_par = list(), m_const = 1,
    m_unique = FALSE, m_source_unique = FALSE,
    m_target_unique = FALSE,
    param = c(1, 1),
    source_param = c(1, 1, 0, 0, 1),
    target_param = c(0, 0, 1, 1, 1),
    group_dist = NA,
    recip_matrix = NA) {
    stopifnot("alpha + beta + bamma + xi + rho must equal to 1." =
        round(alpha + beta + gamma + xi + rho, 10) == 1)
    if (! is.na(group_dist[1])) {
        if (! is.na(recip_matrix[1])) {
            stopifnot("'recip_matrix' or 'group_dist' not valid." = 
                is.matrix(recip_matrix) &
                length(group_dist) == nrow(recip_matrix) & 
                nrow(recip_matrix) == ncol(recip_matrix))
            stopifnot("'recip_matrix' not valid." = 
                all(recip_matrix >= 0) & 
                all(recip_matrix <= 1))
            stopifnot("'group_dist' not valid." = 
                round(sum(group_dist), 10) == 1 & 
                all(group_dist >= 0))
        }
        else {
            stop("'recip_matrix' can not be NA when 'group_dist' is provided.")
        }
    }
    if (any(m_unique, m_source_unique, m_target_unique)) {
        if (all(source_param[c(3, 5)] == 0)) {
            stop("Zero out-strength (out-degree) nodes have no source
                preference, please use a different 'source_param'.")
        }
        if (all(target_param[c(1, 5)] == 0)) {
            stop("Zero in-strength (in-degree) nodes have no target 
                preference, please use a different 'target_param'.")
        }
    }
    if (m_unique) {
      m_source_unique <- m_target_unique <- FALSE
    }
    list(alpha = alpha, beta = beta, gamma = gamma, xi = xi, rho = rho,
        beta_loop = beta_loop,
        w_dist = w_dist, w_par = w_par, w_const = w_const,
        m_dist = m_dist, m_par = m_par, m_const = m_const,
        m_unique = m_unique, m_source_unique = m_source_unique, 
        m_target_unique = m_target_unique,
        param = param,
        source_param = source_param,
        target_param = target_param,
        group_dist = group_dist,
        recip_matrix = recip_matrix)
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
#'   \code{FALSE}, the edge directions are omitted.
#' @param node_group A integer vector presents the group of the nodes from the
#'   seed graph. Only defined for directed networks. If \code{NA}, all the nodes
#'   from the seed graph are labeled as group 1.
#'
#' @return A list with the following components: edgelist, edgeweight, strength
#'   for undirected networks, out- and in-strength for directed networks,
#'   control parameters, node group (if applicable) and edge scenario (1~alpha,
#'   2~beta, 3~gamma, 4~xi, 5~rho, 6~reciprocal). The edges from the seed graph
#'   are denoted as 0.
#' @export
#'
#' @examples
#' set.seed(123)
#' net <- rpanet_general(nstep = 1e5,
#'   control = general.control(alpha = 0.4, beta = 0, gamma = 0.6, param = c(1, 2)),
#'   directed = FALSE)
#' net <- rpanet_general(nstep = 1e5, edgelist = matrix(c(1:8), ncol = 2),
#'   control = general.control(w_dist = stats::runif,
#'     w_par = list(min = 1, max = 10), w_const = 0,
#'     source_param = c(1, 1, 0.1, 1, 2),
#'     target_param = c(0.1, 1, 1, 1, 2)))

rpanet_general <- function(nstep = 10^3, edgelist = matrix(c(1, 2), ncol = 2),
                           edgeweight = NA,
                           node_group = NA,
                           control = general.control(),
                           directed = TRUE) {
    stopifnot("nstep must be greater than 0." = nstep > 0)
    temp <- c(edgelist)
    nnode <- max(temp)
    stopifnot("Nodes' index should be consecutive numbers start from 1." =
              sum(! duplicated(temp)) == nnode)
    if (is.na(node_group[1])) {
        node_group <- rep(0, nnode)
    }
    else {
        stopifnot("Value/length of 'node_group' is not valid." =
            all(is.integer(node_group)) & length(node_group) == nnode)
    }
    nedge <- nrow(edgelist)
    if (is.na(edgeweight[1])) edgeweight[1:nedge] <- 1
    stopifnot(length(edgeweight) == nedge)
    if (! is.numeric(control$m_dist)) {
        m <- do.call(control$m_dist, c(nstep, control$m_par)) + control$m_const
    } else {
        m <- rep(control$m_dist + control$m_const, nstep)
    }
    stopifnot("Number of new edges per step must be positive integers." =
              m %% 1 == 0)
    stopifnot("Number of new edges per step must be positive integers." =
              m > 0)
    sum_m <- sum(m)
    if (! is.numeric(control$w_dist)) {
        w <- do.call(control$w_dist, c(sum_m * 2, control$w_par)) + control$w_const
    } else {
        w <- rep(control$w_dist + control$w_const, sum_m * 2)
    }
    stopifnot("Edgeweight must be greater than 0." = w > 0)
    edgeweight <- c(edgeweight, w)
    node_vec_length <- (sum_m + nedge) * 2
    node_vec1 <- node_vec2 <-  scenario <- integer(node_vec_length)
    node_vec1[1:nedge] <- edgelist[, 1] - 1
    node_vec2[1:nedge] <- edgelist[, 2] - 1
    scenario[1:nedge] <- 0
    seed_strength <- nodeStrength_cpp(edgelist[, 1], edgelist[, 2], edgeweight,
        nnode, weighted = TRUE)
    if (directed) {
        outstrength <- instrength <- double(node_vec_length)
        outstrength[1:nnode] <- seed_strength$outstrength
        instrength[1:nnode] <- seed_strength$instrength
        source_pref <- target_pref <- double(node_vec_length)
        if (is.na(control$group_dist[1])) {
            sample_recip <- FALSE
            control$group_dist <- 1
            control$recip_matrix <- matrix(0)
        }
        else {
            sample_recip <- TRUE
        }
        node_group <- c(node_group, integer(node_vec_length))
        ret_c <- .C("rpanet_general_directed_cpp",
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
            as.double(control$alpha),
            as.double(control$beta),
            as.double(control$gamma),
            as.double(control$xi),
            as.integer(control$beta_loop),
            as.integer(control$m_unique),
            as.integer(control$m_source_unique),
            as.integer(control$m_target_unique),
            as.double(control$source_param),
            as.double(control$target_param),
            as.integer(sample_recip),
            as.double(control$group_dist),
            as.double(t(control$recip_matrix)),
            node_group = as.integer(node_group),
            as.integer(length(control$group_dist)),
            source_pref = as.double(source_pref),
            target_pref = as.double(target_pref),
            PACKAGE = "wdnet")
        control$param <- NULL
    }
    else {
        sample_recip <- FALSE
        strength <- double(node_vec_length)
        strength[1:nnode] <- seed_strength$outstrength + seed_strength$instrength
        pref <- double(node_vec_length)
        ret_c <- .C("rpanet_general_undirected_cpp",
            as.integer(nstep),
            m = as.integer(m),
            nnode = as.integer(nnode),
            nedge = as.integer(nedge),
            node_vec1 = as.integer(node_vec1),
            node_vec2 = as.integer(node_vec2),
            strength = as.double(strength),
            as.double(edgeweight),
            scenario = as.integer(scenario),
            as.double(control$alpha),
            as.double(control$beta),
            as.double(control$gamma),
            as.double(control$xi),
            as.integer(control$beta_loop),
            as.integer(control$m_unique),
            as.double(control$param),
            pref = as.double(pref),
            PACKAGE = "wdnet")
        control$source_param <- control$target_param <- NULL
    }
    nnode <- ret_c$nnode
    nedge <- ret_c$nedge
    node_vec1 <- ret_c$node_vec1[1:nedge] + 1
    node_vec2 <- ret_c$node_vec2[1:nedge] + 1
    scenario <- ret_c$scenario[1:nedge]
    # scenario <- dplyr::recode(scenario, "0" = "NA", "1" = "alpha", "2" = "beta",
    #         "3" = "gamma", "4" = "xi", "5" = "rho", "6" = "reciprocal")
    m <- ret_c$m
    edgeweight <- edgeweight[1:nedge]
    edgelist <- cbind(node_vec1, node_vec2)
    colnames(edgelist) <- NULL
    ret <- list("edgelist" = edgelist,
        "edgeweight" = edgeweight,
        "scenario" = scenario,
        "control" = control,
        "m" = m)
    if (directed) {
        ret$outstrength <- ret_c$outstrength[1:nnode]
        ret$instrength <- ret_c$instrength[1:nnode]
        # ret$source_pref <- ret_c$source_pref[1:nnode]
        # ret$target_pref <- ret_c$target_pref[1:nnode]
        if (ret$control$m_unique) {
            ret$control$m_source_unique <- ret$control$m_target_unique <- NULL
        }
    }
    else {
        ret$strength <- ret_c$strength[1:nnode]
        # ret$pref <- ret_c$pref[1:nnode]
        ret$control$m_source_unique <- ret$control$m_target_unique <- NULL
    }
    if (sample_recip) {
        ret$node_group <- ret_c$node_group[1:nnode] + 1
    }
    else {
        ret$control$group_dist <- ret$control$recip_matrix <- NULL
    }
    return(ret)
}
