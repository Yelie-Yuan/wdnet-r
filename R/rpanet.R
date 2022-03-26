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

#' @importFrom stats runif rpois
NULL

#' Parameter settings for function rpanet.
#'
#' @param alpha The probability of adding an edge from the new node to an
#'   existing node.
#' @param beta The probability of adding an edge between existing nodes.
#' @param gamma The probability of adding an edge from an existing node to a new
#'   node.
#' @param xi The probability of adding an edge between two new nodes.
#' @param rho The probability of introducing a new node with a self looped edge.
#' @param delta_out A tuning parameter related to growth rate for directed
#'   networks. Probability of choosing an existing node as the source node of
#'   the newly added edge is proportional to node outstrength + outdelta.
#' @param delta_in A tuning parameter related to growth rate for directed
#'   networks. Probability of choosing an existing node as the target node of
#'   the newly added edge is proportional to node instrength + indelta.
#' @param delta A tuning parameter related to growth rate for undirected
#'   networks. Probability of choosing an existing node is proportional to node
#'   strength + delta.
#' @param m_dist Distribution function or a constant for number of new
#'   edges per step. The default value is 0.
#' @param m_par Additional parameters passed on to m_dist.
#' @param m_const A constant add to m_dist. The number of newly added edges per step
#'   then follows m_dist(m_par) + m_const. The default value is 1.
#' @param w_dist Dsitribution function or a constant for edge weights. The
#'   default value is 0.
#' @param w_par Additional parameters passed on to w_dist.
#' @param w_const A constant add to w_dist. Weight of new edges follow
#'   distribution w_dist(w_par) + w_const. Default value is 1.
#'
#' @return List of parameters.
#' @export


panet.control  <- function(alpha = 0.5, beta = 0.5, gamma = 0, xi = 0, rho = 0,
                           delta_out = 0.1, delta_in = 0.1, delta = 0.1,
                           m_dist = 0,
                           m_par = list(), m_const = 1,
                           w_dist = 0,
                           w_par = list(), w_const = 1) {
  ## set default value here
  list(alpha = alpha, beta = beta, gamma = gamma, xi = xi, rho = rho,
       delta_out = delta_out, delta_in = delta_in, delta = delta,
       m_dist = m_dist, m_par = m_par, m_const = m_const,
       w_dist = w_dist, w_par = w_par, w_const = w_const)
}



#' Generate a growing preferential attachment network.
#'
#' @param edgelist A two column matrix represents the seed graph.
#' @param edgeweight A vector represents the weight of edges of the seed graph.
#'   Its length equals the number of edges of the seed graph. If NA, all the
#'   edges of the seed graph have weight 1.
#' @param nstep Number of steps when generating a network.
#' @param control A list of parameters to be used when generate network.
#' @param directed Logical, whether to generate directed networks. When FALSE,
#' the edge directions are omitted.
#'
#' @return A list with the following components: edgelist, edgeweight, out- and
#'   in-strength, number of edges per step (m), scenario of each new edge
#'   (1~alpha, 2~beta, 3~gamma, 4~xi, 5~rho). The edges in the seed graph
#'   are denoted as scenario 0.
#' @export
#'
#' @examples
#' net <- rpanet(nstep = 100,
#'         control = panet.control(alpha = 0.4, beta = 0, gamma = 0.6))
#' net <- rpanet(edgelist = matrix(c(1:8), ncol = 2), nstep = 10^5,
#'       control = panet.control(m_dist = stats::rpois,
#'       m_par = list(lambda = 1), m_const = 1,
#'       w_dist = stats::runif, w_par = list(min = 1, max = 10), w_const = 0))

rpanet <- function(nstep = 10^3, edgelist = matrix(c(1, 2), ncol = 2),
                   edgeweight = NA,
                   control = panet.control(),
                   directed = TRUE) {
    if (is.null(control$delta_out)) {
      control$delta_out <- 0
    }
    if (is.null(control$delta_in)) {
      control$delta_in <- 0
    }
    stopifnot("nstep must be greater than 0." = nstep > 0)
    stopifnot("alpha + beta + bamma + xi + rho must equals to 1." =
              round(control$alpha + control$beta + control$gamma +
                    control$xi + control$rho, 2) == 1)
    temp <- c(edgelist)
    ex_node <- max(temp)
    stopifnot("Nodes' index should start from 1." =
              sum(! duplicated(temp)) == ex_node)
    ex_edge <- nrow(edgelist)
    if (is.na(edgeweight[1])) edgeweight[1:ex_edge] <- 1
    stopifnot(length(edgeweight) == ex_edge)
    ex_weight <- sum(edgeweight)
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
        w <- do.call(control$w_dist, c(sum_m, control$w_par)) + control$w_const
    } else {
        w <- rep(control$w_dist + control$w_const, sum_m)
    }
    stopifnot("Edge weight must be greater than 0." = w > 0)

    edgeweight <- c(edgeweight, w)
    scenario <- sample(1:5, size = sum_m, replace = TRUE,
                       prob = c(control$alpha, control$beta,
                       control$gamma, control$xi,
                       control$rho))
    if (! directed) {
        control$delta_out <- control$delta_in <- control$delta / 2
    }

    if (all(edgeweight == edgeweight[1]) & all(m == 1)) {
        control$delta_out <- control$delta_out / edgeweight[1]
        control$delta_in <- control$delta_in / edgeweight[1]
        start_node <- c(edgelist[, 1], rep(0, sum_m))
        end_node <- c(edgelist[, 2], rep(0, sum_m))
        ret <- rpanet_cpp(start_node, end_node,
                          scenario,
                          ex_node, ex_edge,
                          control$delta_out, control$delta_in,
                          directed)
        start_node <- ret$start_node
        end_node <- ret$end_node
        nnode <- ret$nnode
    }
    else {
        scenario1 <- scenario == 1
        scenario4 <- scenario == 4

        no_new_start <- !((scenario > 3) | scenario1)
        no_new_end <- scenario < 3
        total_node <- end_node <- cumsum(c((scenario != 2) + scenario4)) +
            ex_node
        start_node <- total_node - scenario4
        end_node[no_new_end] <- 0
        start_node[no_new_start] <- 0
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
            (total_weight + control$delta_out * total_node)[no_new_start]
        rand_in <- runif(sum(no_new_end)) *
            (total_weight + control$delta_in * total_node)[no_new_end]
        temp_out <- rand_out <= total_weight[no_new_start]
        if (! all(temp_out)) {
            start_node[no_new_start][! temp_out] <- sampleNode_cpp(
              total_node[no_new_start][! temp_out])
        }
        temp_in <- rand_in <= total_weight[no_new_end]
        if (! all(temp_in)) {
            end_node[no_new_end][! temp_in] <- sampleNode_cpp(
              total_node[no_new_end][!temp_in])
        }

      start_node <- c(edgelist[, 1], start_node)
      end_node <- c(edgelist[, 2], end_node)
      start_edge <- findInterval(
          rand_out[temp_out], weight_intv, left.open = TRUE)
      end_edge <- findInterval(rand_in[temp_in], weight_intv, left.open = TRUE)
      if (directed) {
          start_node <- findNode_cpp(start_node, start_edge)
          end_node <- findNode_cpp(end_node, end_edge)
      }
      else {
          ret <- findNode_undirected_cpp(
            start_node, end_node, start_edge, end_edge)
          start_node <- ret$start_node
          end_node <- ret$end_node
      }
    }
    edgelist <- cbind(start_node, end_node)
    strength <- nodeStrength_cpp(start_node, end_node,
                                 edgeweight, nnode, weighted = TRUE)
    colnames(edgelist) <- NULL
    ret <- list(edgelist = edgelist,
                edgeweight = edgeweight,
                scenario = c(rep(0, ex_edge), scenario),
                m = m)
    if (directed) {
        ret$outstrength <- c(strength$outstrength)
        ret$instrength <- c(strength$instrength)
        control$delta <- NULL
    }
    else {
        ret$strength <- c(strength$outstrength) + c(strength$instrength)
        control$delta_out <- control$delta_in <- NULL
    }
    ret$control <- control
    return(ret)
}