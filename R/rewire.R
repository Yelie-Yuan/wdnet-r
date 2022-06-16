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

#' @importFrom stats cor
#' @importFrom CVXR norm2
NULL

#' Degree preserving rewiring.
#'
#' There are two steps in this algorithm. It first solves for an appropriate
#' \code{eta} using the \code{target_assortcoeff}, \code{FUN}, and
#' \code{control}, then proceeds to the rewiring process and rewire the network
#' towards the structure \code{eta}. If \code{eta} is given, the algorithm will
#' skip the first step. The algorithm is prepared for unweighted networks.
#'
#' @param edgelist A two column matrix, each row represents an edge of the
#'   network.
#' @param directed Logical, whether the network is directed or not.
#' @param target_assortcoeff For directed networks, it is a list represents the
#'   predetermined value or range of assortativity coefficients. For undirected
#'   networks, it is a constant between -1 to 1. It will be ignored if
#'   \code{eta} is provided.
#' @param iteration An integer, number of rewiring iterations. Each iteration
#'   consists of \code{nattempts} rewiring attempts. The assortativity
#'   coefficient(s) of the network will be recorded after each iteration.
#' @param eta An matrix represents the target network structure. If specified,
#'   the \code{target_assortcoeff} will be ignored. For directed networks, the
#'   element at row "i-j" and column "k-l" represents the proportion of directed
#'   edges linking a source node with out-degree i and in-degree j to a target
#'   node with out-degree k and in-degree l. For undirected networks, \code{eta}
#'   is symmetric, the summation of the elements at row "i", column "j" and row
#'   "j", column "i" represents the proportion of edges linking to a node with
#'   degree i and a node with degree j.
#' @param history Logical, whether the rewiring attempts should be recorded and
#'   returned.
#' @param control A list of parameters passed to \code{CVXR::solve()} for
#'   solving an appropriate \code{eta} when \code{target_assortcoeff} is
#'   provided. It will be ignored if \code{eta} is provided.
#' @param FUN A convex function of \code{eta} to be minimized when
#'   \code{target_assortcoeff} is provided. Defaults to 0. It will be ignored if
#'   \code{eta} is provided.
#' @param nattempts An integer, number of rewiring attempts for each iteration.
#'   Default value equals the number of rows of \code{edgelist}.
#'
#' @return Rewired \code{edgelist}; assortativity coefficients after each
#'   iteration; rewiring history (including the index of sampled edges and
#'   rewiring result); results from the solver; solved \code{eta} and its
#'   corresponding assortativity coefficient(s), if applicable.
#'
#' @note Each rewiring attempt samples two rows from \code{edgelist}, for
#'   example Edge1:(v_1, v_2) and Edge2:(v_3, v_4). For directed networks, if
#'   the rewiring attempt is accepted, the sampled edges are replaced as (v_1,
#'   v_4), (v_3, v_2); for undirected networks, we try to rewire the sampled
#'   edges as \{v_1, v_4\}, \{v_3, v_2\} (rewire type 1) or \{v_1, v_3\}, \{v_2,
#'   v_4\} (rewire type 2), each with probability 1/2.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' edgelist <- rpanet(2e4, control = scenario.control(
#'    alpha = 0.4, beta = 0.3, gamma = 0.3))$edgelist
#' target_assortcoeff <- list("outout" = -0.2, "outin" = 0.2)
#' ret1 <- rewire(edgelist, directed = TRUE, 
#'                target_assortcoeff = target_assortcoeff,
#'                iteration = 200)
#' plot(ret1$assortcoeff$Iteration, ret1$assortcoeff$"outout")
#' plot(ret1$assortcoeff$Iteration, ret1$assortcoeff$"outin")
#' 
#' ret2 <- assortcoeff_range(edgelist, directed = TRUE, type = "outin", 
#'             constr = list("outout" = c(-0.1, 0.1), "inout" = 0.1))
#' ret3 <- rewire(edgelist, eta = ret2$lbound$eta, iteration = 100)
#' plot(ret3$assortcoeff$Iteration, ret3$assortcoeff$"outin")
#' assortcoeff(ret3$edgelist)$"outin"
#'
#' edgelist <- rpanet(1e4, control = scenario.control(
#'                    alpha = 0.3, beta = 0.1, gamma = 0.3, xi = 0.3), 
#'                    directed = FALSE)$edgelist
#' ret4 <- rewire(edgelist, directed = FALSE, target_assortcoeff = 0.3,
#'                iteration = 100, FUN = CVXR::norm2, history = TRUE)
#' plot(ret4$assortcoeff$Iteration, ret4$assortcoeff$Value)
#' }
#' 


rewire <- function(edgelist, directed = TRUE,
                   iteration = 10, nattempts = NULL, history = FALSE, 
                   target_assortcoeff = list("outout" = NULL, 
                                             "outin" = NULL,
                                             "inout" = NULL, 
                                             "inin" = NULL),
                   control = solver.control(), FUN = function(x) 0,
                   eta = NULL) {
  solver_result <- NULL
  if (is.null(eta)) {
    if (directed) {
      solver_result <- get_jointdist_directed(edgelist = edgelist,
                                    target_assortcoeff = target_assortcoeff,
                                    FUN = FUN, control = control)
    }
    else {
      solver_result <- get_jointdist_undirected(edgelist = edgelist,
                                      target_assortcoeff = target_assortcoeff,
                                      FUN = FUN, control = control)
    }
    eta <- solver_result$eta
    solver_result$e <- NULL
  }
  if (directed) {
    ret <- rewire_directed(edgelist = edgelist,
                           eta = eta, 
                           iteration = iteration,
                           nattempts = nattempts,
                           rewire_history = history)
  }
  else {
    ret <- rewire_undirected(edgelist = edgelist,
                             eta = eta, 
                             iteration = iteration,
                             nattempts = nattempts,
                             rewire_history = history)
  }
  ret$"solver_result" <- solver_result
  ret
}



#' Degree preserving rewiring towards the target structure \code{eta}.
#'
#' @param edgelist A two column matrix, each row represents a directed edge from
#'   the first column to the second column.
#' @param eta An matrix, target structure eta generated by
#'   \code{wdnet::get_jointdist_directed()}.
#' @param iteration An integer, number of rewiring iterations, each iteration
#'   consists of \code{nattempts} rewiring attempts.
#' @param nattempts An integer, number of rewiring attempts for each iteration.
#'   Default value equals the number of rows of edgelist.
#' @param rewire_history Logical, whether the rewiring history should be
#'   returned.
#'
#' @return Rewired edgelist, degree based assortativity coefficients after each
#'   iteration, rewiring history (including the index of sampled edges and
#'   rewiring result). For each rewiring attempt, two rows are sampled form the
#'   edgelist, for example Edge1:(v_1, v_2) and Edge2:(v_3, v_4), if the
#'   rewiring attempt is accepted, the sampled edges are replaced as (v_1, v_4),
#'   (v_3, v_2).
#'
#' @examples
#' set.seed(1234)
#' edgelist <- rpanet(10000, control = scenario.control(
#'    alpha = 0.3, beta = 0.1, gamma = 0.3, xi = 0.3))$edgelist
#' target_assortcoeff <- list("outout" = -0.1, "outin" = 0.3, "inout" = 0.2, "inin" = 0.1)
#' ret1 <- wdnet:::get_jointdist_directed(edgelist, target_assortcoeff = target_assortcoeff)
#' ret2 <- wdnet:::rewire_directed(edgelist, eta = ret1$eta, iteration = 200)
#' plot(ret2$assortcoeff$Iteration, ret2$assortcoeff$"outin")
#' 
rewire_directed <- function(edgelist, eta, 
                            iteration = 1, nattempts = NULL, 
                            rewire_history = FALSE) {
  if (is.null(nattempts)) nattempts <- nrow(edgelist)
  edgelist <- as.matrix(edgelist)
  temp <- range(edgelist)
  sourceNode <- edgelist[, 1]
  targetNode <- edgelist[, 2]
  stopifnot(temp[1] == 1)
  stopifnot(temp[2] == length(unique(c(edgelist))))
  temp <- nodeStrength_cpp(snode = sourceNode, 
                           tnode = targetNode, 
                           nnode = temp[2], 
                           weight = 1,
                           weighted = FALSE)
  outd <- temp$outstrength
  ind <- temp$instrength
  
  df_s <- data.frame(type = rownames(eta), 
                     index = seq_len(nrow(eta)) - 1)
  df_t <- data.frame(type = colnames(eta), 
                     index = seq_len(ncol(eta)) - 1)
  type_s <- paste0(outd[sourceNode], "-", ind[sourceNode], split = "")
  type_t <- paste0(outd[targetNode], "-", ind[targetNode], split = "")
  
  index_s <- df_s[match(type_s, df_s$type), "index"]
  index_t <- df_t[match(type_t, df_t$type), "index"]
  sourceOut <- outd[sourceNode]
  sourceIn <- ind[sourceNode]
  targetOut <- outd[targetNode]
  targetIn <- ind[targetNode]
  # r_sourceOut <- rank(outd[sourceNode], ties.method = "average")
  # r_sourceIn <- rank(ind[sourceNode], ties.method = "average")
  # r_targetOut <- rank(outd[targetNode], ties.method = "average")
  # r_targetIn <- rank(ind[targetNode], ties.method = "average")
  rm(df_s, df_t, type_s, type_t, temp, outd, ind)
  
  ret <- rewire_directed_cpp(iteration, nattempts, 
                             targetNode - 1, 
                             sourceOut, sourceIn,
                             targetOut, targetIn,
                             #  r_sourceOut, r_sourceIn,
                             #  r_targetOut, r_targetIn,
                             index_s, index_t, 
                             eta, rewire_history)
  rho <- data.frame("Iteration" = c(0:iteration), 
                    "outout" = NA, 
                    "outin" = NA, 
                    "inout" = NA, 
                    "inin" = NA)
  # rankRho <- data.frame("Iteration" = c(0:iteration), 
  #                       "r-out-out" = NA, 
  #                       "r-out-in" = NA, 
  #                       "r-in-out" = NA, 
  #                       "r-in-in" = NA)
  rho[1, 2:5] <- c("outout" = stats::cor(sourceOut, targetOut), 
                   "outin" = stats::cor(sourceOut, targetIn), 
                   "inout" = stats::cor(sourceIn, targetOut),
                   "inin" = stats::cor(sourceIn, targetIn))
  rho[2:(iteration + 1), 2] <- ret$out_out
  rho[2:(iteration + 1), 3] <- ret$out_in
  rho[2:(iteration + 1), 4] <- ret$in_out
  rho[2:(iteration + 1), 5] <- ret$in_in
  # rankRho[1, 2:5] <- c("r-out-out" = stats::cor(r_sourceOut, r_targetOut), 
  #                      "r-out-in" = stats::cor(r_sourceOut, r_targetIn), 
  #                      "r-in-out" = stats::cor(r_sourceIn, r_targetOut),
  #                      "r-in-in" = stats::cor(r_sourceIn, r_targetIn))
  # rankRho[2:(iteration + 1), 2] <- ret$r_out_out
  # rankRho[2:(iteration + 1), 3] <- ret$r_out_in
  # rankRho[2:(iteration + 1), 4] <- ret$r_in_out
  # rankRho[2:(iteration + 1), 5] <- ret$r_in_in
  
  colnames(rho) <- c("Iteration", "outout", "outin", "inout", "inin")
  # colnames(rankRho) <- c("Iteration", "r-out-out", "r-out-in", "r-in-out", "r-in-in")
  edgelist[, 2] <- ret$targetNode + 1
  result <- list("assortcoeff" = rho, 
                 #  rankRho = rankRho,
                 "edgelist" = edgelist,
                 "iteration" = iteration,
                 "nattempts" = nattempts)
  if (rewire_history) {
    colnames(ret$history) <- c("Attempt", "Edge1", "Edge2", "Accepted")
    ret$history[, 1:3] <- ret$history[, 1:3] + 1
    result$history <- ret$history
  }
  return(result)
}

#' Degree preserving rewiring towards the target structure \code{e}.
#'
#' @param edgelist A two column matrix, each row represents a directed edge.
#' @param iteration An integer, number of rewiring iterations, each iteration
#'   consists of \code{nattempts} rewiring attempts.
#' @param nattempts An integer, number of rewiring attempts for each iteration.
#'   Default value equals the number of rows of edgelist.
#' @param eta An matrix, target structure \code{eta} generated by
#'   \code{wdnet::get_jointdist_undirected()}.
#' @param rewire_history Logical, whether the rewiring history should be returned.
#' @return Rewired edgelist, assortativity coefficient after each iteration, and
#'   rewiring history (including the index of sampled edges and rewiring
#'   result). For each rewiring attempt, two rows are sampled from the edgelist,
#'   for example Edge1:\{v_1, v_2\} and Edge2:\{v_3, v_4\}, we try to rewire the
#'   sampled edges as \{v_1, v_4\}, \{v_3, v_2\} (rewire type 1) or \{v_1,
#'   v_3\}, \{v_2, v_4\} (rewire type 2) with probability 1/2.
#'
#' @examples
#' set.seed(1234)
#' edgelist <- rpanet(1e4, directed = TRUE)$edgelist
#' ret1 <- wdnet:::get_jointdist_undirected(edgelist)
#' ret2 <- wdnet:::rewire_undirected(edgelist, eta = ret1$lbound$eta, iteration = 500)
#' plot(ret2$assortcoeff$Iteration, ret2$assortcoeff$Value)
#' ret3 <- wdnet:::get_jointdist_undirected(edgelist, target_assortcoeff = 0.5)
#' ret4 <- wdnet:::rewire_undirected(edgelist, eta = ret3$eta, iteration = 500)
#' plot(ret4$assortcoeff$Iteration, ret4$assortcoeff$Value)
#' 
rewire_undirected <- function(edgelist, eta, 
                              iteration = 1, nattempts = NULL, 
                              rewire_history = FALSE) {
  if (is.null(nattempts)) nattempts <- nrow(edgelist)
  
  edgelist <- as.matrix(edgelist)
  temp <- range(edgelist)
  stopifnot(temp[1] == 1)
  stopifnot(temp[2] == length(unique(c(edgelist))))
  degree <- data.frame(table(c(edgelist)))$Freq
  d_df <- data.frame(type = rownames(eta), index = seq_len(nrow(eta)) - 1)
  node1 <- edgelist[, 1]
  node2 <- edgelist[, 2]
  index1 <- d_df[match(degree[node1], d_df$type), "index"]
  index2 <- d_df[match(degree[node2], d_df$type), "index"]
  rm(d_df, temp)
  degree1 <- degree[c(node1, node2)]
  degree2 <- degree[c(node2, node1)]
  ret <- rewire_undirected_cpp(iteration, nattempts, 
                               node1, node2,
                               degree1, degree2,
                               index1, index2,
                               eta, rewire_history)
  rm(node1, node2, degree1, degree2, index1, index2)
  rho <- data.frame("Iteration" = c(0:iteration), "Value" = NA)
  rho[1, 2] <- assortcoeff(edgelist, directed = FALSE)
  rho[2:(iteration + 1), 2] <- ret$rho
  colnames(rho) <- c("Iteration", "Value")
  edgelist <- cbind(ret$node1, ret$node2)
  result <- list("assortcoeff" = rho,
                 "edgelist" = edgelist,
                 "iteration" = iteration,
                 "nattempts" = nattempts)
  if (rewire_history) {
    colnames(ret$history) <- c("Attempt", "Edge1", "Edge2", "RewireType", "Accepted")
    ret$history[, 1:4] <- ret$history[, 1:4] + 1
    result$history <- ret$history
  }
  return(result)
}
