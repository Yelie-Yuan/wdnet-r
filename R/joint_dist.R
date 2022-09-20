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

#' @importFrom CVXR Variable sum_entries Minimize Maximize Problem solve
NULL

#' Get the node-level joint distributions and some empirical distributions with
#' given edgelist.
#'
#' @param edgelist A two column matrix represents the directed edges of a
#'   network.
#' @param directed Logical, whether the network is directed.
#' @param joint_dist Logical, whether to return edge-level distributions.
#'
#' @return A list of distributions and degree vectors.
#' 
#' @keywords internal
#'   
get_dist <- function(edgelist = NA, directed = TRUE, 
                     joint_dist = FALSE) {
  if (! directed) edgelist <- rbind(edgelist, edgelist[, c(2, 1)])
  edgelist <- as.matrix(edgelist)
  temp <- node_strength_cpp(snode = edgelist[, 1], 
                            tnode = edgelist[, 2], 
                            nnode = max(edgelist), 
                            weight = 1,
                            weighted = FALSE)
  
  outd <- temp$outstrength
  ind <- temp$instrength
  nedge <- nrow(edgelist)
  nu <- data.frame("outdegree" = outd, "indegree" = ind)
  nu <- table(nu) / length(outd)
  d_out <- as.numeric(rownames(nu))
  d_in <- as.numeric(colnames(nu))
  
  # get_rank <- function(d, degree_seq) {
  #   degree_seq <- sort(degree_seq)
  #   temp <- data.frame(degree_seq = unique(degree_seq),
  #                      degree_rank = unique(rank(degree_seq, ties.method = "average")))
  #   unlist(sapply(d, function(d1) {
  #     n <- which(temp$degree_seq %in% d1)
  #     ifelse(length(n) == 0, 0, temp$degree_rank[n])
  #   }))
  # }
  # r_s_out <- get_rank(d_out, outd[edgelist[, 1]])
  # r_t_out <- get_rank(d_out, outd[edgelist[, 2]])
  # r_s_in <- get_rank(d_in, ind[edgelist[, 1]])
  # r_t_in <- get_rank(d_in, ind[edgelist[, 2]])
  # denom <- max(r_s_out, r_t_out, r_s_in, r_t_in) / 1e4
  # r_s_out <- r_s_out / denom
  # r_t_out <- r_t_out / denom
  # r_s_in <- r_s_in / denom
  # r_t_in <- r_t_in / denom
  
  p_out <- as.numeric(rowSums(nu))
  p_in <- as.numeric(colSums(nu))
  t1 <- nu * d_out; t1 <- t1 / sum(t1)
  t2 <- t(t(nu) * d_in); t2 <- t2 / sum(t2)
  # source-out
  q_s_out <- rowSums(t1)
  # target-in
  q_t_in <- colSums(t2)
  # source-in
  q_s_in <- colSums(t1)
  # target-out
  q_t_out <- rowSums(t2)
  e <- eta <- NA
  # other joint distributions
  if (joint_dist) {
    e <- list(
      "outout" = table(data.frame(
        "source" = outd[edgelist[, 1]], "target" = outd[edgelist[, 2]])) / nedge,
      "outin" = table(data.frame(
        "source" = outd[edgelist[, 1]], "target" = ind[edgelist[, 2]])) / nedge,
      "inout" = table(data.frame(
        "source" = ind[edgelist[, 1]], "target" = outd[edgelist[, 2]]))/ nedge,
      "inin" = table(data.frame(
        "source" = ind[edgelist[, 1]], "target" = ind[edgelist[, 2]])) / nedge)
    eta <- table(data.frame(
      "source" = paste(outd[edgelist[, 1]], ind[edgelist[, 1]], sep = "-"),
      "target" = paste(outd[edgelist[, 2]], ind[edgelist[, 2]], sep = "-")
    )) / nedge
  }
  list(nu  = nu, e = e, eta = eta,
       d_out = d_out, d_in = d_in,
       p_out = p_out, p_in = p_in, 
       #  r_s_out = r_s_out, r_s_in = r_s_in, 
       #  r_t_in = r_t_in, r_t_out = r_t_out,
       q_s_out = q_s_out, q_s_in = q_s_in,
       q_t_out = q_t_out, q_t_in = q_t_in)
}

#' Get the constraints for the optimization problem. This function is defined
#' for \code{get_eta_directed()}.
#'
#' @param constrs A list of constraints.
#' @param target.assortcoef A list of target assortativity levels.
#' @param rho A list of variable objects.
#'
#' @return A list of constraints.
#' 
#' @keywords internal
#' 
get_constr <- function(constrs, target.assortcoef, rho) {
  for (type in names(target.assortcoef)) {
    if (! is.null(target.assortcoef[[type]])) {
      if (length(target.assortcoef[[type]]) == 1) {
        constrs[[type]] <- rho[[type]] == target.assortcoef[[type]]
      }
      else {
        constrs[[paste0(type, "_max")]] <- rho[[type]] <= max(target.assortcoef[[type]])
        constrs[[paste0(type, "_min")]] <- rho[[type]] >= min(target.assortcoef[[type]])
      }
    }
  }
  return(constrs)
}

#' Get the value of an object from the optimization problem. This function is
#' defined for \code{get_eta_directed()}.
#'
#' @param object An object from the optimization problem.
#' @param result A list returned from \code{CVXR::solve()}.
#' @param mydist A list returned from \code{get_dist()}.
#'
#' @return Value of the object.
#' 
#' @keywords internal
#'
get_values <- function(object, result, mydist) {
  # if ("r-out-out" %in% names(object)) {
  #   return(list(
  #     "r-out-out" = result$getValue(object[["r-out-out"]]),
  #     "r-out-in" = result$getValue(object[["r-out-in"]]),
  #     "r-in-out" = result$getValue(object[["r-in-out"]]),
  #     "r-in-in" = result$getValue(object[["r-in-in"]])
  #   ))
  # }
  out_out <- result$getValue(object[["outout"]])
  out_in <- result$getValue(object[["outin"]])
  in_out <- result$getValue(object[["inout"]])
  in_in <- result$getValue(object[["inin"]])
  if (deparse(substitute(object)) == "e" & 
      ! any(is.na(out_out), is.na(out_in), 
            is.na(in_out), is.na(in_in))) {
    rownames(out_out) <- rownames(out_in) <- mydist$d_out
    colnames(in_out) <- colnames(out_out) <- mydist$d_out
    rownames(in_out) <- rownames(in_in) <- mydist$d_in
    colnames(out_in) <- colnames(in_in) <- mydist$d_in
  }
  list("outout" = out_out, "outin" = out_in,
       "inout" = in_out, "inin" = in_in)
}

#' Parameters passed to CVXR::solver().
#'
#' Defined for the convex optimization problems for solving \code{eta}. The
#' control list is passed to \code{dprewire} and \code{dprewire.range}.
#'
#' @param solver (Optional) A string indicating the solver to use. Defaults to
#'   "ECOS".
#' @param ignore_dcp (Optional) A logical value indicating whether to override
#'   the DCP check for a problem.
#' @param warm_start (Optional) A logical value indicating whether the previous
#'   solver result should be used to warm start.
#' @param verbose (Optional) A logical value indicating whether to print
#'   additional solver output.
#' @param parallel (Optional) A logical value indicating whether to solve in
#'   parallel if the problem is separable.
#' @param gp (Optional) A logical value indicating whether the problem is a
#'   geometric program. Defaults to FALSE.
#' @param feastol The feasible tolerance on the primal and dual residual.
#' @param reltol The relative tolerance on the duality gap.
#' @param abstol The absolute tolerance on the duality gap.
#' @param num_iter The maximum number of iterations.
#' @param ... Additional options that will be passed to the specific solver. In
#'   general, these options will override any default settings imposed by CVXR.
#'
#' @return A list containing the parameters.
#' @export
#'
#' @examples
#' control <- cvxr.control(solver = "OSQP", abstol = 1e-5)
cvxr.control <- function(solver = "ECOS", 
                         ignore_dcp = FALSE,
                         warm_start = FALSE,
                         verbose = FALSE,
                         parallel = FALSE,
                         gp = FALSE,
                         feastol = NULL,
                         reltol = NULL,
                         abstol = NULL,
                         num_iter = NULL,
                         ...) {
  return(list(solver = solver,
              ignore_dcp = ignore_dcp,
              warm_start = warm_start,
              verbose = verbose,
              parallel = parallel,
              gp = gp,
              feastol = feastol,
              reltol = reltol,
              abstol = abstol,
              num_iter = num_iter, 
              ...))
}

#' Compute edge-level distributions for directed networks with respect to
#' desired assortativity level(s).
#'
#' @param edgelist A two column matrix represents the directed edges of a
#'   network.
#' @param target.assortcoef List, represents the predetermined value or range
#'   of assortativity coefficients.
#' @param eta.obj A convex function of \code{eta} to be minimized when
#'   \code{which.range} is \code{NULL}. Defaults to 0.
#' @param which.range Character, "outout", "outin", "inout" or "inin".
#'   Represents the interested degree based assortativity coefficient. Default
#'   is \code{NA}.
#' @param control A list of parameters passed to \code{CVXR::solve()} when
#'   solving for \code{eta} or computing the range of assortativity coefficient.
#' @return Assortativity coefficients and joint distributions. If
#'   \code{which.range} is specified, the range of the interested coefficient
#'   and the corresponding joint distributions will be returned, provided the
#'   predetermined \code{target.assortcoef} is satisfied.
#' 
#' @keywords internal
#' 
get_eta_directed <- function(edgelist, 
                             target.assortcoef = list("outout" = NULL, "outin" = NULL,
                                                      "inout" = NULL, "inin" = NULL),
                             # target_rank_assortcoef = list("r-out-out" = NULL, "r-out-in" = NULL,
                             #                      "r-in-out" = NULL, "r-in-in" = NULL),
                             eta.obj = function(x) 0, which.range = NULL, 
                             control = cvxr.control()) {
  # stopifnot(all(names(target_rank_assortcoef) %in% c("r-out-out", "r-out-in", 
  #                                           "r-in-out", "r-in-in")))
  stopifnot(all(names(target.assortcoef) %in% c("outout", "outin", 
                                                "inout", "inin")))
  mydist <- get_dist(edgelist = edgelist, directed = TRUE)
  m <- length(mydist$d_out)
  n <- length(mydist$d_in)
  
  s_outin <- c(t(mydist$nu * mydist$d_out))
  s_outin <- s_outin / sum(s_outin)
  t_outin <- c(t(mydist$nu) * mydist$d_in)
  t_outin <- t_outin / sum(t_outin)
  index_s <- s_outin != 0
  index_t <- t_outin != 0
  eMat <- CVXR::Variable(sum(index_s), sum(index_t), nonneg = TRUE)
  constrs <- list("rowSum" = CVXR::sum_entries(eMat, 1) == s_outin[index_s],
                  "colSum" = CVXR::sum_entries(eMat, 2) == t_outin[index_t])
  rm(s_outin, t_outin)
  
  mat1 <- kronecker(diag(rep(1, m)), t(rep(1, n)))
  mat2 <- kronecker(rep(1, m), diag(rep(1, n)))
  e <- list("outout" = mat1[, index_s] %*% eMat %*% t(mat1[, index_t]),
            "outin"  = mat1[, index_s] %*% eMat %*% mat2[index_t, ],
            "inout"  = t(mat2[index_s, ]) %*% eMat %*% t(mat1[, index_t]),
            "inin"   = t(mat2[index_s, ]) %*% eMat %*% mat2[index_t, ])
  rm(mat1, mat2, m, n)
  
  my_sigma <- function(j, q) {
    (sum(j^2 * q) - sum(j * q)^2)^0.5
  }
  sig <- list(s_out = my_sigma(mydist$d_out, mydist$q_s_out),
              s_in  = my_sigma(mydist$d_in, mydist$q_s_in),
              t_out = my_sigma(mydist$d_out, mydist$q_t_out),
              t_in  = my_sigma(mydist$d_in, mydist$q_t_in))
  # rankSig <- list(s_out = my_sigma(mydist$r_s_out, mydist$q_s_out),
  #                 s_in  = my_sigma(mydist$r_s_in, mydist$q_s_in),
  #                 t_out = my_sigma(mydist$r_t_out, mydist$q_t_out),
  #                 t_in  = my_sigma(mydist$r_t_in, mydist$q_t_in))
  
  rho <- list(
    "outout" = t(mydist$d_out) %*% 
      (e$"outout" - mydist$q_s_out %*% t(mydist$q_t_out)) %*% 
      mydist$d_out / sig$s_out / sig$t_out, 
    "outin"  = t(mydist$d_out) %*% 
      (e$"outin" - mydist$q_s_out %*% t(mydist$q_t_in)) %*% 
      mydist$d_in / sig$s_out / sig$t_in, 
    "inout"  = t(mydist$d_in) %*% 
      (e$"inout" - mydist$q_s_in %*% t(mydist$q_t_out)) %*% 
      mydist$d_out / sig$s_in / sig$t_out, 
    "inin"   = t(mydist$d_in) %*% 
      (e$"inin" - mydist$q_s_in %*% t(mydist$q_t_in)) %*% 
      mydist$d_in / sig$s_in / sig$t_in)
  
  # rankRho <- list(
  #   "r-out-out" = t(mydist$r_s_out) %*% 
  #     (e$"outout" - mydist$q_s_out %*% t(mydist$q_t_out)) %*% 
  #     mydist$r_t_out / rankSig$s_out / rankSig$t_out, 
  #   "r-out-in"  = t(mydist$r_s_out) %*% 
  #     (e$"outin" - mydist$q_s_out %*% t(mydist$q_t_in)) %*% 
  #     mydist$r_t_in / rankSig$s_out / rankSig$t_in, 
  #   "r-in-out"  = t(mydist$r_s_in) %*% 
  #     (e$"inout" - mydist$q_s_in %*% t(mydist$q_t_out)) %*% 
  #     mydist$r_t_out / rankSig$s_in / rankSig$t_out, 
  #   "r-in-in"   = t(mydist$r_s_in) %*% 
  #     (e$"inin" - mydist$q_s_in %*% t(mydist$q_t_in)) %*% 
  #     mydist$r_t_in / rankSig$s_in / rankSig$t_in)
  
  name_eMat <- function(eMat, a = mydist$d_out, b = mydist$d_in, 
                        index_a = index_s, index_b = index_t) {
    temp <- paste0(rep(a, each = length(b)), "-",
                   rep(b, length(a)), split = "")
    colnames(eMat) <- temp[index_b]
    rownames(eMat) <- temp[index_a]
    names(attributes(eMat)$dimnames) <- c("source", "target")
    eMat
  }
  constrs <- get_constr(constrs, target.assortcoef, rho)
  # constrs <- get_constr(constrs, target_rank_assortcoef, rankRho)
  retitems <- c("value", "status", "solver", "solve_time", "setup_time", "num_iters")
  if (is.null(which.range)) {
    problem <- CVXR::Problem(CVXR::Minimize(do.call(eta.obj, list(eMat))), constrs)
    result <- do.call(CVXR::solve, c(list(problem), control))
    ret <- result[retitems]
    if (result$status == "solver_error" | result$status == "infeasible") {
      warning(paste0("Solver status: ", result$status))
      return(ret)
    }
    ret$assortcoef <- get_values(rho, result, mydist)
    # ret$e <- get_values(e, result, mydist)
    ret$eta <- name_eMat(result$getValue(eMat))
    return(ret)
  } else {
    # tempRho <- append(rho, rankRho)
    tempRho <- rho
    stopifnot("'which.range' is not valid." = which.range %in% names(tempRho))
    problem1 <- CVXR::Problem(CVXR::Minimize(tempRho[[which.range]]), constrs)
    result1 <- do.call(CVXR::solve, c(list(problem1), control))
    if (result1$status == "solver_error" | result1$status == "infeasible") {
      warning(paste0("Lower bound solver status: ", result1$status))
    }
    
    problem2 <- CVXR::Problem(CVXR::Maximize(tempRho[[which.range]]), constrs)
    result2 <- do.call(CVXR::solve, c(list(problem2), control))
    if (result2$status == "solver_error" | result2$status == "infeasible") {
      warning(paste0("Upper bound solver status: ", result2$status))
    }
    return(list("range" = c(result1$getValue(tempRho[[which.range]]), 
                            result2$getValue(tempRho[[which.range]])),
                "lbound.solver.result" = result1[retitems], 
                "ubound.solver.result" = result2[retitems]))
  }
}

#' Compute edge-level distribution for undirected networks with respect to
#' desired assortativity level.
#'
#' @param edgelist A two column matrix represents the undirected edges of a
#'   network.
#' @param target.assortcoef Numeric, represents the predetermined assortativity
#'   coefficient. If \code{NA}, the range of assortativity coefficient and
#'   corresponding joint distribution are returned.
#' @param eta.obj A convex function of \code{eta} to be minimized when
#'   \code{target.assortcoef} is not \code{NA}. Defaults to 0.
#' @param control A list of parameters passed to \code{CVXR::solve()} when
#'   solving for \code{eta} or computing the range of assortativity coefficient.
#'
#' @return Assortativity level and corresponding edge-level distribution.
#' 
#' @keywords internal
#'
get_eta_undirected <- function(edgelist, target.assortcoef = NULL, 
                               eta.obj = function(x) 0,
                               control = cvxr.control()) {
  stopifnot((target.assortcoef <= 1 & target.assortcoef >= -1) | is.null(target.assortcoef))
  mydist <- get_dist(edgelist = edgelist, directed = FALSE)
  k <- mydist$d_out
  q_k <- mydist$q_s_out
  rm(mydist)
  name_eMat <- function(eMat, k) {
    colnames(eMat) <- rownames(eMat) <- k
    eMat
  }
  if (! is.null(target.assortcoef)) {
    if (target.assortcoef == 0) {
      return(list("assortcoef" = 0, 
                  "eta" = name_eMat(q_k %*% t(q_k), k)))
    }
  }
  n <- length(k)
  sig2 <- sum(k^2 * q_k) - (sum(k * q_k))^2
  eMat <- CVXR::Variable(n, n, nonneg = TRUE)
  rho <- t(k) %*% (eMat - q_k %*% t(q_k)) %*% k / sig2
  constrs <- list(CVXR::sum_entries(eMat, 1) == q_k, 
                  eMat == t(eMat))
  retitems <- c("value", "status", "solver", "solve_time", "setup_time", "num_iters")
  if (! is.null(target.assortcoef)) {
    constrs$"rho" <- rho == target.assortcoef
    problem <- CVXR::Problem(CVXR::Minimize(do.call(eta.obj, list(eMat))), constrs)
    result <- do.call(CVXR::solve, c(list(problem), control))
    ret <- result[retitems]
    if (result$status == "solver_error" | result$status == "infeasible") {
      warning(paste0("Solver status: ", result$status))
      return(ret)
    }
    ret$assortcoef <- result$getValue(rho)
    ret$eta <- name_eMat(result$getValue(eMat), k)
    return(ret)
  } else {
    # constrs$"rho" <- rho <= 1
    problem1 <- CVXR::Problem(CVXR::Minimize(rho), constrs)
    result1 <- do.call(CVXR::solve, c(list(problem1), control))
    if (result1$status == "solver_error" | result1$status == "infeasible") {
      warning(paste0("Lower bound solver status: ", result1$status))
    }
    problem2 <- CVXR::Problem(CVXR::Maximize(rho), constrs)
    result2 <- do.call(CVXR::solve, c(list(problem2), control))
    if (result2$status == "solver_error" | result2$status == "infeasible") {
      warning(paste0("Upper bound solver status: ", result2$status))
    }
    return(list("range" = c(result1$getValue(rho), result2$getValue(rho)),
                "lbound.solver.result" = result1[retitems], 
                "ubound.solver.result" = result2[retitems]))
  }
}
