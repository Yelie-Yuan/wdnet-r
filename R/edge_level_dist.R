#' @importFrom CVXR Variable sum_entries Minimize Maximize Problem solve
NULL

#' Get the nodel-level joint distributions and some empirical distributions with
#' given edgelist.
#'
#' @param edgelist A two column matrix represents the directed edges of a
#'   network.
#' @param directed Logical, whether the network is directed.
#' @param joint_dist Logical, whether to return edge-level distributions.
#'
#' @return A list of distributions and degree vectors.
#'   
get_dist <- function(edgelist = NA, directed = TRUE, 
                     joint_dist = FALSE) {
  if (! directed) edgelist <- rbind(edgelist, edgelist[, c(2, 1)])
  edgelist <- as.matrix(edgelist)
  temp <- nodeStrength_cpp(startNode = edgelist[, 1], 
                           endNode = edgelist[, 2], 
                           nNodes = max(edgelist), 
                           weight = 1,
                           weighted = FALSE)
  outd <- temp$outstrength
  ind <- temp$instrength
  nedge <- nrow(edgelist)
  n_jk <- data.frame('outdegree' = outd, 'indegree' = ind)
  n_jk <- table(n_jk) / length(outd)
  d_out <- as.numeric(rownames(n_jk))
  d_in <- as.numeric(colnames(n_jk))
  p_out <- as.numeric(rowSums(n_jk))
  p_in <- as.numeric(colSums(n_jk))
  t1 <- n_jk * d_out; t1 <- t1 / sum(t1)
  t2 <- t(t(n_jk) * d_in); t2 <- t2 / sum(t2)
  # source-out
  q_s_out <- rowSums(t1)
  # target-in
  q_t_in <- colSums(t2)
  # source-in
  q_s_in <- colSums(t1)
  # target-out
  q_t_out <- rowSums(t2)
  e <- joint_e <- NA
  # other joint distributions
  if (joint_dist) {
    e <- list(
      'out-out' = table(data.frame(
        'source' = outd[edgelist[, 1]], 'target' = outd[edgelist[, 2]])) / nedge,
      'out-in' = table(data.frame(
        'source' = outd[edgelist[, 1]], 'target' = ind[edgelist[, 2]])) / nedge,
      'in-out' = table(data.frame(
        'source' = ind[edgelist[, 1]], 'target' = outd[edgelist[, 2]]))/ nedge,
      'in-in' = table(data.frame(
        'source' = ind[edgelist[, 1]], 'target' = ind[edgelist[, 2]])) / nedge)
    joint_e <- table(data.frame(
      'source' = paste(outd[edgelist[, 1]], ind[edgelist[, 1]], sep = '-'),
      'target' = paste(outd[edgelist[, 2]], ind[edgelist[, 2]], sep = '-')
    )) / nedge
  }
  list(n_jk  = n_jk, e = e, joint_e = joint_e,
       d_out = d_out, d_in = d_in,
       p_out = p_out, p_in = p_in,
       q_s_out = q_s_out, q_s_in = q_s_in,
       q_t_out = q_t_out, q_t_in = q_t_in)
}

#' Get the constraints for the optimization problem. This function is defined
#' for \code{directed_edge_level_dist}.
#'
#' @param constrs A list of constraints.
#' @param targetRho A list of target assortativity levels.
#' @param rho A list of variable objects.
#'
#' @return A list of constraints.
#'
get_constrs <- function(constrs, targetRho, rho) {
  for (type in names(targetRho)) {
    if (! is.null(targetRho[[type]])) {
      if (length(targetRho[[type]]) == 1) {
        constrs[[type]] <- rho[[type]] == targetRho[[type]]
      }
      else {
        constrs[[paste0(type, '_max')]] <- rho[[type]] <= max(targetRho[[type]])
        constrs[[paste0(type, '_min')]] <- rho[[type]] >= min(targetRho[[type]])
      }
    }
  }
  return(constrs)
}

#' Get the value of an object from the optimization problem. This function is
#' defined for \code{directed_edge_level_dist}.
#'
#' @param object An object from the optimization problem.
#' @param result A list returned from \code{CVXR::solve()}.
#' @param mydist A list of distributions returned from \code{wdnet:::get_dist()}.
#'
#' @return Value of the object.
#'   
get_values <- function(object, result, mydist) {
  out_out = result$getValue(object[['out-out']])
  out_in = result$getValue(object[['out-in']])
  in_out = result$getValue(object[['in-out']])
  in_in = result$getValue(object[['in-in']])
  if (deparse(substitute(object)) == 'e' & 
      ! any(is.na(out_out), is.na(out_in), 
            is.na(in_out), is.na(in_in))) {
    rownames(out_out) <- rownames(out_in) <- mydist$d_out
    colnames(in_out) <- colnames(out_out) <- mydist$d_out
    rownames(in_out) <- rownames(in_in) <- mydist$d_in
    colnames(out_in) <- colnames(in_in) <- mydist$d_in
  }
  list("out-out" = out_out, "out-in" = out_in,
       "in-out" = in_out, "in-in" = in_in)
}

#' Parameters passed to CVXR::solver().
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
#' control <- solver.control(solver = "OSQP", abstol = 1e-5)
solver.control <- function(solver = "ECOS", 
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

#' Edge-level distributions for directed networks with respect to desired
#' assortativity level(s).
#'
#' @param edgelist A two column matrix represents the directed edges of a
#'   network.
#' @param targetRho List, represents the specific values or ranges of
#'   assortativity coefficients.
#' @param f The convex function of the joint edge-level distribution to be
#'   minimized when \code{whichRange} is \code{NA}. Defaults to 0.
#' @param whichRange Character, "out-out", "out-in", "in-out" or "in-in". The
#'   range of interested assortativity level provide \code{targetRho} is
#'   satisfied.
#' @param control A list of parameters passed to \code{CVXR::solve()}.
#'   predetermined assortativity level(s) are satisfied.
#' @return Assortativity levels, edge-level distributions and the joint
#'   distribution of source nodes' out- and in-degrees and target nodes' out-
#'   and in-degrees.
#' @export
#'
#' @examples
#' edgelist <- rpanet(3000,
#'     control = panet.control(alpha = 0.3, beta = 0.1,
#'     gamma = 0.3, xi = 0.3, delta_out = 1, delta_in = 1))$edgelist
#' edge_assort(edgelist)
#' r <- list('out-out' = -0.1, 'out-in' = c(-0.2, 0.3), 'in-out' = 0.4)
#' ret1 <- directed_edge_level_dist(edgelist, targetRho = r, 
#'     whichRange = 'in-in')
#' r$'in-in' <- 0.4
#' r$'out-in' <- 0.3
#' ret2 <- directed_edge_level_dist(edgelist, targetRho = r, f = CVXR::norm2)
#' 
directed_edge_level_dist <- function(edgelist, 
                                     targetRho = list('out-out' = NULL, 'out-in' = NULL,
                                                      'in-out' = NULL, 'in-in' = NULL),
                                     f = function(x) 0,
                                     whichRange = NA,
                                     control = solver.control()) {
  mydist <- get_dist(edgelist = edgelist, directed = TRUE)
  m <- length(mydist$d_out)
  n <- length(mydist$d_in)
  
  s_outin <- c(t(mydist$n_jk * mydist$d_out))
  s_outin <- s_outin / sum(s_outin)
  t_outin <- c(t(mydist$n_jk) * mydist$d_in)
  t_outin <- t_outin / sum(t_outin)
  index_s <- s_outin != 0
  index_t <- t_outin != 0
  eMat <- CVXR::Variable(sum(index_s), sum(index_t), nonneg = TRUE)
  constrs <- list('rowSum' = CVXR::sum_entries(eMat, 1) == s_outin[index_s],
                  'colSum' = CVXR::sum_entries(eMat, 2) == t_outin[index_t])
  rm(s_outin, t_outin)
  
  mat1 <- matrix(0, m, m*n)
  for (i in 1:m) {
    temp <- c(1:n) + n * (i - 1)
    mat1[i, temp] <- 1
  }
  mat2 <- matrix(0, m*n, n)
  for (i in 1:n) {
    temp <- seq(0, m*n - 1, n) + i
    mat2[temp, i] <- 1  
  }
  e <- list("out-out" = mat1[, index_s] %*% eMat %*% t(mat1[, index_t]),
            "out-in"  = mat1[, index_s] %*% eMat %*% mat2[index_t, ],
            "in-out"  = t(mat2[index_s, ]) %*% eMat %*% t(mat1[, index_t]),
            "in-in"   = t(mat2[index_s, ]) %*% eMat %*% mat2[index_t, ])
  rm(mat1, mat2, m, n, i, temp)
  
  my_sigma <- function(j, q) {
    (sum(j^2 * q) - sum(j * q)^2)^0.5
  }
  sig <- list(s_out = my_sigma(mydist$d_out, mydist$q_s_out),
              s_in  = my_sigma(mydist$d_in, mydist$q_s_in),
              t_out = my_sigma(mydist$d_out, mydist$q_t_out),
              t_in  = my_sigma(mydist$d_in, mydist$q_t_in))
  
  rho <- list(
    "out-out" = t(mydist$d_out) %*% 
      (e$"out-out" - mydist$q_s_out %*% t(mydist$q_t_out)) %*% 
      mydist$d_out / sig$s_out / sig$t_out, 
    "out-in"  = t(mydist$d_out) %*% 
      (e$"out-in" - mydist$q_s_out %*% t(mydist$q_t_in)) %*% 
      mydist$d_in / sig$s_out / sig$t_in, 
    "in-out"  = t(mydist$d_in) %*% 
      (e$"in-out" - mydist$q_s_in %*% t(mydist$q_t_out)) %*% 
      mydist$d_out / sig$s_in / sig$t_out, 
    "in-in"   = t(mydist$d_in) %*% 
      (e$"in-in" - mydist$q_s_in %*% t(mydist$q_t_in)) %*% 
      mydist$d_in / sig$s_in / sig$t_in)
  
  # constrs$'out-out' <- rho$`out-out` <= 1
  # constrs$'out-in' <- rho$`out-in` <= 1
  # constrs$'in-out' <- rho$`in-out` <= 1
  # constrs$'in-in' <- rho$`in-in` <= 1
  name_eMat <- function(eMat, a = mydist$d_out, b = mydist$d_in, 
                        index_a = index_s, index_b = index_t) {
    temp <- paste0(rep(a, each = length(b)), '-',
                   rep(b, length(a)), split = '')
    colnames(eMat) <- temp[index_b]
    rownames(eMat) <- temp[index_a]
    names(attributes(eMat)$dimnames) <- c('source', 'target')
    eMat
  }
  constrs <- get_constrs(constrs, targetRho, rho)
  if (is.na(whichRange)) {
    problem <- CVXR::Problem(CVXR::Minimize(do.call(f, list(eMat))), constrs)
    result <- do.call(CVXR::solve, c(list(problem), control))
    if (result$status == 'solver_error') stop('SOLVER ERROR.')
    if (result$status == 'infeasible') stop('PROBLEM IS INFEASIBLE.')
    
    return(list(rho = get_values(rho, result, mydist),
                e = get_values(e, result, mydist),
                joint_e = name_eMat(result$getValue(eMat)))) 
  } else {
    whichRange <- switch (whichRange,
                          'out-out' = 1, 'out-in' = 2, 
                          'in-out' = 3, 'in-in' = 4)
    problem1 <- CVXR::Problem(CVXR::Minimize(rho[[whichRange]]), constrs)
    result1 <- do.call(CVXR::solve, c(list(problem1), control))
    if (result1$status == 'solver_error') stop('SOLVER ERROR.')
    if (result1$status == 'infeasible') stop('PROBLEM IS INFEASIBLE.')
    
    problem2 <- CVXR::Problem(CVXR::Maximize(rho[[whichRange]]), constrs)
    result2 <- do.call(CVXR::solve, c(list(problem2), control))

    if (result2$status == 'solver_error') stop('SOLVER ERROR.')
    if (result2$status == 'infeasible') stop('PROBLEM IS INFEASIBLE.')
    
    return(list(range = c(result1$getValue(rho[[whichRange]]), 
                          result2$getValue(rho[[whichRange]])),
                lbound = list(rho = get_values(rho, result1, mydist),
                              e = get_values(e, result1, mydist),
                              joint_e = name_eMat(result1$getValue(eMat))),
                ubound = list(rho = get_values(rho, result2, mydist),
                              e = get_values(e, result2, mydist),
                              joint_e = name_eMat(result2$getValue(eMat)))))
  }
}

#' Edge-level distribution for undirected networks with respect to desired
#' assortativity level.
#'
#' @param edgelist A two column matrix represents the undirected edges of a
#'   network.
#' @param targetRho Numeric, represents the predetermined assortativity
#'   coefficient. If \code{NA}, the range of assortativity coefficient and
#'   corresponding edge-level distribution are returned.
#' @param f The convex function of the edge-level distribution to be minimized
#'   when \code{targetRho} is not \code{NA}. Default is 0.
#' @param control A list of parameters passed to \code{CVXR::solve()}.
#'
#' @return Assortativity level and corresponding edge-level distribution.
#' @export
#'
#' @examples
#' set.seed(1234)
#' edgelist <- matrix(sample(1:10, 500, replace = TRUE), ncol = 2)
#' ret1 <- undirected_edge_level_dist(edgelist, f = CVXR::norm2)
#' ret2 <- undirected_edge_level_dist(edgelist, targetRho = 0.6, f = CVXR::norm2)
#' 
undirected_edge_level_dist <- function(edgelist, targetRho = NA, 
                                       f = function(x) 0,
                                       control = solver.control()) {
  stopifnot((targetRho <= 1 & targetRho >= -1) | is.na(targetRho))
  mydist <- get_dist(edgelist = edgelist, directed = FALSE)
  k <- mydist$d_out
  q_k <- mydist$q_s_out
  rm(mydist)
  name_eMat <- function(eMat, k) {
    colnames(eMat) <- rownames(eMat) <- k
    eMat
  }
  if ((! is.na(targetRho)) & targetRho == 0) {
    return(list(rho = 0, 
                e = name_eMat(q_k %*% t(q_k), k)))
  }
  n <- length(k)
  sig2 <- sum(k^2 * q_k) - (sum(k * q_k))^2
  eMat <- CVXR::Variable(n, n, nonneg = TRUE)
  rho <- t(k) %*% (eMat - q_k %*% t(q_k)) %*% k / sig2
  constrs <- list(CVXR::sum_entries(eMat, 1) == q_k, 
                  eMat == t(eMat))
  
  if (! is.na(targetRho)) {
    constrs$'rho' <- rho == targetRho
    problem <- CVXR::Problem(CVXR::Minimize(do.call(f, list(eMat))), constrs)
    result <- do.call(CVXR::solve, c(list(problem), control))
    if (result$status == 'solver_error') stop('SOLVER ERROR.')
    if (result$status == 'infeasible') stop('PROBLEM IS INFEASIBLE.')
    return(list(rho = result$getValue(rho),
                e = name_eMat(result$getValue(eMat), k)))
  } else {
    # constrs$'rho' <- rho <= 1
    problem1 <- CVXR::Problem(CVXR::Minimize(rho), constrs)
    result1 <- do.call(CVXR::solve, c(list(problem1), control))
    if (result1$status == 'solver_error') stop('SOLVER ERROR.')
    if (result1$status == 'infeasible') stop('PROBLEM IS INFEASIBLE.')
    
    problem2 <- CVXR::Problem(CVXR::Maximize(rho), constrs)
    result2 <- do.call(CVXR::solve, c(list(problem2), control))
    if (result2$status == 'solver_error') stop('SOLVER ERROR.')
    if (result2$status == 'infeasible') stop('PROBLEM IS INFEASIBLE.')
    
    return(list(range = c(result1$getValue(rho), result2$getValue(rho)),
                lbound = list(rho = result1$getValue(rho),
                              e = name_eMat(result1$getValue(eMat), k)),
                ubound = list(rho = result2$getValue(rho),
                              e = name_eMat(result2$getValue(eMat), k))))
  }
}
