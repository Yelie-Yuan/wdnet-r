#' @importFrom stats runif rpois
NULL

#' Generate a growing preferential attachment network.
#'
#' @param edgelist A two column matrix represents the seed graph.
#' @param edgeweight A vector represents the weight of edges of the seed graph.
#'   Its length equals the number of edges of the seed graph. If NA, all the
#'   edges of the seed graph have weight 1.
#' @param nsteps Number of steps when generating a network.
#' @param control A list of parameters to be used when generate network.
#'
#' @return A list with the following components: edgelist, edgeweight, out- and
#'   in-strength, number of edges per step (m), scenario of each new edge
#'   (1~alpha, 2~beta, 3~gamma, 4~xi, 5~rho). The edges in the seed graph 
#'   are denoted as scenario 0.
#' @export
#'
#' @examples
#' net <- rpanet_simple(nsteps = 100, 
#'         control = panet.control(alpha = 0.4, beta = 0, gamma = 0.6))
#' net <- rpanet_simple(edgelist = matrix(c(1:8), ncol = 2), nsteps = 10^3,
#'       control = panet.control(mdist = stats::rpois,
#'       mpar = list(lambda = 1), m_c = 1,
#'       wdist = stats::runif, wpar = list(min = 1, max = 10), w_c = 0))


rpanet_simple <- function(nsteps = 10^3, edgelist = matrix(c(1, 2), ncol = 2), 
                          edgeweight = NA, control = panet.control()) {
  stopifnot("nsteps must be greater than 0." = nsteps > 0)
  stopifnot("alpha + beta + bamma + xi + rho must be less or equal to 1." = 
              control$alpha + control$beta + control$gamma + 
              control$xi + control$rho <= 1)
  temp <- c(edgelist)
  nnode <- max(temp)
  stopifnot("Nodes' index should be consecutive natural numbers start from 1." =
              sum(! duplicated(temp)) == nnode)
  if (is.na(edgeweight[1])) edgeweight[1:nrow(edgelist)] <- 1
  stopifnot(length(edgeweight) == nrow(edgelist))
  
  if (! is.numeric(control$mdist)) {
    m <- do.call(control$mdist, c(nsteps, control$mpar)) + control$m_c
  }
  else m <- rep(control$mdist + control$m_c, nsteps)
  stopifnot("Number of new edges per step must be positive integers." = m %% 1 == 0)
  stopifnot("Number of new edges per step must be positive integers." = m > 0)
  sum_m <- sum(m)
  if (! is.numeric(control$wdist)) {
    w <- do.call(control$wdist, c(sum_m, control$wpar)) + control$w_c
  }
  else w <- rep(control$wdist + control$w_c, sum_m)
  stopifnot("Edge weight must be greater than 0." = w > 0)
  
  strength <- nodeStrength_cpp(edgelist[, 1], edgelist[, 2], 
                               edgeweight, nnode, weighted = TRUE)
  outstrength <- instrength <- rep(0, 2 * sum_m + nnode)
  outstrength[1:nnode] <- strength$outstrength
  instrength[1:nnode] <- strength$instrength
  sumstrength <- sum(strength$outstrength)
  
  control_cpp <- c(control$alpha, control$beta, control$gamma, control$xi, 
                   control$delta_out, control$delta_in)
  
  ret <- rpanet_simple_cpp(nsteps, control_cpp, m, w,
                           outstrength, instrength, sumstrength,
                           nnode)
  control$delta <- NULL
  list("edgelist" = rbind(edgelist, cbind(ret$startnode, ret$endnode)), 
       "edgeweight" = c(edgeweight, w),
       "out-strength" = ret$outstrength, 
       "in-strength" = ret$instrength, 
       "edgescenario" = c(rep(0, nrow(edgelist)), ret$edgescenario),
       "m" = m, 
       "control" = control)
}