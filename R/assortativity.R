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

#' @importFrom stats weighted.mean
#' @importFrom wdm wdm
NULL

## Directed assortativity coefficient

#' Compute the assortativity coefficient of a weighted and directed network.
#'
#' @param adj is an adjacency matrix of a weighted and directed network.
#' @param type which type of assortativity coefficient to compute: "outin"
#'   (default), "inin", "outout" or "inout"?
#'
#' @return a scalar of assortativity coefficient
#'
#' @references \itemize{ \item Foster, J.G., Foster, D.V., Grassberger, P. and
#' Paczuski, M. (2010). Edge direction and the structure of networks.
#' \emph{Proceedings of the National Academy of Sciences of the United States},
#' 107(24), 10815--10820. \item Yuan, Y. Zhang, P. and Yan, J. (2021).
#' Assortativity coefficients for weighted and directed networks. \emph{Journal
#' of Complex Networks}, 9(2), cnab017. }
#'
#' @note When the adjacency matrix is binary (i.e., directed but unweighted
#' networks), \code{dw_assort} returns the assortativity coefficient proposed in
#' Foster et al. (2010).
#'
#' @examples
#' ## Generate a network according to the Erd\"{o}s-Renyi model of order 20
#' ## and parameter p = 0.3
#' edge_ER <- rbinom(400,1,0.3)
#' weight_ER <- sapply(edge_ER, function(x) x*sample(3,1))
#' adj_ER <- matrix(weight_ER,20,20)
#' system.time(myassort <- wdnet:::dw_assort(adj_ER, type = "out-in"))
#' myassort
#' 

dw_assort <- function(adj, type = c("out-in", "in-in", "out-out", "in-out")) {
  stopifnot(dim(adj)[1] == dim(adj)[2])
  ## determine the location of edges in the network
  in_str <- colSums(adj)
  out_str <- rowSums(adj)
  vert_from <- unlist(apply(adj, 2, function(x){which(x > 0)}))
  number_to <- apply(adj, 2, function(x){length(which(x > 0) == TRUE)})
  temp_to <- cbind(seq(1:dim(adj)[1]),number_to)
  vert_to <- rep(temp_to[,1],temp_to[,2])
  weight <- adj[which(adj > 0)]
  type  <- match.arg(type)
  .type <- unlist(strsplit(type, "-"))
  x <- switch(.type[1], "out" = out_str, "in" = in_str)[vert_from]
  y <- switch(.type[2], "out" = out_str, "in" = in_str)[vert_to]
  weighted.cor <- function(x, y, w) {
    mean_x <- stats::weighted.mean(x, w)
    mean_y <- stats::weighted.mean(y, w)
    var_x <- sum((x - mean_x)^2 * w)
    var_y <- sum((y - mean_y)^2 * w)
    return(sum(w * (x - mean_x) * (y - mean_y)) / 
             sqrt(var_x * var_y))
  }
  return(weighted.cor(x, y, weight))
}

#' Compute the assortativity coefficient for a network.
#'
#' @param edgelist A two column matrix represents edges. If \code{NULL},
#'   \code{edgelist} and \code{edgeweight} will be extracted from the adjacency
#'   matrix \code{adj}.
#' @param edgeweight A vector represents the weight of edges. If \code{edgelist}
#'   is provided and \code{edgeweight} is \code{NULL}, all the edges will be
#'   considered have weight 1.
#' @param adj An adjacency matrix of a network.
#' @param directed Logical. Whether the edges will be considered as directed.
#' @param f1 A vector, represents the first feature of existing nodes. Number of
#'   nodes \code{= length(f1) = length(f2)}. Defined for directed networks. If
#'   \code{NULL}, out-strength will be used.
#' @param f2 A vector, represents the second feature of existing nodes. Defined
#'   for directed networks. If \code{NULL}, in-strength will be used.
#'
#' @return Assortativity coefficient for undirected networks, or four
#'   assortativity coefficients for directed networks.
#'
#' @references \itemize{ \item Foster, J.G., Foster, D.V., Grassberger, P. and
#'   Paczuski, M. (2010). Edge direction and the structure of networks.
#'   \emph{Proceedings of the National Academy of Sciences of the United
#'   States}, 107(24), 10815--10820. \item Yuan, Y. Zhang, P. and Yan, J.
#'   (2021). Assortativity coefficients for weighted and directed networks.
#'   \emph{Journal of Complex Networks}, 9(2), cnab017.}
#'
#' @note When the adjacency matrix is binary (i.e., directed but unweighted
#'   networks), \code{assort} returns the assortativity coefficient proposed in
#'   Foster et al. (2010).
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' control <- edgeweight.control(distribution = rgamma,
#'     dparams = list(shape = 5, scale = 0.2), shift = 0)
#' netwk <- rpanet(nstep = 10^4, control = control)
#' result <- assortcoeff(netwk$edgelist, edgeweight = netwk$edgeweight, directed = TRUE)
#' 
assortcoeff <- function(edgelist = NULL, edgeweight = NULL, adj = NULL, directed = TRUE, 
                        f1 = NULL, f2 = NULL) {
  if (is.null(edgelist)) {
    if (is.null(adj)) {
      stop('"edgelist" and "adj" can not both be NULL.')
    }
    temp <- adj_to_edges(adj = adj, directed = directed)
    edgelist <- temp$edgelist
    edgeweight <- temp$edgeweight
    rm(temp)
  }
  if (is.null(edgeweight)) {
    edgeweight <- rep(1, nrow(edgelist))
  }
  
  temp <- range(c(edgelist))
  stopifnot("Node index should start from 1." = temp[1] == 1)
  nnode <- temp[2]
  rm(temp)
  
  if ((! is.null(f1)) | (! is.null(f2))) {
    if (! directed) {
      stop("Node feature based assortativity coefficients are 
           defined for DIRECTED networks.")
    }
    return(dw_feature_assort(edgelist = edgelist, edgeweight = edgeweight, 
                             f1 = f1, f2 = f2))
  }
  
  if (! directed) {
    edgelist <- rbind(edgelist, edgelist[, c(2, 1)])
    edgeweight <- c(edgeweight, edgeweight)
  }
  sourceNode <- edgelist[, 1]
  targetNode <- edgelist[, 2]
  temp <- nodeStrength_cpp(snode = sourceNode, 
                           tnode = targetNode, 
                           nnode = nnode, 
                           weight = edgeweight,
                           weighted = TRUE)
  outs <- temp$outstrength
  ins <- temp$instrength
  rm(temp)
  sourceOut <- outs[sourceNode]
  targetIn <- ins[targetNode]
  if (! directed) {
    return(wdm::wdm(x = sourceOut, y = targetIn, 
                    weights = edgeweight, method = "pearson"))
  }
  sourceIn <- ins[sourceNode]
  targetOut <- outs[targetNode]
  return(list("outout" = wdm::wdm(x = sourceOut, y = targetOut,
                                   weights = edgeweight, method = "pearson"), 
              "outin" = wdm::wdm(x = sourceOut, y = targetIn, 
                                  weights = edgeweight, method = "pearson"), 
              "inout" = wdm::wdm(x = sourceIn, y = targetOut,
                                  weights = edgeweight, method = "pearson"),
              "inin" = wdm::wdm(x = sourceIn, y = targetIn,
                                 weights = edgeweight, method = "pearson")))
}

#' Node feature based assortativity coefficients of a weighted and directed
#' network.
#'
#' @param edgelist A two column matrix represents edges. If \code{NULL},
#'   \code{edgelist} and \code{edgeweight} will be extracted from the adjacency
#'   matrix \code{adj}.
#' @param edgeweight A vector represents the weight of edges. If \code{NULL},
#'   all the edges are considered have weight 1.
#' @param f1 A vector, represents the first feature of existing nodes. Number of
#'   nodes \code{= length(f1) = length(f2)}. Defined for directed networks. If
#'   \code{NULL}, out-strength will be used.
#' @param f2 A vector, represents the second feature of existing nodes. Defined
#'   for directed networks. If \code{NULL}, in-strength will be used.
#'
#' @return Directed weighted assortativity coefficients between source nodes'
#'   \code{f1} (or \code{f2}) and target nodes' \code{f2}(or \code{f1}).
#'
#' @examples
#' set.seed(123)
#' adj <- matrix(rbinom(400, 1, 0.2) * sample(1:3, 400, replace = TRUE), 20, 20)
#' f1 <- runif(20)
#' f2 <- abs(rnorm(20))
#' ret <- assortcoeff(adj = adj, f1 = f1, f2 = f2)
#' 
dw_feature_assort <- function(edgelist, edgeweight, f1, f2) {
  nnode <- max(edgelist)
  sourceNode <- edgelist[, 1]
  targetNode <- edgelist[, 2]
  if (is.null(f1) | is.null(f2)) {
    temp <- nodeStrength_cpp(snode = sourceNode, 
                             tnode = targetNode, 
                             weight = edgeweight, 
                             nnode = nnode, 
                             weighted = TRUE)
  }
  if (is.null(f1)) {
    f1 <- temp$outstrength
  }
  if (is.null(f2)) {
    f2 <- temp$instrength
  }
  stopifnot('Length of "f1" must equal number of nodes' = 
              length(f1) == nnode)
  stopifnot('Length of "f2" must equal number of nodes' = 
              length(f2) == nnode)
  sourceF1 <- f1[sourceNode]
  sourceF2 <- f2[sourceNode]
  targetF1 <- f1[targetNode]
  targetF2 <- f2[targetNode]
  ret <- list()
  ret$"f1-f1" <- wdm::wdm(x = sourceF1, y = targetF1,
                          weights = edgeweight, method = "pearson")
  ret$"f1-f2" <- wdm::wdm(x = sourceF1, y = targetF2,
                          weights = edgeweight, method = "pearson")
  ret$"f2-f1" <- wdm::wdm(x = sourceF2, y = targetF1,
                          weights = edgeweight, method = "pearson")
  ret$"f2-f2" <- wdm::wdm(x = sourceF2, y = targetF2,
                          weights = edgeweight, method = "pearson")
  return(ret)
}


#' Range of assortativity coefficient of a given unweighted network.
#'
#' Compute range of assortativity coefficients through convex optimization. The
#' problems are defined and solved via the \code{R} package \code{CVXR}. For
#' undirected networks, the function returns the range of the assortativity
#' coefficient. For directed networks, the function computes the range of an
#' assortativity coefficient while other assortativity coefficients are fixed.
#'
#' @param edgelist A two column matrix, each row represents an edge of the
#'   network.
#' @param directed Logical, whether the network is directed or not.
#' @param type The type of interested assortativity coefficient. For directed
#'   networks, it takes one of the values: "outout", "outin", "inout" and
#'   "inin". It will be ignored if the network is undirected.
#' @param constr A list represents the predetermined value or range imposed on
#'   other assortativity coefficients. It will be ignored if the network is 
#'   undirected.
#' @param control A list of parameters passed to \code{CVXR::solve()} for
#'   solving an appropriate \code{eta} when \code{target_assortcoeff} is
#'   provided. It will be ignored if \code{eta} is provided.
#'
#' @return Range of the interested assortativity coefficient; solved \code{eta}
#'   and its corresponding assortativity coefficients.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' edgelist <- rpanet(5e3, control = 
#'         scenario.control(alpha = 0.5, beta = 0.5))$edgelist
#' ret1 <- assortcoeff_range(edgelist, directed = TRUE, type = "outin", 
#'         constr = list("outout" = c(-0.3, 0.3), "inout" = 0.1))
#' ret2 <- rewire(edgelist, eta = ret1$lbound$eta, iteration = 100)
#' plot(ret2$assortcoeff$Iteration, ret2$assortcoeff$"outin")
#' ret3 <- rewire(edgelist, eta = ret1$ubound$eta, iteration = 500)
#' plot(ret3$assortcoeff$Iteration, ret3$assortcoeff$"outin")
#' }
#' 
assortcoeff_range <- function(edgelist, directed = TRUE, 
                              type = c("outout", "outin", "inout", "inin"),
                              control = solver.control(),
                              constr = list("outout" = NULL,
                                            "outin" = NULL,
                                            "inout" = NULL,
                                            "inin" = NULL)) {
  type <- match.arg(type)
  if (directed) {
    result <- get_jointdist_directed(edgelist = edgelist, 
                                     target_assortcoeff = constr,
                                     whichRange = type,
                                     control = control)
  }
  else {
    result <- get_jointdist_undirected(edgelist = edgelist,
                                       control = control)
  }
  result$lbound$e <- result$ubound$e <- NULL
  result
}