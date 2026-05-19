##
## wdnet: Weighted directed network
## Copyright (C) 2026  Yelie Yuan, Tiandong Wang, Jun Yan and Panpan Zhang
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
#'   Paczuski, M. (2010). Edge direction and the structure of networks.
#'   \emph{Proceedings of the National Academy of Sciences of the United
#'   States}, 107(24), 10815--10820. \item Yuan, Y. Zhang, P. and Yan, J.
#'   (2021).
#' Assortativity coefficients for weighted and directed networks. \emph{Journal
#' of Complex Networks}, 9(2), cnab017. }
#'
#' @note When the adjacency matrix is binary (i.e., directed but unweighted
#'   networks), \code{dw_assort} returns the assortativity coefficient proposed
#'   in Foster et al. (2010).
#'
#' @keywords internal
#'

dw_assort <- function(adj, type = c("outin", "inin", "outout", "inout")) {
  stopifnot(dim(adj)[1] == dim(adj)[2])
  ## determine the location of edges in the network
  in_str <- colSums(adj)
  out_str <- rowSums(adj)
  vert_from <- unlist(apply(adj, 2, function(x) {
    which(x > 0)
  }))
  number_to <- apply(adj, 2, function(x) {
    length(which(x > 0) == TRUE)
  })
  temp_to <- cbind(seq(1:dim(adj)[1]), number_to)
  vert_to <- rep(temp_to[, 1], temp_to[, 2])
  weight <- adj[which(adj > 0)]
  type <- match.arg(type)
  .type <- unlist(strsplit(type, "-"))
  x <- switch(.type[1],
    "out" = out_str,
    "in" = in_str
  )[vert_from]
  y <- switch(.type[2],
    "out" = out_str,
    "in" = in_str
  )[vert_to]
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


#' Compute the assortativity coefficient(s) for a network.
#'
#' @param netwk A \code{wdnet} object that represents the network. If
#'   \code{NULL}, the function will compute the coefficient using either
#'   \code{edgelist} and \code{edgeweight}, or \code{adj}.
#' @param edgelist A two-column matrix representing edges.
#' @param edgeweight A numeric vector of edge weights with the same length as
#'   the number of rows in edgelist. If \code{NULL}, all edges will be assigned
#'   weight 1.
#' @param adj The adjacency matrix of a network.
#' @param directed Logical. Indicates whether the edges in \code{edgelist} or
#'   \code{adj} are directed. It will be omitted if \code{netwk} is provided.
#' @param f1 A vector representing the first feature of existing nodes. The
#'   number of nodes should be equal to the length of both \code{f1} and
#'   \code{f2}. Defined for directed networks. If \code{NULL}, out-strength will
#'   be used.
#' @param f2 A vector representing the second feature of existing nodes. Defined
#'   for directed networks. If \code{NULL}, in-strength will be used.
#' @param weighted.rank Logical. If supplied, assortativity is computed using
#'   ranked node strengths; \code{TRUE} uses weighted midranks,
#'   \code{FALSE} uses average ranks.
#'
#' @return Assortativity coefficient for undirected networks, or a list of four
#'   assortativity coefficients for directed networks.
#'
#' @references \itemize{ \item Foster, J.G., Foster, D.V., Grassberger, P. and
#'   Paczuski, M. (2010). Edge direction and the structure of networks.
#'   \emph{Proceedings of the National Academy of Sciences of the United
#'   States}, 107(24), 10815--10820. \item Yuan, Y. Zhang, P. and Yan, J.
#'   (2021). Assortativity coefficients for weighted and directed networks.
#'   \emph{Journal of Complex Networks}, 9(2), cnab017.
#'   \item Shen, A., Feng, Q., Yan, J. and Zhang, P. (2025).
#'   Rank-based assortativity for weighted, directed networks.
#'   \emph{Journal of Complex Networks}, 13(2), cnaf002.}
#'
#' @note When the adjacency matrix is binary (i.e., directed but unweighted
#'   networks), \code{assortcoef} returns the assortativity coefficient proposed
#'   in Foster et al. (2010).
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' control <- rpa_control_edgeweight(
#'   sampler = function(n) rgamma(n, shape = 5, scale = 0.2)
#' )
#' netwk <- rpanet(nstep = 10^4, control = control)
#' assortcoef(netwk)
#' assortcoef(
#'   edgelist = netwk$edgelist,
#'   edgeweight = netwk$edge.attr$weight,
#'   directed = TRUE
#' )
#'
assortcoef <- function(
    netwk,
    edgelist,
    edgeweight,
    adj,
    directed,
    f1,
    f2,
    weighted.rank) {
  netwk <- create_wdnet(
    netwk = netwk,
    edgelist = edgelist,
    edgeweight = edgeweight,
    directed = directed,
    adj = adj,
    weighted = TRUE
  )

  edgelist <- netwk$edgelist
  edgeweight <- netwk$edge.attr$weight
  directed <- netwk$directed
  nnode <- max(edgelist)

  if ((!missing(f1)) || (!missing(f2))) {
    if (!directed) {
      stop("Node feature based assortativity coefficients are defined for directed networks.")
    }
    return(dw_feature_assort(
      netwk,
      f1 = f1, f2 = f2
    ))
  }

  if (!directed) {
    edgelist <- rbind(edgelist, edgelist[, c(2, 1)])
    edgeweight <- c(edgeweight, edgeweight)
  }

  snode <- edgelist[, 1]
  tnode <- edgelist[, 2]
  temp <- node_strength_cpp(
    snode = snode,
    tnode = tnode,
    nnode = nnode,
    weight = edgeweight,
    weighted = TRUE
  )
  outs <- temp$outs
  ins <- temp$ins
  rm(temp)

  if (!missing(weighted.rank)) {
    if (!directed) {
      stop("Rank-based assortativity coefficients are defined for directed networks.")
    }
    if (weighted.rank) {
      f1 <- weighted_rank(outs, outs)
      f2 <- weighted_rank(ins, ins)
    } else {
      f1 <- rank(outs, ties.method = "average")
      f2 <- rank(ins, ties.method = "average")
    }
    # if (!directed) {
    #   return(dw_feature_assort(netwk, f1 = f1, f2 = f2)$"f1-f2")
    # }
    result <- dw_feature_assort(netwk, f1 = f1, f2 = f2)
    names(result) <- c("outout", "outin", "inout", "inin")
    return(result)
  }

  sout <- outs[snode]
  tin <- ins[tnode]
  if (!directed) {
    return(wdm::wdm(
      x = sout, y = tin,
      weights = edgeweight, method = "pearson"
    ))
  }
  sin <- ins[snode]
  tout <- outs[tnode]
  return(list(
    "outout" = wdm::wdm(
      x = sout, y = tout,
      weights = edgeweight, method = "pearson"
    ),
    "outin" = wdm::wdm(
      x = sout, y = tin,
      weights = edgeweight, method = "pearson"
    ),
    "inout" = wdm::wdm(
      x = sin, y = tout,
      weights = edgeweight, method = "pearson"
    ),
    "inin" = wdm::wdm(
      x = sin, y = tin,
      weights = edgeweight, method = "pearson"
    )
  ))
}

#' Feature based assortativity coefficient
#'
#' Node feature based assortativity coefficients for weighted and directed
#' networks.
#'
#' @param netwk A \code{wdnet} object that represents the network.
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
#' assortcoef(adj = adj, f1 = f1, f2 = f2)
#'
#' @keywords internal
#'
dw_feature_assort <- function(netwk, f1, f2) {
  nnode <- max(netwk$edgelist)
  snode <- netwk$edgelist[, 1]
  tnode <- netwk$edgelist[, 2]
  edgeweight <- netwk$edge.attr$weight
  if (is.null(f1)) {
    f1 <- netwk$node.attr$outs
  }
  if (is.null(f2)) {
    f2 <- netwk$node.attr$ins
  }
  stopifnot(
    'Length of "f1" must equal number of nodes.' =
      length(f1) == nnode
  )
  stopifnot(
    'Length of "f2" must equal number of nodes.' =
      length(f2) == nnode
  )
  sf1 <- f1[snode]
  sf2 <- f2[snode]
  tf1 <- f1[tnode]
  tf2 <- f2[tnode]
  ret <- list()
  ret$"f1-f1" <- wdm::wdm(
    x = sf1, y = tf1,
    weights = edgeweight, method = "pearson"
  )
  ret$"f1-f2" <- wdm::wdm(
    x = sf1, y = tf2,
    weights = edgeweight, method = "pearson"
  )
  ret$"f2-f1" <- wdm::wdm(
    x = sf2, y = tf1,
    weights = edgeweight, method = "pearson"
  )
  ret$"f2-f2" <- wdm::wdm(
    x = sf2, y = tf2,
    weights = edgeweight, method = "pearson"
  )
  return(ret)
}

#' Weighted rank
#'
#' Computes weighted ranks for a numeric vector using weighted midranks for
#' tied observations. Missing values in either \code{x} or \code{weight} are
#' excluded from the calculation and returned as \code{NA} in the output.
#'
#' @param x A numeric vector to be ranked.
#'
#' @param weight A numeric vector of nonnegative weights with the same length as
#' \code{x}. Defaults to \code{x} (i.e. self-weighted ranks).
#'
#' @return A numeric vector of weighted ranks.
#'
#' @keywords internal

weighted_rank <- function(x, weight = x) {
  stopifnot(length(x) == length(weight))
  n <- length(x)
  ## identify valid observations (exclude NA in x or weight)
  ok <- !(is.na(x) | is.na(weight))
  ## initialize output vector with NA
  result <- rep(NA_real_, n)
  ## if all values are NA return NA vector
  if (!any(ok)) {
    return(result)
  }
  ## keep only valid observations
  x.valid <- x[ok]
  w.valid <- weight[ok]
  ## order observations by x (ascending rank order)
  ord <- order(x.valid)
  x.sorted <- x.valid[ord]
  w.sorted <- w.valid[ord]
  ## identify tie groups of identical x values
  rl <- rle(x.sorted)
  ## number of elements in each tie group
  g.size <- rl$lengths
  ## value of x for each tie group
  g.val <- rl$values
  ## group id for each observation in sorted order
  g.id <- rep.int(seq_along(g.size), g.size)
  ## sum of weights within each tie group, corresponds to sum_{k: s_k = s_i} w_k
  if (identical(w.valid, x.valid)) {
    g.sum <- g.val * g.size
  } else {
    g.sum <- rowsum(as.numeric(w.sorted), g.id, reorder = FALSE)[,1]
  }
  ## cumulative weight of strictly lower groups, corresponds to sum_{k: s_k < s_i} w_k
  before <- c(0, cumsum(g.sum))[seq_along(g.sum)]
  ## weighted midrank for tied observations, analogous to (t_i + 2)/2 * s_i in Eq. (2.2)
  mid <- (g.size + 1) / (2 * g.size) * g.sum
  ## weighted rank value for each group
  g.rank <- before + mid
  ## expand group ranks back to individual observations
  rank.sorted <- rep.int(g.rank, g.size)
  ## restore original order
  rank.ok <- numeric(sum(ok))
  rank.ok[ord] <- rank.sorted
  ## place results into full output vector
  result[ok] <- rank.ok
  return(result)
}
