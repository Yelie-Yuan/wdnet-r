##
## wdnet: Weighted directed network
## Copyright (C) 2023  Yelie Yuan, Tiandong Wang, Jun Yan and Panpan Zhang
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

#' @importFrom Matrix Matrix t colSums rowSums diag
NULL

#' Directed clustering coefficient
#'
#' Compute the clustering coefficient of a weighted and directed network.
#'
#' @param netwk A \code{wdnet} object that represents the network. If
#'   \code{NULL}, the function will compute the coefficient using either
#'   \code{edgelist}, \code{edgeweight}, or \code{adj}.
#' @param edgelist  A two-column matrix, each row represents a directed edge of
#'   the network.
#' @param edgeweight A vector representing the weight of edges.
#' @param adj An adjacency matrix of a weighted and directed network.
#' @param directed Logical. Indicates whether the edges in \code{edgelist} or
#'   \code{adj} are directed.
#' @param method Which method used to compute clustering coefficients: Clemente
#'   and Grassi (2018) or Fagiolo (2007).
#' @param isolates Binary, defines how to treat vertices with degree zero and
#'   one. If 0, then their clustering coefficient is returned as 0 and are
#'   included in the averaging. Otherwise, their clustering coefficient is \code{NaN}
#'   and are excluded in the averaging. Default value is 0.
#'
#' @return Lists of local clustering coefficients (in terms of a vector), global
#'   clustering coefficient (in terms of a scalar) and number of weighted
#'   directed triangles (in terms of a vector) based on \code{total}, \code{in},
#'   \code{out}, middleman (\code{middle}), or \code{cycle} triplets.
#'
#' @references
#' \itemize{
#'   \item Barrat, A., Barth\'{e}lemy, M., Pastor-Satorras,
#'   R. and Vespignani, A. (2004). The architecture of complex weighted
#'   networks. \emph{Proceedings of National Academy of Sciences of the United
#'   States of America}, 101(11), 3747--3752.
#'   \item Clemente, G.P. and Grassi,
#'   R. (2018). Directed clustering in weighted networks: A new perspective.
#'   \emph{Chaos, Solitons & Fractals}, 107, 26--38.
#'   \item Fagiolo, G. (2007).
#'   Clustering in complex directed networks. \emph{Physical Review E}, 76,
#'   026107.
#' }
#'
#' @note Self-loops (if exist) are removed prior to the computation of
#'   clustering coefficient. When the adjacency matrix is symmetric (i.e.,
#'   undirected but possibly unweighted networks), \code{clustcoef} returns
#'   local and global clustering coefficients proposed by Barrat et al. (2010).
#'
#' @examples
#' ## Generate a network according to the Erd\"{o}s-Renyi model of order 20
#' ## and parameter p = 0.3
#' edge_ER <- rbinom(400, 1, 0.3)
#' weight_ER <- sapply(edge_ER, function(x) x * sample(3, 1))
#' adj_ER <- matrix(weight_ER, 20, 20)
#' mycc <- clustcoef(adj = adj_ER, method = "Clemente")
#' system.time(mycc)
#'
#' @export
#' 

clustcoef <- function(
    netwk,
    edgelist,
    edgeweight,
    adj,
    directed = TRUE,
    method = c("Clemente", "Fagiolo"),
    isolates = 0) {
  if (missing(adj)) {
    netwk <- create_wdnet(
      netwk = netwk,
      edgelist = edgelist,
      edgeweight = edgeweight,
      directed = directed
    )
    # stopifnot(
    #   "Network must be directed." = netwk$directed
    # )
    adj <- edgelist_to_adj(
      edgelist = netwk$edgelist,
      edgeweight = netwk$edge.attr$weight,
      directed = netwk$directed
    )
  }
  stopifnot(dim(adj)[1] == dim(adj)[2])
  method <- match.arg(method)
  ## Force to remove self-loops.
  diag(adj) <- 0
  ## Extract the unweighted adjacency matrix
  adj <- Matrix::Matrix(adj, sparse = TRUE)
  A <- Matrix::Matrix(adj > 0, sparse = TRUE)
  ## Compute strength vector
  s_in <- Matrix::colSums(adj)
  s_out <- Matrix::rowSums(adj)
  s_tot <- s_in + s_out
  s_bil <- (Matrix::colSums(Matrix::t(adj) * A) +
    Matrix::colSums(Matrix::t(A) * adj)) / 2
  ## Compute degree vector
  d_in <- Matrix::colSums(A)
  d_out <- Matrix::rowSums(A)
  d_tot <- d_in + d_out
  A_A <- A %*% A
  d_bil <- Matrix::diag(A_A)
  if (method == "Clemente") {
    A_At <- A %*% Matrix::t(A)
    At_A <- Matrix::t(A) %*% A
    W_A_A <- Matrix::colSums(Matrix::t(adj) * A_A)
    W_A_At <- Matrix::colSums(Matrix::t(adj) * A_At)
    W_At_At <- Matrix::colSums(Matrix::t(adj) * Matrix::t(A_A))
    Wt_A_A <- Matrix::colSums(adj * A_A)
    Wt_A_At <- Matrix::colSums(adj * A_At)
    Wt_At_At <- Matrix::colSums(adj * Matrix::t(A_A))
    Wt_At_A <- Matrix::colSums(adj * At_A)
    W_At_A <- Matrix::colSums(Matrix::t(adj) * At_A)

    denomTotal <- s_tot * (d_tot - 1) - 2 * s_bil
    denomIn <- s_in * (d_in - 1)
    denomOut <- s_out * (d_out - 1)
    denomMiddle <- (s_in * d_out + s_out * d_in) / 2 - s_bil
    numTriangles <- list(
      "total" = (W_A_A + W_A_At + W_At_A + W_At_At +
        Wt_A_A + Wt_A_At + Wt_At_A + Wt_At_At) / 2,
      "in" = (Wt_A_A + Wt_At_A) / 2,
      "out" = (W_A_At + W_At_At) / 2,
      "middle" = (Wt_A_At + W_At_A) / 2,
      "cycle" = (W_A_A + Wt_At_At) / 2
    )
  }
  if (method == "Fagiolo") {
    # W is adjhat
    W <- (adj / max(adj))^(1 / 3)
    W_W <- W %*% W
    W_W_W <- Matrix::colSums(Matrix::t(W) * W_W)
    W_W_Wt <- Matrix::rowSums(W_W * W)
    W_Wt_W <- Matrix::colSums(Matrix::t(W) * (Matrix::t(W) %*% W))
    Wt_W_W <- Matrix::colSums(W * W_W)

    denomTotal <- d_tot * (d_tot - 1) - 2 * d_bil
    denomIn <- d_in * (d_in - 1)
    denomOut <- d_out * (d_out - 1)
    denomMiddle <- d_in * d_out - d_bil

    numTriangles <- list(
      "total" = (W_W_W + W_W_Wt + W_Wt_W + Wt_W_W),
      "in" = Wt_W_W,
      "out" = W_W_Wt,
      "middle" = W_Wt_W,
      "cycle" = W_W_W
    )
  }
  localcc <- list(
    "total" = numTriangles$"total" / denomTotal,
    "in" = numTriangles$"in" / denomIn,
    "out" = numTriangles$"out" / denomOut,
    "middle" = numTriangles$"middle" / denomMiddle,
    "cycle" = numTriangles$"cycle" / denomMiddle
  )
  if (isolates == "zero") {
    cat('Argument "isolates" has been revised; use "isolates = 0" instead.\n"')
    isolates <- 0
  }
  if (isolates == 0) {
    localcc <- rapply(localcc, function(i) ifelse(is.na(i), 0, i),
      how = "replace"
    )
  }
  globalcc <- list(
    "total" = mean(localcc$"total", na.rm = TRUE),
    "in" = mean(localcc$"in", na.rm = TRUE),
    "out" = mean(localcc$"out", na.rm = TRUE),
    "middle" = mean(localcc$"middle", na.rm = TRUE),
    "cycle" = mean(localcc$"cycle", na.rm = TRUE)
  )
  return(list(
    "total" = list(
      "localcc" = localcc$"total",
      "globalcc" = globalcc$"total",
      "numtriangles" = numTriangles$"total"
    ),
    "out" = list(
      "localcc" = localcc$"out",
      "globalcc" = globalcc$"out",
      "numtriangles" = numTriangles$"out"
    ),
    "in" = list(
      "localcc" = localcc$"in",
      "globalcc" = globalcc$"in",
      "numtriangles" = numTriangles$"in"
    ),
    "middle" = list(
      "localcc" = localcc$"middle",
      "globalcc" = globalcc$"middle",
      "numtriangles" = numTriangles$"middle"
    ),
    "cycle" = list(
      "localcc" = localcc$"cycle",
      "globalcc" = globalcc$"cycle",
      "numtriangles" = numTriangles$"cycle"
    )
  ))
}
