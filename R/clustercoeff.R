##
## wdnet: Weighted directed network
## Copyright (C) 2020  Panpan Zhang and Jun Yan
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

#' @importFrom matrixcalc matrix.power
#' @importFrom stats na.omit
NULL

#' Directed clustering coefficient
#' 
#' Compute the clustering coefficient of a weighted and directed network.
#'
#' @usage
#' dw_clustcoeff(adj, method = c("Clemente","Fagiolo"),
#'               mode = c("total","in","out","middle","cycle"))
#' 
#'
#' @param adj is an adjacency matrix of an weighted and directed network.
#' @param method which method used to compute clustering coefficients:
#' Clemente and Grassi (2018) or Fagiolo (2007)
#' @param mode what kind of triplets considered for computation:
#' \code{total}, \code{in}, \code{out}, 
#' middleman (\code{middle}), or \code{cycle}?
#'
#' @return a list of local clustering coefficients (in terms of a vector) and global clustering
#' coefficient (in terms of a scalar)
#'
#' @references
#' \itemize{
#' \item Barrat, A., Barth\'{e}lemy, M., Pastor-Satorras, R. and Vespignani, A. (2004). The 
#' architecture of complex weighted networks. \emph{Proceddings of National Academy of Sciences
#' of the United States of America}, 101(11), 3747--3752. 
#' \item Clemente, G.P. and Grassi, R. (2018). Directed clustering in weighted networks: 
#' A new perspective. \emph{Chaos, Solitons & Fractals}, 107, 26--38.
#' \item Fagiolo, G. (2007). Clustering in complex directed networks. \emph{Physical Review E},
#' 76, 026107.
#' }
#'
#' @note 
#' Self-loops (if exist) are removed prior to the computation of clustering oefficient.
#' When the adjacency matrix is symmetric (i.e., undirected but possibly unweighted networks), 
#' \code{dw_clustcoeff} returns local and global clustering coefficients proposedy by 
#' Barrat et al. (2010).
#'
#' @examples
#' ## Generate a network according to the Erd\"{o}s-Renyi model of order 20
#' ## and parameter p = 0.3
#' edge_ER <- rbinom(400,1,0.3)
#' weight_ER <- sapply(edge_ER, function(x) x*sample(3,1))
#' adj_ER <- matrix(weight_ER,20,20)
#' mycc_total <- dw_clustcoeff(adj_ER, method = "Clemente", mode = "total")
#' system.time(mycc_total)
#' 
#' @export

dw_clustcoeff <- function(adj, method = c("Clemente","Fagiolo"), 
                          mode = c("total","in","out","middle","cycle")){
  stopifnot(dim(adj)[1] == dim(adj)[2])
  ## Force to remove self-loops.
  if (all(diag(adj) == 0) == FALSE){
    n <- dim(adj)[1]
    adj <- adj - diag(diag(adj),n)
  }
  ## Extract the unweighted adjacency matrix
  uw_adj <- adj
  uw_adj[uw_adj > 0] <- 1
  ## Compute strength vector
  s_in <- colSums(adj)
  s_out <- rowSums(adj)
  s_tot <- s_in + s_out
  s_bil <- diag(adj %*% uw_adj + uw_adj %*% adj)/2
  ## Compute degee vector
  d_in <- colSums(uw_adj)
  d_out <- rowSums(uw_adj)
  d_tot <- d_in + d_out
  d_bil <- diag(matrixcalc::matrix.power(uw_adj, 2))
  if (method == "Clemente"){
    if (mode == "total"){
      local_cc <- diag((adj + t(adj)) %*% matrixcalc::matrix.power((uw_adj + t(uw_adj)), 2))/
        2/(s_tot*(d_tot - 1) - 2*s_bil)
      global_cc <- mean(stats::na.omit(local_cc))
    }
    if (mode == "in"){
      local_cc <- diag(t(adj) %*% (uw_adj + t(uw_adj)) %*% uw_adj)/
        2/(s_in*(d_in - 1))
      global_cc <- mean(na.omit(local_cc))
    }
    if (mode == "out"){
      local_cc <- diag(adj %*% (uw_adj + t(uw_adj)) %*% t(uw_adj))/
        2/(s_out*(d_out - 1))
      global_cc <- mean(na.omit(local_cc))
    }
    if (mode == "middle"){
      local_cc <- diag(t(adj) %*% uw_adj %*% t(uw_adj) + adj %*% t(uw_adj) %*% uw_adj)/
        2/((s_in*d_out + s_out*d_in)/2 - s_bil)
      global_cc <- mean(na.omit(local_cc))
    }
    if (mode == "cycle"){
      local_cc <- diag(adj %*% matrixcalc::matrix.power(uw_adj, 2) 
                       + t(adj) %*% matrixcalc::matrix.power(t(uw_adj), 2))/
        2/((s_in*d_out + s_out*d_in)/2 - s_bil)
      global_cc <- mean(na.omit(local_cc))
    }
  }
  if (method == "Fagiolo"){
    adjhat <- (adj/max(adj))^(1/3)
    if (mode == "total"){
      local_cc <- diag(matrixcalc::matrix.power(adjhat + t(adjhat), 3))/
        2/(d_tot*(d_tot - 1) - 2*d_bil)
      global_cc <- mean(na.omit(local_cc))
    }
    if (mode == "in"){
      local_cc <- diag(t(adjhat) %*% matrixcalc::matrix.power(adjhat, 2))/
        (d_in*(d_in - 1))
      global_cc <- mean(na.omit(local_cc))
    }
    if (mode == "out"){
      local_cc <- diag(matrixcalc::matrix.power(adjhat, 2) %*% t(adjhat))/
        (d_out*(d_out - 1))
      global_cc <- mean(na.omit(local_cc))
    }
    if (mode == "middle"){
      local_cc <- diag(adjhat %*% t(adjhat) %*% adjhat)/
        (d_in*d_out - d_bil)
      global_cc <- mean(na.omit(local_cc))
    }
    if (mode == "cycle"){
      local_cc <- diag(matrixcalc::matrix.power(adjhat, 3))/
        (d_in*d_out - d_bil)
      global_cc <- mean(na.omit(local_cc))
    }
  }
  return(list(localcc = local_cc, globalcc = global_cc))
}
