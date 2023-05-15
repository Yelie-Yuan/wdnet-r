##
## wdnet: Weighted directed network
## Copyright (C) 2023  Yelie Yuan, Tiandong Wang, Jun Yan and Panpan Zhang
## Yelie Yuan <yelie.yuan@uconn.edu>
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

#' wdnet: Weighted and Directed Networks
#'
#' This package provides functions to conduct network analysis
#' \itemize{
#' \item Assortativity, centrality, clustering coefficient
#' for weighted and directed networks
#' \item Rewire an unweighted network with given assortativity coefficient(s)
#' \item Preferential attachment (PA) network generation
#' }
#'
#' @section wdnet networks: wdnet networks have a class \code{wdnet}. It is a
#'   list containing the following components:
#' \itemize{
#' \item A logical value \code{directed} indicating if the network is directed.
#' \item A logical value \code{weighted} indicating if the network is weighted.
#' \item A two-column matrix representing the edges.
#' \item A data frame \code{node.attr} that includes node attributes,
#' such as node strengths.
#' \item A data frame \code{edge.attr} that includes edge attributes,
#' such as edge weights.
#' }
#'
#' @section Creating a \code{wdnet} Object:
#' \itemize{
#' \item To generate a preferential attachment (PA) network,
#' use \code{rpanet()}.
#' \item To create a \code{wdnet} object from an edge list
#' and edge weights, use \code{edgelist_to_wdnet()}.
#' \item To create a \code{wdnet} object from an adjacency
#' matrix, use \code{adj_to_wdnet()}.
#' \item To convert an \code{igraph} object to a \code{wdnet}
#' object, use \code{igraph_to_wdnet()}.
#' }
#'
#' @section Further information: The development version of this package is
#'   available on Gitlab (\url{https://gitlab.com/wdnetwork/wdnet}).
#'
#' @docType package
#' @name wdnet
#' @useDynLib wdnet, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp

NULL
