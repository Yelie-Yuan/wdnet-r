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

#' Simulating a Preferential Attachment Network
#'
#' @param alpha Scalar probability of adding an edge from the new node
#'     to an existing node
#' @param beta Scalar probability of adding an edge between two
#'     existing nodes.
#' @param gamma Scalar probability of adding an edge from an existing
#'     node to a new node.
#' @param xi Scalar probability of ...
#' @param delta_in Growth rate parameter for nodes' instrength
#' @param delta_out Growth rate parameter for nodes' outstrength
#' @param nedge The number of edges to be generated
#' 
#' @return A list with the following components: in_degree,
#'     out_degree, edge_start, edge_end, evolution
#'
#' @keywords internal
#' 


rpanet_wan <- function(alpha, beta, gamma, xi,
                       delta_in, delta_out, nedge) {
    in_degree <- integer(nedge * 2)
    out_degree <- integer(nedge * 2)
    edge_start <- integer(nedge)
    edge_end <- integer(nedge)
    evolution <- integer(nedge)
    ret <- .C("netSim",
              as.double(alpha), as.double(beta), as.double(gamma),
              as.double(xi), as.double(delta_in),as.double(delta_out),
              as.integer(nedge),
              ## output
              ind  = as.integer(in_degree), outd = as.integer(out_degree),
              start = as.integer(edge_start), end = as.integer(edge_end),
              evol = as.integer(evolution),
              PACKAGE = "wdnet")
    list(in_degree = ret$ind, out_degree = ret$outd,
         edge_start = ret$start, edge_end = ret$end,
         evolution = ret$evol)
}
