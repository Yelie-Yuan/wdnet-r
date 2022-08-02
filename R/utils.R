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

#' @importFrom igraph graph_from_adjacency_matrix as_edgelist E
#' @importFrom Rcpp sourceCpp
NULL

#' Convert adjacency matrix to edgelist and edgeweight.
#'
#' @param adj Adjacency matrix of a network.
#' @param directed Logical, whether the network is directed. Passed to
#'   \code{igraph::graph_from_adjacency_matrix}.
#' @param weighted Passed to \code{igraph::graph_from_adjacency_matrix}. This
#'   argument specifies whether to create a weighted graph from an adjacency
#'   matrix. If it is NULL then an unweighted graph is created and the elements
#'   of the adjacency matrix gives the number of edges between the vertices. If
#'   it is TRUE then a weighted graph is created and the name of the edge
#'   attribute will be weight.
#'
#' @return A list of edgelist and edgeweight.
#'   
adj_to_edge <- function(adj, directed = TRUE, weighted = TRUE) {
  if (! directed) {
    stopifnot('"adj" must be symmetric if the network is undirected.' = 
                isSymmetric(adj))
  }
  mode <- ifelse(directed, "directed", "undirected")
  g <- igraph::graph_from_adjacency_matrix(adj, mode = mode, 
                                           weighted = weighted, diag = TRUE)
  edgelist <- igraph::as_edgelist(g)
  edgeweight <- igraph::E(g)$weight
  return(list("edgelist" = edgelist, 
              "edgeweight" = edgeweight))
}

#' Convert edgelist and edgeweight to adjacency matrix.
#'
#' @param edgelist A two column matrix represents edges.
#' @param edgeweight A vector represents the weight of edges. If \code{NULL},
#'   all the edges are considered have weight 1.
#' @param directed Logical, whether the network is directed.
#'
#' @return An adjacency matrix.
#'
edge_to_adj <- function(edgelist, edgeweight = NULL, directed = TRUE) {
  nnode <- max(edgelist)
  adj <- matrix(0, nrow = nnode, ncol = nnode)
  if (is.null(edgeweight)) {
    edgeweight <- rep(1, nrow(edgelist))
  }
  adj <- fill_weight_cpp(adj, edgelist - 1, edgeweight)
  if (! directed) {
    adj <- adj + t(adj)
    diag(adj) <- diag(adj) / 2
  }
  return(adj)
}


#' Compile preference functions via \code{Rcpp}.
#' 
#' @param preference A list for defining the preference functions.
#' @return External pointers of preference functions.
compile_pref_func <- function(preference) {
  compile_spref <- compile_tpref <- compile_pref <- FALSE
  cpp_code <- "// [[Rcpp::depends(RcppArmadillo)]]
    #include <RcppArmadillo.h>
    #include <math.h>
    using namespace Rcpp;
    typedef double (*funcPtrD)(double x, double y);"
  if (typeof(preference$spref) == "character") {
    compile_spref <- TRUE
    cpp_code <- paste(cpp_code, "\n",
                      "// [[Rcpp::export]]
      double spref_func(double outs, double ins) {
        return ",
                      preference$spref,
                      ";
      }
      // [[Rcpp::export]]
      XPtr<funcPtrD> put_spref_XPtr() {
        return(XPtr<funcPtrD>(new funcPtrD(spref_func)));
      }", collapse = "\n")
  }
  else if (typeof(preference$spref) != "externalptr") {
    stop('Type of "spref" must be "externalptr" or "character".')
  }
  
  if (typeof(preference$tpref) == "character") {
    compile_tpref <- TRUE
    cpp_code <- paste(cpp_code, "\n",
                      "// [[Rcpp::export]]
      double tpref_func(double outs, double ins) {
        return ",
                      preference$tpref,
                      ";
      }
      // [[Rcpp::export]]
      XPtr<funcPtrD> put_tpref_XPtr() {
        return(XPtr<funcPtrD>(new funcPtrD(tpref_func)));
      }", collapse = "\n")
  }
  else if (typeof(preference$tpref) != "externalptr") {
    stop('Type of "tpref" must be "externalptr" or "character".')
  }
  
  if (typeof(preference$pref) == "character") {
    compile_pref <- TRUE
    cpp_code <- paste(cpp_code, "\n",
                      "typedef double (*funcPtrUnd)(double x);
      // [[Rcpp::export]]
      double pref_func(double s) {
        return ",
                      preference$pref,
                      ";
      }
      // [[Rcpp::export]]
      XPtr<funcPtrUnd> put_pref_XPtr() {
        return(XPtr<funcPtrUnd>(new funcPtrUnd(pref_func)));
      }", collapse = "\n")
  }
  else if (typeof(preference$pref) != "externalptr") {
    stop('Type of "pref" must be "externalptr" or "character".')
  }
  
  if (compile_spref | compile_tpref | compile_pref) {
    cat("Compiling preference function(s)...\n")
    Rcpp::sourceCpp(code = cpp_code)
    # cat("Done.\n")
  }
  preference$spref.pointer <- ifelse(compile_spref, 
                                  put_spref_XPtr(),
                                  preference$spref)
  preference$tpref.pointer <- ifelse(compile_tpref, 
                                  put_tpref_XPtr(),
                                  preference$tpref)
  preference$pref.pointer <- ifelse(compile_pref, 
                                 put_pref_XPtr(),
                                 preference$pref)
  # test functions
  test_pref_func_directed(preference$spref.pointer, 1, 1)
  test_pref_func_directed(preference$tpref.pointer, 1, 1)
  test_pref_func_undirected(preference$pref.pointer, 1)
  preference
}
