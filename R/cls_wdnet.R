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

#' @importFrom igraph graph_from_edgelist E V plot.igraph is.igraph as_edgelist
#' @importFrom igraph is.directed vertex.attributes edge.attributes
#' @importFrom utils head
NULL

#' Creates a \code{wdnet} object using \code{edgelist}.
#'
#' @param edgelist A two-column matrix representing the edges.
#' @param edgeweight A numeric vector of edge weights with the same length as
#'   the number of rows in edgelist. If \code{NULL}, all edges will be assigned
#'   weight 1.
#' @param directed Logical, whether the network is directed (TRUE) or undirected
#'   (FALSE).
#' @param nodegroup A numeric vector of node groups.
#' @param ... Additional components to be added to the \code{wdnet} object.
#' @return A \code{wdnet} object with the specified \code{edgelist},
#'   \code{edgeweight} and \code{directed}.
#' @examples
#' edgelist <- matrix(c(1, 2, 2, 3, 3, 1), ncol = 2, byrow = TRUE)
#' edgeweight <- c(1, 2, 3)
#' nodegroup <- c(1, 1, 2)
#' netwk <- edgelist_to_wdnet(
#'   edgelist = edgelist,
#'   edgeweight = edgeweight,
#'   directed = TRUE,
#'   nodegroup = nodegroup
#' )
#'
#' @export
#'
edgelist_to_wdnet <- function(
    edgelist,
    edgeweight,
    directed,
    nodegroup,
    ...) {
  if (missing(directed) || is.null(directed)) {
    # cat("Assume the network is directed.\n\n")
    directed <- TRUE
  }
  stopifnot(is.logical(directed))
  if (missing(edgelist) || is.null(edgelist)) {
    stop('Please provide "edgelist".')
  }

  if (missing(edgeweight) || is.null(edgeweight)) {
    edgeweight <- rep(1, nrow(edgelist))
  }

  if (nrow(edgelist) != length(edgeweight)) {
    stop("Number of rows in 'edgelist' must match the length of 'edgeweight'.")
  }

  if (ncol(edgelist) != 2) {
    stop('"edgelist" must have exactly 2 columns.')
  }

  if (missing(nodegroup)) {
    nodegroup <- NULL
  }

  max_node <- max(edgelist)

  if (!is.null(nodegroup)) {
    if (length(nodegroup) != max_node) {
      stop('Length of "nodegroup" must match the number of nodes in "edgelist".')
    }
  }
  mode(edgelist) <- "integer"

  weighted <- any(edgeweight != 1)
  netwk <- structure(
    list(
      "edgelist" = edgelist,
      "directed" = directed,
      "weighted" = weighted,
      "edge.attr" = data.frame("weight" = edgeweight)
    ),
    class = "wdnet"
  )

  tmp <- node_strength_cpp(
    snode = edgelist[, 1],
    tnode = edgelist[, 2],
    weight = edgeweight,
    nnode = max_node,
    weighted = weighted
  )

  if (directed) {
    netwk$node.attr <- data.frame(
      "outs" = tmp$outs,
      "ins" = tmp$ins
    )
  } else {
    tmp$s <- tmp$outs + tmp$ins
    netwk$node.attr <- data.frame(
      "s" = tmp$s
    )
  }
  netwk$node.attr$group <- nodegroup

  additional_components <- list(...)
  if (length(additional_components) > 0) {
    netwk <- c(netwk, additional_components)
    class(netwk) <- "wdnet"
  }

  if (!is_wdnet(netwk)) {
    stop('Failed to create a valid "wdnet" object.')
  }

  return(netwk)
}

#' Creates a \code{wdnet} object using an adjacency matrix
#'
#' @param adj An adjacency matrix used to extract \code{edgelist} and
#'   \code{edgeweight} using \code{igraph}.
#' @param directed Logical, whether the network is directed (TRUE) or undirected
#'   (FALSE). If \code{adj} is asymmetric, the network is directed.
#' @param weighted Logical, whether the network is weighted (TRUE) or unweighted
#'   (FALSE).
#' @param nodegroup A numeric vector of node groups.
#' @param ... Additional components to be added to the \code{wdnet} object.
#' @return A \code{wdnet} object with the specified \code{adj}.
#' @export
#' @examples
#' adj <- matrix(c(0, 1, 2, 0), nrow = 2, ncol = 2, byrow = TRUE)
#' adj_to_wdnet(adj = adj, directed = TRUE, weighted = FALSE)
#'
adj_to_wdnet <- function(
    adj,
    directed = TRUE,
    weighted = TRUE,
    nodegroup,
    ...) {
  if (missing(adj) || is.null(adj)) {
    stop('Please provide "adj".')
  }
  if (missing(directed) || is.null(directed)) {
    # cat("Assume the network is directed.\n\n")
    directed <- TRUE
  } else if (!directed) {
    if (!isSymmetric(adj)) {
      directed <- TRUE
      cat('Network is directed because "adj" is asymmetric.\n\n')
    }
  }

  if (missing(weighted) || is.null(weighted)) {
    weighted <- any(adj[adj > 0] != 1)
    if (weighted) {
      cat("Assume the network is weighted.\n")
    }
  }

  stopifnot(is.logical(directed))
  stopifnot(is.logical(weighted))

  if (missing(nodegroup)) {
    nodegroup <- NULL
  }

  if (!is.null(nodegroup)) {
    if (length(nodegroup) != nrow(adj)) {
      stop('Length of "nodegroup" must match the number of nodes in "adj".')
    }
  }

  tmp <- adj_to_edgelist(
    adj = adj,
    directed = directed,
    weighted = weighted
  )
  edgelist <- tmp$edgelist
  edgeweight <- tmp$edgeweight
  directed <- tmp$directed
  rm(tmp)

  edgelist_to_wdnet(
    edgelist = edgelist,
    edgeweight = edgeweight,
    directed = directed,
    nodegroup = nodegroup,
    ...
  )
}


#' Creates a \code{wdnet} object from input data.
#'
#' This function creates a \code{wdnet} object from \code{edgelist} and
#' \code{edgeweight} or \code{adj} or returns the existing \code{wdnet} object.
#' For internal usage.
#'
#' @param netwk A \code{wdnet} object. If \code{NULL}, the function will use the
#'   provided \code{edgelist} and \code{edgeweight}, or \code{adj} parameters to
#'   create a new \code{wdnet} object.
#' @param edgelist A two-column matrix representing edges.
#' @param edgeweight A vector representing the weights of the edges.
#' @param nodegroup A numeric vector of node groups.
#' @param directed A logical value indicating whether the network is directed.
#'   Required if \code{netwk} is \code{NULL}.
#' @param adj An adjacency matrix.
#' @param weighted A logical value indicating whether the network is weighted.
#' @param ... Additional components to be added to the wdnet list.
#'
#' @return A \code{wdnet} object.
#'
#' @keywords internal
#'
create_wdnet <- function(
    netwk,
    edgelist,
    edgeweight,
    nodegroup,
    directed,
    adj,
    weighted,
    ...) {
  if (missing(netwk) || is.null(netwk) || !is_wdnet(netwk)) {
    if (missing(edgelist) || is.null(edgelist)) {
      if (missing(adj) || is.null(adj)) {
        stop('Please provide either "edgelist" or "adj".')
      } else {
        netwk <- adj_to_wdnet(
          adj = adj,
          directed = directed,
          weighted = weighted,
          nodegroup = nodegroup,
          ...
        )
      }
    } else {
      colnames(edgelist) <- NULL
      netwk <- edgelist_to_wdnet(
        edgelist = edgelist,
        edgeweight = edgeweight,
        directed = directed,
        nodegroup = nodegroup,
        ...
      )
    }
  } else {
    if (!missing(directed)) {
      if (directed != netwk$directed) {
        cat(
          'The "directed" argument is omitted since "netwk" is ', ifelse(
            netwk$directed, "directed.", "undirected."
          ), "\n"
        )
      }
    }
    additional_components <- list(...)
    if (length(additional_components) > 0) {
      netwk <- c(netwk, additional_components)
      class(netwk) <- "wdnet"
    }
    if (!is.null(netwk$nodegroup)) {
      netwk$node.attr$group <- netwk$nodegroup
      netwk$nodegroup <- NULL
    }
  }
  invisible(netwk)
}

#' Checks if the input is a \code{wdnet} object
#'
#' @param netwk A \code{wdnet} object.
#' @return Logical, \code{TRUE} if argument netwk is a \code{wdnet} object.
#' @export
#' @examples
#' netwk <- rpanet(nstep = 1e3)
#' is_wdnet(netwk)
#'
is_wdnet <- function(netwk) {
  valid_attrs <- c("edgelist", "directed", "weighted", "edge.attr", "node.attr")
  valid_cols <- ifelse(netwk$directed, c("outs", "ins"), c("s"))
  max_node <- max(netwk$edgelist)
  return(all(
    is.logical(netwk$directed),
    is.logical(netwk$weighted),
    ifelse(
      netwk$weighted,
      netwk$edge.attr$weight > 0,
      all(netwk$edge.attr$weight == 1)
    ),
    inherits(netwk, "wdnet"),
    "weight" %in% names(netwk$edge.attr),
    all(valid_attrs %in% names(netwk)),
    all(valid_cols %in% colnames(netwk$node.attr)),
    nrow(netwk$edgelist) == nrow(netwk$edge.attr),
    is.integer(netwk$edgelist),
    max_node == nrow(netwk$node.attr)
  ))
}

#' Converts a \code{wdnet} object to an \code{igraph} object
#'
#' @param netwk A \code{wdnet} object.
#'
#' @return An \code{igraph} object.
#'
#' @export
#' @examples
#' netwk <- rpanet(nstep = 1e3)
#' g <- wdnet_to_igraph(netwk)
#'
wdnet_to_igraph <- function(netwk) {
  stopifnot(is_wdnet(netwk))
  g <- igraph::graph_from_edgelist(netwk$edgelist,
    directed = netwk$directed
  )
  if (is.data.frame(netwk$edge.attr)) {
    for (each in colnames(netwk$edge.attr)) {
      g <- igraph::set_edge_attr(
        g, each,
        value = netwk$edge.attr[[each]]
      )
    }
  }
  if (is.data.frame(netwk$node.attr)) {
    for (each in colnames(netwk$node.attr)) {
      g <- igraph::set_vertex_attr(
        g, each,
        value = netwk$node.attr[[each]]
      )
    }
  }
  return(g)
}

#' Converts an \code{igraph} object to a \code{wdnet} object
#'
#' @param g An \code{igraph} object.
#'
#' @return A \code{wdnet} object.
#'
#' @export
#' @examples
#' g <- igraph::sample_pa(50)
#' netwk <- igraph_to_wdnet(g)
#'
igraph_to_wdnet <- function(g) {
  stopifnot(igraph::is.igraph(g))

  edgelist <- igraph::as_edgelist(g, names = FALSE)
  mode(edgelist) <- "integer"
  edgeweight <- igraph::E(g)$weight
  directed <- igraph::is.directed(g)

  netwk <- create_wdnet(
    edgelist = edgelist,
    edgeweight = edgeweight,
    directed = directed
  )

  nattr <- igraph::vertex.attributes(g)
  if (length(nattr) > 0) {
    nattr <- as.data.frame(nattr)
    for (each in colnames(nattr)) {
      if (each %in% colnames(netwk$node.attr)) {
        next
      }
      netwk$node.attr[[each]] <- nattr[[each]]
    }
  }

  eattr <- igraph::edge.attributes(g)
  if (length(eattr) > 0) {
    eattr <- as.data.frame(eattr)
    for (each in colnames(eattr)) {
      if (each %in% colnames(netwk$edge.attr)) {
        next
      }
      netwk$edge.attr[[each]] <- eattr[[each]]
    }
  }

  is_wdnet(netwk)
  netwk
}

#' Plots the input network
#'
#' Plots the input network via \code{igraph::plot.igraph()}.
#'
#' @param x A \code{wdnet} object.
#' @param ... Additional parameters passed to \code{igraph::plot.igraph()}.
#' @return Returns \code{NULL}, invisibly.
#' @method plot wdnet
#' @export
#'
plot.wdnet <- function(x, ...) {
  stopifnot(is_wdnet(x))
  igraph::plot.igraph(wdnet_to_igraph(x), ...)
  invisible(NULL)
}

#' Prints the input network
#'
#' These functions print a network to the terminal.
#'
#' \code{summary.wdnet} prints the number of nodes and edges, preference
#' functions, and whether the network is directed, weighted. \code{print.wdnet}
#' prints the same information, and also lists some edges and node attributes,
#' if available. Edge scenarios are 0: from initial network; 1: \code{alpha}; 2:
#' \code{beta}; 3: \code{gamma}; 4: \code{xi}; 5; \code{rho}; 6: reciprocal.
#'
#' @param x A \code{wdnet} object.
#' @param node.attrs Logical, whether to print node attributes, if available.
#' @param edge.attrs Logical, whether to print edge attributes, if available.
#' @param max.lines Integer, the maximum number of lines of edgelist and node
#'   attributes to print. The rest of the output will be truncated.
#' @param object The graph of which the summary will be printed.
#' @param ... Additional arguments.
#' @rdname print.wdnet
#' @method print wdnet
#' @export
#'
print.wdnet <- function(
    x,
    node.attrs = TRUE,
    edge.attrs = TRUE,
    max.lines = 5,
    ...) {
  summary.wdnet(x)
  nedge <- nrow(x$edgelist)
  nnode <- nrow(x$node.attr)
  n <- min(max.lines, nedge)

  cat("\nEdges:\n")
  tmp <- data.frame(x$edgelist[seq_len(n), ])
  if (x$directed) {
    colnames(tmp) <- c("source", "target")
  } else {
    colnames(tmp) <- c("i", "j")
  }
  if (edge.attrs) {
    if (is.null(x$edge.attr)) {
      cat("No edge attributes\n")
    } else {
      tmp <- cbind(tmp, x$edge.attr[seq_len(n), ])
      colnames(tmp)[3:ncol(tmp)] <- colnames(x$edge.attr)
    }
  }
  print(tmp)

  if (n < nedge) cat("...omitted remaining edges\n")

  n <- min(max.lines, nnode)
  cat("\nNode attributes:\n")
  if (node.attrs) {
    if (is.null(x$node.attr)) {
      cat("No available node attributes to print.")
    } else {
      print(utils::head(x$node.attr, n))
      if (n < max(x$edgelist)) cat("...omitted remaining nodes\n")
    }
  }
  invisible(x)
}


#' @rdname print.wdnet
#' @method summary wdnet
#' @export
#'
summary.wdnet <- function(object, ...) {
  stopifnot(is_wdnet(object))
  cat(
    "Weighted: ", object$weighted, "\n",
    "Directed: ", object$directed, "\n",
    "Number of edges: ", nrow(object$edge.attr), "\n",
    "Number of nodes: ", nrow(object$node.attr), "\n",
    sep = ""
  )

  invisible(object)
}
