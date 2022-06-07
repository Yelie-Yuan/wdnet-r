##
## wdnet: Weighted directed network
## Copyright (C) 2022  Yelie Yuan, Tiandong Wang, Jun Yan and Panpan Zhang
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

#' @importFrom utils modifyList
NULL

#' Add components to the control list
#'
#' `+` is used to combine components to control the PA network generation
#' process. Available components are \code{scenario.control()},
#' \code{edgeweight.control()}, \code{newedge.control()},
#' \code{preference.control()} and \code{reciprocal.control()}.
#'
#' @param e1 An object of class \code{panet.control}.
#' @param e2 An object of class \code{panet.control}.
#'
#' @export
#'
#' @examples
#' scenario.control(alpha = 0.5, beta = 0.5) +
#'     preference.control(sparams = c(1, 1, 0, 0, 1),
#'         tparams = c(0, 0, 1, 1, 1))
#'
#' scenario.control(alpha = 1) +
#'     edgeweight.control(distribution = rgamma,
#'         dparams = list(shape = 5, scale = 0.2),
#'         shift = 1)
"+.panet.control" <- function(e1, e2) {
  e1 <- structure(modifyList(e1, e2, keep.null = TRUE), class = "panet.control")
  if (is.list(e2$edgeweight$dparams)) {
    e1$edgeweight$dparams <- e2$edgeweight$dparams
  }
  if (is.list(e2$newedge$dparams)) {
    e1$newedge$dparams <- e2$newedge$dparams
  }
  e1
}

#' Set parameters for controlling the probability of edge scenarios
#'
#' @param alpha Probability of adding an edge from a new node to an existing
#'   node.
#' @param beta Probability of adding an edge between existing nodes.
#' @param gamma Probability of adding an edge from an existing node to a new
#'   node.
#' @param xi Probability of adding an edge between two new nodes. \code{rho = 1
#'   - alpha - beta - gamma - xi} represents the probability of introducing a
#'   new node with a self looped edge.
#' @param beta.loop Logical, whether self-loops are allowed under beta scenario.
#'   Default value is \code{TRUE}.
#' @param source.first Logical. Defined for \code{beta} scenario edges of
#'   directed networks. If \code{TRUE}, the source node of a new edge is sampled
#'   from existing nodes before the target node is sampled; if \code{FALSE}, the
#'   target node is sampled from existing nodes before the source node is
#'   sampled. Default value is \code{TRUE}.
#'
#' @export
#'
#' @examples
#' scenario.control(alpha = 0.5, beta = 0.5, beta.loop = FALSE)
#' 
scenario.control <- function(alpha = 1, beta = 0, gamma = 0, xi = 0,
                             beta.loop = TRUE, source.first = TRUE) {
  stopifnot('"alpha + beta + gamma + xi" must be smaller or equal to 1.' =
              alpha + beta + gamma + xi <= 1)
  scenario <- list("alpha" = alpha, 
                   "beta" = beta,
                   "gamma" = gamma, 
                   "xi" = xi,
                   "rho" = 1 - alpha - beta - gamma - xi,
                   "beta.loop" = beta.loop, 
                   "source.first" = source.first)
  structure(list("scenario" = scenario),
            class = "panet.control")
}

#' Set parameters for controlling weight of new edges
#'
#' @param distribution Distribution function for edge weights. Default is
#'   \code{NA}. If specified, its first argument must be the number of
#'   observations.
#' @param dparams Additional parameters passed on to \code{distribution}. The
#'   name of parameters must be specified.
#' @param shift A constant add to the specified distribution. Default value is
#'   1.
#'
#' @export
#'
#' @examples
#' # Edge weight follows Gamma(5, 0.2).
#' edgeweight.control(distribution = rgamma,
#'     dparams = list(shape = 5, scale = 0.2),
#'     shift = 0)
#'
#' # Constant edge weight
#' edgeweight.control(shift = 2)
#' 
edgeweight.control <- function(distribution = NA,
                               dparams = list(),
                               shift = 1) {
  edgeweight <- list("distribution" = distribution,
                     "dparams" = dparams,
                     "shift" = shift)
  if (length(edgeweight$dparams) > 0) {
    stopifnot("Please specify the name of distribution parameters." = 
                all(! is.null(names(edgeweight$dparams))))
    stopifnot("Please provide the distribution function." = 
                is.function(distribution))
  }
  structure(list("edgeweight" = edgeweight),
            class = "panet.control")
}

#' Set parameters for controlling new edges in each step
#'
#' @param distribution Distribution function for number of new edges. Default is
#'   \code{NA}. If specified, its first argument must be the number of
#'   observations.
#' @param dparams Additional parameters passed on to \code{distribution}. The
#'   name of parameters must be specified.
#' @param shift A constant add to the specified distribution. Default value is
#'   1.
#' @param snode.replace Logical, whether the source nodes in the same step
#'   should be sampled with replacement. Defined for directed networks.
#' @param tnode.replace Logical, whether the target nodes in the same step
#'   should be sampled with replacement. Defined for directed networks.
#' @param node.replace Logical, whether the nodes in the same step should be
#'   sampled with replacement. Defined for undirected and directed networks. For
#'   directed networks, when \code{node.replace} is \code{FALSE}, sampled
#'   source and target nodes in the same step are all different from each other,
#'   \code{beta.loop} will be set as \code{FALSE}.
#'
#' @export
#'
#' @examples
#' newedge.control(distribution = rpois,
#'     dparams = list(lambda = 2),
#'     shift = 1,
#'     node.replace = FALSE)
newedge.control <- function(distribution = NA,
                            dparams = list(),
                            shift = 1,
                            snode.replace = TRUE,
                            tnode.replace = TRUE,
                            node.replace = TRUE) {
  if (! node.replace) {
    snode.replace <- tnode.replace <- FALSE
  }
  newedge <- list("distribution" = distribution,
                  "dparams" = dparams,
                  "shift" = shift,
                  "snode.replace" = snode.replace,
                  "tnode.replace" = tnode.replace,
                  "node.replace" = node.replace)
  if (length(newedge$dparams) > 0) {
    stopifnot("Please specify the name of distribution parameters" = 
                all(! is.null(names(newedge$dparams))))
  }
  structure(list("newedge" = newedge), class = "panet.control")
}

#' Set parameters for source and target preference function
#'
#' @param sparams Parameters of the source preference function for directed
#'   networks. Probability of choosing an exising node as the source node is
#'   proportional to \code{sparams[1] * out-strength^sparams[2] + sparams[3] *
#'   in-strength^sparams[4] + sparams[5]}.
#' @param tparams Parameters of the target preference function for directed
#'   networks. Probability of choosing an exising node as the source node is
#'   proportional to \code{tparams[1] * out-strength^tparams[2] + tparams[3] *
#'   in-strength^tparams[4] + tparams[5]}.
#' @param params Parameters of the preference function for undirected networks.
#'   Probability of choosing an exising node is proportional to
#'   \code{strength^param[1] + param[2]}.
#'
#' @export
#'
#' @examples
#' preference.control(sparams = c(1, 2, 0, 0, 0.1),
#'     tparams = c(0, 0, 1, 2, 0.1))
preference.control <- function(sparams = c(1, 1, 0, 0, 1),
                               tparams = c(0, 0, 1, 1, 1),
                               params = c(1, 1)) {
  preference <- list("sparams" = sparams,
               "tparams" = tparams,
               "params" = params)
  structure(list("preference" = preference), 
            class = "panet.control")
}

#' Set parameters for controlling reciprocal edges
#'
#' @param group.prob A vector of probability weights for sampling the group of
#'   new nodes. Defined for directed networks. Groups are from 1 to
#'   \code{length(group.prob)}. Its length must equal to the number of rows of
#'   \code{recip.prob}.
#' @param recip.prob A square matrix giving the probability of adding a
#'   reciprocal edge after a new edge is introduced. Defined for directed
#'   networks. Its element \code{p_{ij}} represents the probability of adding a
#'   reciprocal edge from node \code{A}, which belongs to group \code{i}, to
#'   node \code{B}, which belongs to group \code{j}, immediately after a
#'   directed edge from \code{B} to \code{A} is added.
#' @param selfloop.recip Logical, whether reciprocal edge of self-loops are
#'   allowed.
#'
#' @export
#'
#' @examples
#' reciprocal.control(group.prob = c(0.4, 0.6),
#'     recip.prob = matrix(runif(4), ncol = 2))
reciprocal.control <- function(group.prob = NULL,
                               recip.prob = NULL, 
                               selfloop.recip = FALSE) {
  if (! is.null(group.prob)) {
    stopifnot('"group.prob" must sum to 1.' = 
                round(sum(group.prob), 10) == 1)
    if (! is.null(recip.prob)) {
      recip.prob <- as.matrix(recip.prob)
      stopifnot('"recip.prob" or "group.prob" is not valid.' =
                  length(group.prob) == nrow(recip.prob) &
                  nrow(recip.prob) == ncol(recip.prob))
      stopifnot('"recip.prob" is not valid.' =
                  all(recip.prob >= 0) &
                  all(recip.prob <= 1))
      stopifnot('"group.prob" is not valid.' =
                  round(sum(group.prob), 10) == 1 &
                  all(group.prob >= 0))
    }
    else {
      stop('"recip.prob" can not be NA when "group.prob" is specified.')
    }
  }
  if (is.null(group.prob)) {
    if (! is.null(recip.prob)) {
      stop('Please specify "group.prob".')
    }
  }
  reciprocal <- list("group.prob" = group.prob,
                "recip.prob" = recip.prob, 
                "selfloop.recip" = selfloop.recip)
  structure(list("reciprocal" = reciprocal),
            class = "panet.control")
}
