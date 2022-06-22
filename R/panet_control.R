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

#' @importFrom utils modifyList
NULL

#' Add components to the control list
#'
#' `+` is used to combine components to control the PA network generation
#' process. Available components are \code{rpactl.scenario()},
#' \code{rpactl.edgeweight()}, \code{rpactl.newedge()},
#' \code{rpactl.preference()} and \code{rpactl.reciprocal()}.
#'
#' @param e1 A list of class \code{rpactl}.
#' @param e2 A list of class \code{rpactl}.
#'
#' @return A list of class \code{rpactl} with components from
#'   \code{e1} and \code{e2}.
#' @export
#'
#' @examples
#' control <- rpactl.scenario(alpha = 0.5, beta = 0.5) +
#'     rpactl.preference(sparams = c(1, 1, 0, 0, 1),
#'         tparams = c(0, 0, 1, 1, 1))
#'
#' control <- rpactl.scenario(alpha = 1) +
#'     rpactl.edgeweight(distribution = rgamma,
#'         dparams = list(shape = 5, scale = 0.2),
#'         shift = 1)
"+.rpactl" <- function(e1, e2) {
  e1 <- structure(utils::modifyList(e1, e2, keep.null = TRUE), class = "rpactl")
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
#' @param xi Probability of adding an edge between two new nodes.
#' @param rho Probability of adding a new node with a self-loop.
#' @param beta.loop Logical, whether self-loops are allowed under beta scenario.
#'   Default value is \code{TRUE}.
#' @param source.first Logical. Defined for \code{beta} scenario edges of
#'   directed networks. If \code{TRUE}, the source node of a new edge is sampled
#'   from existing nodes before the target node is sampled; if \code{FALSE}, the
#'   target node is sampled from existing nodes before the source node is
#'   sampled. Default value is \code{TRUE}.
#' 
#' @return A list of class \code{rpactl} with components \code{alpha},
#'   \code{beta}, \code{gamma}, \code{xi}, \code{rho}, \code{beta.loop} and
#'   \code{source.first} with meanings as explained under 'Arguments'.
#'
#' @export
#'
#' @examples
#' control <- rpactl.scenario(alpha = 0.5, beta = 0.5, beta.loop = FALSE)
#' 
rpactl.scenario <- function(alpha = 1, beta = 0, gamma = 0, xi = 0, rho = 0,
                             beta.loop = TRUE, source.first = TRUE) {
  stopifnot('"alpha + beta + gamma + xi + rho" must equal to 1.' =
            round(alpha + beta + gamma + xi + rho, 10) == 1)
  scenario <- list("alpha" = alpha, 
                   "beta" = beta,
                   "gamma" = gamma, 
                   "xi" = xi,
                   "rho" = rho,
                   "beta.loop" = beta.loop, 
                   "source.first" = source.first)
  structure(list("scenario" = scenario),
            class = "rpactl")
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
#' @return A list of class \code{rpactl} with components 
#'   \code{distribution}, \code{dparams}, and \code{shift} with meanings as 
#'   explained under 'Arguments'.
#' 
#' @export
#'
#' @examples
#' # Edge weight follows Gamma(5, 0.2).
#' control <- rpactl.edgeweight(distribution = rgamma,
#'     dparams = list(shape = 5, scale = 0.2),
#'     shift = 0)
#'
#' # Constant edge weight
#' control <- rpactl.edgeweight(shift = 2)
#' 
rpactl.edgeweight <- function(distribution = NA,
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
            class = "rpactl")
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
#'   directed networks, when \code{node.replace} is \code{FALSE}, sampled source
#'   and target nodes in the same step are all different from each other,
#'   \code{beta.loop}, \code{snode.replace} and \code{tnode.replace} will be set
#'   as \code{FALSE}.
#' 
#' @return A list of class \code{rpactl} with components \code{distribution},
#'   \code{dparams}, \code{shift}, \code{snode.replace}, \code{tnode.replace} and 
#'   \code{node.replace} with meanings as explained under 'Arguments'.
#'
#' @export
#'
#' @examples
#' control <- rpactl.newedge(distribution = rpois,
#'     dparams = list(lambda = 2),
#'     shift = 1,
#'     node.replace = FALSE)
rpactl.newedge <- function(distribution = NA,
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
  structure(list("newedge" = newedge), class = "rpactl")
}

#' Set parameters for source and target preference function
#'
#' @param sparams Parameters of the source preference function for directed
#'   networks. Probability of choosing an existing node as the source node is
#'   proportional to \code{sparams[1] * out-strength^sparams[2] + sparams[3] *
#'   in-strength^sparams[4] + sparams[5]}.
#' @param tparams Parameters of the target preference function for directed
#'   networks. Probability of choosing an existing node as the source node is
#'   proportional to \code{tparams[1] * out-strength^tparams[2] + tparams[3] *
#'   in-strength^tparams[4] + tparams[5]}.
#' @param params Parameters of the preference function for undirected networks.
#'   Probability of choosing an existing node is proportional to
#'   \code{strength^param[1] + param[2]}.
#' 
#' @return A list of class \code{rpactl} with components \code{sparams},
#'   \code{tparams}, and \code{params} with meanings as explained under 
#'   'Arguments'.
#'
#' @export
#'
#' @examples
#' control <- rpactl.preference(sparams = c(1, 2, 0, 0, 0.1),
#'     tparams = c(0, 0, 1, 2, 0.1))
rpactl.preference <- function(sparams = c(1, 1, 0, 0, 1),
                               tparams = c(0, 0, 1, 1, 1),
                               params = c(1, 1)) {
  preference <- list("sparams" = sparams,
               "tparams" = tparams,
               "params" = params)
  stopifnot(sparams[5] >= 0 & tparams[5] >= 0 & params[2] >= 0)
  structure(list("preference" = preference), 
            class = "rpactl")
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
#' @return A list of class \code{rpactl} with components 
#'   \code{group.prob}, \code{recip.prob}, and \code{selfloop.recip} with
#'   meanings as explained under 'Arguments'.
#'
#' @export
#'
#' @examples
#' control <- rpactl.reciprocal(group.prob = c(0.4, 0.6),
#'     recip.prob = matrix(runif(4), ncol = 2))
rpactl.reciprocal <- function(group.prob = NULL,
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
            class = "rpactl")
}
