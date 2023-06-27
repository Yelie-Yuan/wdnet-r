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

#' @importFrom utils modifyList
#' @importFrom Rcpp sourceCpp
NULL


#' rpacontrol: Controls the Preferential Attachment (PA) Network Generation
#' Process
#'
#' The \code{rpacontrol} object is designed to control the Preferential
#' Attachment (PA) network generation process within the \code{rpanet()}
#' function. It can have the following components:
#' \itemize{
#' \item \code{scenario}: controls the edge scenarios
#' at each step. For more information, please refer to
#' \code{rpa_control_scenario()}.
#' \item \code{edgeweight}: controls the weight of
#' the edges; see \code{rpa_control_edgeweight()} for details.
#' \item \code{newedge}: controls the creation of
#' new edges at each step; see \code{rpa_control_newedge()}
#' for details.
#' \item \code{preference}: sets preference functions; see
#' \code{rpa_control_preference()} for details.
#' \item \code{reciprocal}: controls the creation of reciprocal
#' edges; see \code{rpa_control_reciprocal()} for details.
#' }
#' @name rpacontrol
#'
NULL


#' Add components to the control list
#'
#' `+` is used to combine components to control the PA network generation
#' process. Available components are \code{rpa_control_scenario()},
#' \code{rpa_control_edgeweight()}, \code{rpa_control_newedge()},
#' \code{rpa_control_preference()} and \code{rpa_control_reciprocal()}.
#'
#' @param e1 A list of class \code{rpacontrol}.
#' @param e2 A list of class \code{rpacontrol}.
#'
#' @return A list of class \code{rpacontrol} with components from \code{e1} and
#'   \code{e2}.
#' @export
#'
#' @examples
#' \donttest{
#' control <- rpa_control_scenario(alpha = 0.5, beta = 0.5) +
#'   rpa_control_preference(
#'     ftype = "customized",
#'     spref = "pow(outs, 2) + 1",
#'     tpref = "pow(ins, 2) + 1"
#'   )
#' }
#'
#' control <- rpa_control_scenario(alpha = 1) +
#'   rpa_control_edgeweight(
#'     sampler = function(n) rgamma(n, shape = 5, scale = 0.2)
#'   )
"+.rpacontrol" <- function(e1, e2) {
  e1 <- structure(
    utils::modifyList(e1, e2, keep.null = TRUE),
    class = "rpacontrol"
  )
  if (is.list(e2$edgeweight$dparams)) {
    e1$edgeweight$dparams <- e2$edgeweight$dparams
  }
  if (is.list(e2$newedge$dparams)) {
    e1$newedge$dparams <- e2$newedge$dparams
  }
  e1
}

#' Control edge scenarios. Defined for \code{rpanet}.
#'
#' @param alpha Probability of adding an edge from a new node to an existing
#'   node.
#' @param beta Probability of adding an edge between existing nodes.
#' @param gamma Probability of adding an edge from an existing node to a new
#'   node.
#' @param xi Probability of adding an edge between two new nodes.
#' @param rho Probability of adding a new node with a self-loop.
#' @param beta.loop Logical. Determines whether self-loops are allowed under the
#'   beta scenario. Default value is \code{TRUE}.
#' @param source.first Logical. Defined for \code{beta} scenario edges of
#'   directed networks. If \code{TRUE}, the source node of a new edge is sampled
#'   from existing nodes before the target node is sampled; if \code{FALSE}, the
#'   target node is sampled from existing nodes before the source node is
#'   sampled. Default value is \code{TRUE}.
#'
#' @return A list of class \code{rpacontrol} with components \code{alpha},
#'   \code{beta}, \code{gamma}, \code{xi}, \code{rho}, \code{beta.loop} and
#'   \code{source.first} with meanings as explained under 'Arguments'.
#'
#' @export
#'
#' @examples
#' control <- rpa_control_scenario(alpha = 0.5, beta = 0.5, beta.loop = FALSE)
#'
rpa_control_scenario <- function(
    alpha = 1, beta = 0, gamma = 0, xi = 0, rho = 0,
    beta.loop = TRUE, source.first = TRUE) {
  stopifnot(
    '"alpha + beta + gamma + xi + rho" must be equal to 1.' =
      round(alpha + beta + gamma + xi + rho, 10) == 1
  )
  scenario <- list(
    "alpha" = alpha,
    "beta" = beta,
    "gamma" = gamma,
    "xi" = xi,
    "rho" = rho,
    "beta.loop" = beta.loop,
    "source.first" = source.first
  )
  structure(list("scenario" = scenario),
    class = "rpacontrol"
  )
}

#' Control weight of new edges. Defined for \code{rpanet}.
#'
#' @param sampler A function used for sampling edge weights. If \code{NULL}, all
#'   new edges will default to a weight of 1. If a function is provided, it must
#'   accept a single argument, \code{n}, and return a vector of length \code{n}
#'   that represents the sampled edge weights.
#'
#' @return A list of class \code{rpacontrol} containing the \code{sampler}
#'   function.
#'
#' @export
#'
#' @examples
#' # Sample edge weights from Gamma(5, 0.2).
#' my_gamma <- function(n) rgamma(n, shape = 5, scale = 0.2)
#' control <- rpa_control_edgeweight(
#'   sampler = my_gamma
#' )
#'
rpa_control_edgeweight <- function(
    sampler = NULL) {
  if (!is.null(sampler)) {
    tryCatch(do.call(sampler, list(5)), error = function(e) {
      message('Invalid "sampler" for rpa_control_edgeweight().')
      stop(e)
    })
  }

  edgeweight <- list(
    "sampler" = sampler
  )
  structure(list("edgeweight" = edgeweight),
    class = "rpacontrol"
  )
}

#' Control new edges in each step. Defined for \code{rpanet}.
#'
#' @param sampler A function used for sampling the number of new edges to be
#'   added at each step. If \code{NULL}, one new edge will be added at each
#'   step. If a function is provided, it must accept a single argument,
#'   \code{n}, and return a vector of length \code{n} that represents the
#'   sampled number of new edges.
#' @param snode.replace Logical. Determines whether the source nodes in the same
#'   step should be sampled with replacement. Defined for directed networks.
#' @param tnode.replace Logical. Determines whether the target nodes in the same
#'   step should be sampled with replacement. Defined for directed networks.
#' @param node.replace Logical. Determines whether the nodes in the same step
#'   should be sampled with replacement. Defined for undirected networks. If
#'   FALSE, self-loops will not be allowed under beta scenario.
#'
#' @return A list of class \code{rpacontrol} with components \code{sampler},
#'   \code{snode.replace}, \code{tnode.replace} and \code{node.replace} with
#'   meanings as explained under 'Arguments'.
#'
#' @export
#'
#' @examples
#' my_rpois <- function(n) rpois(n, lambda = 2) + 1
#' control <- rpa_control_newedge(
#'   sampler = my_rpois,
#'   node.replace = FALSE
#' )
rpa_control_newedge <- function(
    sampler = NULL,
    snode.replace = TRUE,
    tnode.replace = TRUE,
    node.replace = TRUE) {
  if (!is.null(sampler)) {
    tryCatch(do.call(sampler, list(5)), error = function(e) {
      message('Invalid "sampler" for rpa_control_newedge().')
      stop(e)
    })
  }

  newedge <- list(
    "sampler" = sampler,
    "snode.replace" = snode.replace,
    "tnode.replace" = tnode.replace,
    "node.replace" = node.replace
  )
  structure(list("newedge" = newedge), class = "rpacontrol")
}

#' Set preference function(s). Defined for \code{rpanet}.
#'
#' @param ftype Preference function type. Either "default" or "customized".
#'   "customized" preference functions require "binary" or "linear" generation
#'   methods. If using default preference functions, \code{sparams},
#'   \code{tparams} and \code{params} must be specified. If using customized
#'   preference functions, \code{spref}, \code{tpref} and \code{pref} must be
#'   specified.
#' @param sparams A numerical vector of length 5 giving the parameters of the
#'   default source preference function. Defined for directed networks.
#'   Probability of choosing an existing node as the source node is proportional
#'   to \code{sparams[1] * out-strength^sparams[2] + sparams[3] *
#'   in-strength^sparams[4] + sparams[5]}.
#' @param tparams A numerical vector of length 5 giving the parameters of the
#'   default target preference function. Defined for directed networks.
#'   Probability of choosing an existing node as the target node is proportional
#'   to \code{tparams[1] * out-strength^tparams[2] + tparams[3] *
#'   in-strength^tparams[4] + tparams[5]}.
#' @param params A numerical vector of length 2 giving the parameters of the
#'   default preference function. Defined for undirected networks. Probability
#'   of choosing an existing node is proportional to \code{strength^params[1] +
#'   params[2].}
#' @param spref Character expression or an object of class \code{XPtr} giving
#'   the customized source preference function. Defined for directed networks.
#'   Default value is \code{"outs + 1"}, i.e., node out-strength + 1. See
#'   Details and Examples for more information.
#' @param tpref Character expression or an object of class \code{XPtr} giving
#'   the customized target preference function. Defined for directed networks.
#'   Default value is \code{"ins + 1"}, i.e., node in-strength + 1.
#' @param pref Character expression or an object of class \code{XPtr} giving the
#'   customized preference function. Defined for undirected networks. Default
#'   value is \code{"s + 1"}, i.e, node strength + 1.
#'
#' @details If choosing customized preference functions, \code{spref},
#'   \code{tpref} and \code{pref} will be used and the network generation method
#'   must be "binary" or "linear". \code{spref} (\code{tpref}) defines the
#'   source (target) preference function, it can be a character expression or an
#'   object of class \code{XPtr}. \itemize{ \item{Character expression: } {it
#'   must be a one-line \code{C++} style expression of \code{outs}
#'   (node out-strength) and
#'   \code{ins} (node in-strength). For example, \code{"pow(outs, 2) + 1"},
#'   \code{"pow(outs, 2) + pow(ins, 2) + 1"}, etc. The expression will be used
#'   to define an \code{XPtr} via \code{RcppXPtrUtils::cppXPtr}. The \code{XPtr}
#'   will be passed to the network generation function. The expression must not
#'   have variables other than \code{outs} and \code{ins}.} \item{\code{XPtr}: }
#'   {an external pointer wrapped in an object of class \code{XPtr} defined via
#'   \code{RcppXPtrUtils::cppXPtr}. An example for defining an \code{XPtr} with
#'   \code{C++} source code is included in Examples. For more information
#'   about passing function pointers, see
#'   \url{https://gallery.rcpp.org/articles/passing-cpp-function-pointers-rcppxptrutils/}.
#'   Please note the supplied \code{C++} function accepts two \code{double}
#'   arguments and returns a \code{double}. The first and second arguments
#'   represent node out- and in-strength, respectively. Note that the \code{XPtr} will
#'   be invalid and cannot be used to control network generation
#'   in another separate R session. Therefore, we recommend preserving the source code of your
#'   preference function for future use.}}
#'
#'   \code{pref} is defined analogously. If using character expression, it must
#'   be a one-line \code{C++} style expression of \code{s} (node strength). If
#'   using \code{XPtr}, the supplied \code{C++} function accepts only one
#'   \code{double} argument and returns a \code{double}.
#'
#' @return A list of class \code{rpacontrol} with components \code{ftype},
#'   \code{sparams}, \code{tparams}, \code{params} or \code{ftype},
#'   \code{spref}, \code{tpref}, \code{pref} with function pointers
#'   \code{spref.pointer}, \code{tpref.pointer}, \code{pref.pointer}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Set source preference as out-strength^2 + in-strength + 1,
#' # target preference as out-strength + in-strength^2 + 1.
#' # 1. use default preference functions
#' ctr1 <- rpa_control_preference(
#'   ftype = "default",
#'   sparams = c(1, 2, 1, 1, 1), tparams = c(1, 1, 1, 2, 1)
#' )
#' # 2. use character expressions
#' ctr2 <- rpa_control_preference(
#'   ftype = "customized",
#'   spref = "pow(outs, 2) + ins + 1", tpref = "outs + pow(ins, 2) + 1"
#' )
#' # 3. define XPtr's with C++ source code
#' spref.pointer <- RcppXPtrUtils::cppXPtr(
#'   code =
#'     "double spref(double outs, double ins) {return pow(outs, 2) + ins + 1;}"
#' )
#' tpref.pointer <- RcppXPtrUtils::cppXPtr(
#'   code =
#'     "double tpref(double outs, double ins) {return outs + pow(ins, 2) + 1;}"
#' )
#' ctr3 <- rpa_control_preference(
#'   ftype = "customized",
#'   spref = spref.pointer,
#'   tpref = tpref.pointer
#' )
#' ret <- rpanet(1e5, control = ctr3)
#' }
rpa_control_preference <- function(
    ftype = c("default", "customized"),
    sparams = c(1, 1, 0, 0, 1),
    tparams = c(0, 0, 1, 1, 1),
    params = c(1, 1),
    spref = "outs + 1",
    tpref = "ins + 1",
    pref = "s + 1") {
  ftype <- match.arg(ftype)
  if (ftype == "default") {
    stopifnot(
      "Length or type of parameter is not valid" =
        all(
          length(sparams) == 5,
          length(tparams) == 5,
          length(params) == 2,
          is.numeric(sparams),
          is.numeric(tparams),
          is.numeric(params)
        )
    )
    preference <- list(
      "ftype" = ftype,
      "sparams" = sparams,
      "tparams" = tparams,
      "params" = params
    )
  } else {
    preference <- list(
      "ftype" = ftype,
      "spref" = spref,
      "tpref" = tpref,
      "pref" = pref
    )
    preference <- compile_pref_func(preference)
  }
  structure(list("preference" = preference),
    class = "rpacontrol"
  )
}

#' Control reciprocal edges. Defined for \code{rpanet}.
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
#' @return A list of class \code{rpacontrol} with components \code{group.prob},
#'   \code{recip.prob}, and \code{selfloop.recip} with meanings as explained
#'   under 'Arguments'.
#'
#' @export
#'
#' @examples
#' control <- rpa_control_reciprocal(
#'   group.prob = c(0.4, 0.6),
#'   recip.prob = matrix(runif(4), ncol = 2)
#' )
rpa_control_reciprocal <- function(
    group.prob = NULL,
    recip.prob = NULL,
    selfloop.recip = FALSE) {
  if (!is.null(group.prob)) {
    stopifnot(
      '"group.prob" must sum to 1.' =
        round(sum(group.prob), 10) == 1
    )
    if (!is.null(recip.prob)) {
      recip.prob <- as.matrix(recip.prob)
      stopifnot(
        '"recip.prob" or "group.prob" is not valid.' =
          length(group.prob) == nrow(recip.prob) &
            nrow(recip.prob) == ncol(recip.prob)
      )
      stopifnot(
        '"recip.prob" is not valid.' =
          all(recip.prob >= 0) &
            all(recip.prob <= 1)
      )
      stopifnot(
        '"group.prob" is not valid.' =
          round(sum(group.prob), 10) == 1 &
            all(group.prob >= 0)
      )
    } else {
      stop('"recip.prob" can not be NULL when "group.prob" is specified.')
    }
  }
  if (is.null(group.prob)) {
    if (!is.null(recip.prob)) {
      stop('"group.prob" is not valid.')
    }
  }
  reciprocal <- list(
    "group.prob" = group.prob,
    "recip.prob" = recip.prob,
    "selfloop.recip" = selfloop.recip
  )
  structure(list("reciprocal" = reciprocal),
    class = "rpacontrol"
  )
}

#' Default controls for \code{rpanet}
#'
#' @return Returns a list of default controls.
#' @keywords internal
#'
rpa_control_default <- function() {
  rpa_control_scenario() +
    rpa_control_edgeweight() +
    rpa_control_newedge() +
    rpa_control_reciprocal() +
    rpa_control_preference()
}
