##
## wdnet: Weighted directed network
## Copyright (C) 2024  Yelie Yuan, Tiandong Wang, Jun Yan and Panpan Zhang
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

#' @importFrom igraph indent_print
#' @importFrom utils capture.output
NULL


#' Checks whether the input is a \code{rpacontrol} object
#'
#' @param control A \code{rpacontrol} object.
#' @return Logical, \code{TRUE} if the input is a \code{rpacontrol} object.
#' @keywords internal
#'
is_rpacontrol <- function(control) {
    return(inherits(control, "rpacontrol"))
}


#' Prints \code{rpa_control_scenario()} in terminal
#'
#' @param control A list of control parameters for
#'   \code{rpa_control_scenario()}.
#'
#' @return Returns \code{NULL} invisibly.
#' @keywords internal
#'
print_control_scenario <- function(control) {
    cat("Edge scenarios:\n")
    cat(" - alpha: ", control$alpha, "\n", sep = "")
    cat(" - beta: ", control$beta, "\n", sep = "")
    cat(" - gamma: ", control$gamma, "\n", sep = "")
    cat(" - xi: ", control$xi, "\n", sep = "")
    cat(" - rho: ", control$rho, "\n", sep = "")
    cat(" - beta.loop: ", control$beta.loop, "\n", sep = "")
    cat(" - source.first: ", control$source.first, "\n", sep = "")
    invisible(NULL)
}


#' Prints \code{rpa_control_edgeweight()} in terminal
#'
#' @param control A list of control parameters for
#'   \code{rpa_control_edgeweight()}.
#'
#' @return Returns \code{NULL} invisibly.
#' @keywords internal
#'
print_control_edgeweight <- function(control) {
    cat("Edge weights:\n")

    cat(" - sampler: ")
    if (is.null(control$sampler)) {
        cat("NULL; all new edges have weight 1\n")
    } else {
        cat("\n")
        igraph::indent_print(control$sampler, .indent = "   ")
    }

    invisible(NULL)
}


#' Prints \code{rpa_control_newedge()} in terminal
#'
#' @param control A list of control parameters for
#'   \code{rpa_control_newedge()}.
#'
#' @return Returns \code{NULL} invisibly.
#' @keywords internal
#'
print_control_newedge <- function(control) {
    cat("New edges in each step:\n")

    cat(" - sampler: ")
    if (is.null(control$sampler)) {
        cat("NULL; add one new edge at each step\n")
    } else {
        cat("\n")
        igraph::indent_print(control$sampler, .indent = "   ")
    }

    cat(" - snode.replace: ", control$snode.replace, "\n", sep = "")
    cat(" - tnode.replace: ", control$tnode.replace, "\n", sep = "")
    cat(" - node.replace: ", control$node.replace, "\n", sep = "")
    invisible(NULL)
}


#' Prints \code{rpa_control_reciprocal()} in terminal
#'
#' @param control A list of control parameters for
#'   \code{rpa_control_reciprocal()}.
#' @return Returns \code{NULL} invisibly.
#' @keywords internal
#'
print_control_reciprocal <- function(control) {
    cat("Reciprocal edges:\n")

    cat(" - group.prob: ")
    if (is.null(control$group.prob)) {
        cat("NULL\n")
    } else {
        cat(control$group.prob, "\n", sep = " ")
    }

    cat(" - recip.prob: ")
    if (is.null(control$recip.prob)) {
        cat("NULL; no immediate reciprocal edges\n")
    } else {
        cat("\n")
        igraph::indent_print(control$recip.prob, .indent = "   ")
    }
    invisible(NULL)
}


#' Prints \code{rpa_control_preference()} in terminal
#'
#' @param control A list of control parameters for
#'   \code{rpa_control_preference()}.
#' @param directed Logical, whether to print preference functions for directed
#'   networks only. If missing, print preference functions for both directed
#'   and undirected networks.
#'
#' @return Returns \code{NULL} invisibly.
#' @keywords internal
#'
print_control_preference <- function(control, directed) {
    cat("Preference functions:\n")
    cat(" - ftype: ", control$ftype, "\n", sep = "")
    if (control$ftype == "default") {
        spref <- paste0(
            " - sparams: ",
            paste(control$sparams, collapse = " ")
        )
        tpref <- paste0(
            " - tparams: ",
            paste(control$tparams, collapse = " ")
        )
        pref <- paste0(
            " - params: ",
            paste(control$params, collapse = " ")
        )
    } else if (control$ftype == "customized") {
        my_print <- function(pref, type) {
            if (inherits(pref, "XPtr")) {
                tmp <- utils::capture.output(pref)
                if (grepl("pointer:\\ \\(nil\\)", tmp)) {
                    pref <- paste0(" - ", type,
                        ": XPtr; not valid, please recompile the C++ code.",
                        sep = ""
                    )
                } else {
                    pref <- paste0(" - ", type, ": XPtr;", tmp, sep = "")
                }
            } else {
                pref <- paste0(" - ", type, ": ", pref)
            }
        }
        spref <- my_print(control$spref, "spref")
        tpref <- my_print(control$tpref, "tpref")
        pref <- my_print(control$pref, "pref")
    } else {
        stop("Preference function type is not valid.")
    }
    if (missing(directed)) {
        cat(spref, tpref, pref, sep = "\n")
    } else if (directed) {
        cat(spref, tpref, sep = "\n")
    } else {
        cat(pref, "\n")
    }

    invisible(NULL)
}


#' Prints \code{rpacontrol} in terminal
#'
#' @param x An object of class \code{rpacontrol}.
#' @param control_name A string, the name of the control component.
#'
#' @return Returns \code{NULL} invisibly.
#' @keywords internal
#'
print_control_details <- function(x, control_name) {
    switch(control_name,
        "scenario" = print_control_scenario(control = x$scenario),
        "edgeweight" = print_control_edgeweight(control = x$edgeweight),
        "newedge" = print_control_newedge(control = x$newedge),
        "reciprocal" = print_control_reciprocal(control = x$reciprocal),
        "preference" = print_control_preference(control = x$preference)
    )
    cat("\n")
    invisible(NULL)
}


#' Prints \code{rpacontrol} objects
#'
#' These functions print \code{rpacontrol} objects in the terminal.
#' \code{print.rpacontrol()} shows only the current controls, whereas
#' \code{summary.rpacontrol()} includes both specified controls and the
#' unspecified controls that use default values.
#'
#' @param x An object of class \code{rpacontrol}.
#' @param object An object of class \code{rpacontrol}.
#' @param ... Additional arguments.
#' @return Returns the controls invisibly.
#' @rdname print.rpacontrol
#' @method print rpacontrol
#' @export
#' @examples
#'
#' control <- rpa_control_scenario()
#' print(control)
#'
print.rpacontrol <- function(x, ...) {
    control_names <- names(x)

    for (each in control_names) {
        print_control_details(x, each)
    }
    invisible(x)
}

#' @rdname print.rpacontrol
#' @method summary rpacontrol
#' @export
#'
summary.rpacontrol <- function(object, ...) {
    control_default <- rpa_control_default()
    object <- control_default + object

    control_names <- names(rpa_control_default())

    count <- 0
    cat("Specified control(s):\n")
    cat("--------------------\n")
    for (each in control_names) {
        if (!identical(control_default[[each]], object[[each]])) {
            print_control_details(object, each)
            count <- 1
        }
    }
    if (count == 0) cat("None\n\n")

    count <- 0
    cat("\nDefault (unspecified) controls:\n")
    cat("------------------------------\n")
    for (each in control_names) {
        if (identical(control_default[[each]], object[[each]])) {
            print_control_details(object, each)
            count <- 1
        }
    }
    if (count == 0) cat("None\n")

    invisible(object)
}
