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

#' Checks if the input is a \code{rpacontrol} object
#'
#' @param control An \code{rpacontrol} object.
#' @return Logical, \code{TRUE} if the input is a \code{rpacontrol} object.
#' @keywords internal
#' 
is_rpacontrol <- function(control) {
    return(inherits(control, "rpacontrol"))
}

#' Prints preference functions in terminal
#' 
#' @param control An object of class \code{rpacontrol}.
#' @return Returns \code{NULL} invisibly.
#' @keywords internal
#' 
print_control_preference <- function(control, directed = NULL) {
    cat("Preference functions:\n")
    spref <- " - spref: "
    tpref <- " - tpref: "
    pref <- " - pref: "
    if (control$preference$ftype == "default") {
        tmp <- control$preference$sparams
        spref <- paste0(
            spref,
            tmp[1], " * outs^", tmp[2], " + ", tmp[3], " * ins^", tmp[4],
            " + ", tmp[5],
            sep = ""
        )
        tmp <- control$preference$tparams
        tpref <- paste0(
            tpref,
            tmp[1], " * outs^", tmp[2], " + ", tmp[3], " * ins^", tmp[4],
            " + ", tmp[5],
            sep = ""
        )
        tmp <- control$preference$params
        pref <- paste0(
            pref,
            "s^", tmp[1], " + ", tmp[2],
            sep = ""
        )
    } else if (control$preference$ftype == "customized") {
        spref <- paste0(
            spref,
            control$preference$spref
        )
        tpref <- paste0(
            tpref,
            control$preference$tpref
        )
        pref <- paste0(
            pref,
            control$preference$pref
        )
    } else {
        stop("Preference function type is not valid.")
    }
    if (is.null(directed)) {
        cat(spref, tpref, pref, sep = "\n")
    } else if (directed) {
        cat(spref, tpref, sep = "\n")
    } else {
        cat(pref, "\n")
    }
    cat("\n")

    invisible(NULL)
}

#' Prints \code{rpacontrol} in terminal
#'
#' @param x An object of class \code{rpacontrol}.
#' @param control_name A string of control name.
#' @param control_description A list of control descriptions.
#' @return Returns \code{NULL} invisibly.
#' @keywords internal
#' 
print_control_details <- function(x, control_name, control_description) {
    if (control_name == "preference") {
        print_control_preference(control = x, directed = NULL)
        return(invisible(NULL))
    }
    cat(control_description[[control_name]], ":\n", sep = "")
    for (name in names(x[[control_name]])) {
        value <- x[[control_name]][[name]]
        if (is.function(value)) {
            value <- as.character(substitute(rgamma))
        } else if (is.list(value) && length(value) > 0) {
            value <- sapply(
                seq_along(value),
                function(i) paste0(names(value)[i], "=", value[i], " ")
            )
        }
        cat(paste0(" - ", name, ": "), value, "\n", sep = "")
    }
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
    tmp <- rpa_control_list()
    control_names <- names(x)
    control_descriptions <- tmp$control_descriptions
    rm(tmp)

    for (each in control_names) {
        print_control_details(x, each, control_descriptions)
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

    tmp <- rpa_control_list()
    control_names <- tmp$control_names
    control_descriptions <- tmp$control_descriptions
    rm(tmp)

    count <- 0
    cat("Specified control(s):\n")
    cat("--------------------\n")
    for (each in control_names) {
        if (!identical(control_default[[each]], object[[each]])) {
            print_control_details(object, each, control_descriptions)
            count <- 1
        }
    }
    if (count == 0) cat("None\n\n")

    count <- 0
    cat("\nDefault (unspecified) controls:\n")
    cat("------------------------------\n")
    for (each in control_names) {
        if (identical(control_default[[each]], object[[each]])) {
            print_control_details(object, each, control_descriptions)
            count <- 1
        }
    }
    if (count == 0) cat("None\n")

    invisible(object)
}
