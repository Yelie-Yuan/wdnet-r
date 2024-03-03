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

#' @importFrom RcppXPtrUtils checkXPtr cppXPtr
NULL

#' Compile preference functions via \code{RcppXPtrUtils}.
#'
#' @param preference A list for defining the preference functions.
#' @param directed Logical, whether to compile the preference functions for
#'   directed networks. If missing, the preference functions
#'   for both directed and undirected networks will be compiled.
#'
#' @return Returns the input list and their corresponding external pointers.
#'
#' @keywords internal
#'
compile_pref_func <- function(preference, directed) {
  if (missing(directed) || directed) {
    if (inherits(preference$spref, "character")) {
      tmp <- paste("double spref(double outs, double ins) { return ",
        preference$spref, ";}",
        sep = ""
      )
      preference$spref.pointer <- RcppXPtrUtils::cppXPtr(code = tmp)
      rm(tmp)
    } else if (inherits(preference$spref, "XPtr")) {
      tryCatch(
        {
          RcppXPtrUtils::checkXPtr(
            ptr = preference$spref,
            type = "double",
            args = c("double", "double")
          )
        },
        error = function() {
          stop('Incorrect argument or return type for "spref"; all should be "double".')
        }
      )
      preference$spref.pointer <- preference$spref
      tmp <- utils::capture.output(preference$spref.pointer)
      if (grepl("pointer:\\ \\(nil\\)", tmp)) {
        stop('"XPtr" for "spref" is not valid, please recompile the C++ code.')
      }
      rm(tmp)
    } else {
      stop('The class of "spref" must be either "XPtr" or "character".')
    }
    if (inherits(preference$tpref, "character")) {
      tmp <- paste("double tpref(double outs, double ins) { return ",
        preference$tpref, ";}",
        sep = ""
      )
      preference$tpref.pointer <- RcppXPtrUtils::cppXPtr(code = tmp)
      rm(tmp)
    } else if (inherits(preference$tpref, "XPtr")) {
      tryCatch(
        {
          RcppXPtrUtils::checkXPtr(
            ptr = preference$tpref,
            type = "double",
            args = c("double", "double")
          )
        },
        error = function() {
          stop('Incorrect argument or return type for "tpref"; all should be "double".')
        }
      )
      preference$tpref.pointer <- preference$tpref
      tmp <- utils::capture.output(preference$tpref.pointer)
      if (grepl("pointer:\\ \\(nil\\)", tmp)) {
        stop('"XPtr" for "tpref" is not valid, please recompile the C++ code.')
      }
      rm(tmp)
    } else {
      stop('The class of "tpref" must be either "XPtr" or "character".')
    }
  }
  if (missing(directed) || !directed) {
    if (inherits(preference$pref, "character")) {
      tmp <- paste("double pref(double s) { return ",
        preference$pref, ";}",
        sep = ""
      )
      preference$pref.pointer <- RcppXPtrUtils::cppXPtr(code = tmp)
      rm(tmp)
    } else if (inherits(preference$pref, "XPtr")) {
      tryCatch(
        {
          RcppXPtrUtils::checkXPtr(
            ptr = preference$pref,
            type = "double",
            args = "double"
          )
        },
        error = function() {
          stop('Incorrect argument or return type for "pref"; all should be "double".')
        }
      )
      preference$pref.pointer <- preference$pref
      tmp <- utils::capture.output(preference$pref.pointer)
      if (grepl("pointer:\\ \\(nil\\)", tmp)) {
        stop('"XPtr" for "pref" is not valid, please recompile the C++ code.')
      }
      rm(tmp)
    } else {
      stop('The class of "pref" must be either "XPtr" or "character".')
    }
  }
  return(preference)
}
