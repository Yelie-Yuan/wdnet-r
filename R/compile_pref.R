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

#' @importFrom RcppXPtrUtils checkXPtr cppXPtr
NULL

#' Compile preference functions via \code{Rcpp}.
#'
#' @param preference A list for defining the preference functions.
#'
#' @return Preference functions and external pointers.
#' 
#' @keywords internal
#' 
compile_pref_func <- function(preference) {
  if (inherits(preference$spref, "character")) {
    temp <- paste("double spref(double outs, double ins) { return ",
                  preference$spref, ";}", sep = "")
    preference$spref.pointer <- RcppXPtrUtils::cppXPtr(code = temp)
    rm(temp)
  }
  else if (inherits(preference$spref, "XPtr")) {
    RcppXPtrUtils::checkXPtr(ptr = preference$spref,
                             type = "double",
                             args = c("double", "double"))
    preference$spref.pointer <- preference$spref
  }
  else {
    stop('Class of "spref" must be "XPtr" or "character".')
  }
  if (inherits(preference$tpref, "character")) {
    temp <- paste("double tpref(double outs, double ins) { return ",
                  preference$tpref, ";}", sep = "")
    preference$tpref.pointer <- RcppXPtrUtils::cppXPtr(code = temp)
    rm(temp)
  }
  else if (inherits(preference$tpref, "XPtr")) {
    RcppXPtrUtils::checkXPtr(ptr = preference$tpref,
                             type = "double",
                             args = c("double", "double"))
    preference$tpref.pointer <- preference$tpref
  }
  else {
    stop('Class of "tpref" must be "externalptr" or "character".')
  }
  if (inherits(preference$pref, "character")) {
    temp <- paste("double pref(double s) { return ",
                  preference$pref, ";}", sep = "")
    preference$pref.pointer <- RcppXPtrUtils::cppXPtr(code = temp)
    rm(temp)
  }
  else if (inherits(preference$pref, "XPtr")) {
    RcppXPtrUtils::checkXPtr(ptr = preference$pref,
                             type = "double",
                             args = "double")
    preference$pref.pointer <- preference$pref
  }
  else {
    stop('Class of "pref" must be "externalptr" or "character".')
  }
  preference
}