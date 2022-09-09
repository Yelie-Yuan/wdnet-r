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

#' @importFrom Rcpp sourceCpp
NULL

#' Compile preference functions via \code{Rcpp}.
#'
#' @param preference A list for defining the preference functions.
#' @param fname A list of function names: \code{spref_func} represents the name
#'   of source preference function; \code{spref_XPtr} represents the function
#'   for returning the pointer of \code{spref_func}; \code{tpref_func},
#'   \code{tpref_XPtr}, \code{pref_func} and \code{pref_XPtr} are defined
#'   analogously. \code{pref_func} and \code{pref_XPtr} are defined for
#'   undirected networks.
#'
#' @return Preference functions and external pointers.
compile_pref_func <- function(preference,
                              fname = list(
                                spref_func = "spref_func",
                                tpref_func = "tpref_func",
                                pref_func = "pref_func",
                                spref_XPtr = "spref_XPtr",
                                tpref_XPtr = "tpref_XPtr",
                                pref_XPtr = "pref_XPtr"
                              )) {
  compile_spref <- compile_tpref <- compile_pref <- FALSE
  cpp_code <- paste("// [[Rcpp::depends(RcppArmadillo)]]\n",
    "#include <RcppArmadillo.h>\n",
    "#include <math.h>\n",
    "using namespace Rcpp;\n",
    "typedef double (*funcPtrD)(double x, double y);\n", sep = "")
  get_func_d <- function(fname, fname_XPtr, ret) {
    return(paste("\n", 
                 "// [[Rcpp::export]]\n",
                 "double ", fname, "(double outs, double ins) {\n",
                 "return ", ret, ";\n",
                 "}\n", 
                 "// [[Rcpp::export]]\n", 
                 "XPtr<funcPtrD> ", fname_XPtr, "() {\n",
                 "return(XPtr<funcPtrD>(new funcPtrD(", fname, ")));\n",
                 "}", sep = ""))
  }
  get_func_und <- function(fname, fname_XPtr, ret) {
    return(paste("\n",
                 "typedef double (*funcPtrUnd)(double x);\n", 
                 "// [[Rcpp::export]]\n",
                 "double ", fname, "(double s) {\n",
                 "return ", ret, ";\n",
                 "}\n", 
                 "// [[Rcpp::export]]\n", 
                 "XPtr<funcPtrUnd> ", fname_XPtr, "() {\n",
                 "return(XPtr<funcPtrUnd>(new funcPtrUnd(", fname, ")));\n",
                 "}", sep = ""))
  }
  if (typeof(preference$spref) == "character") {
    compile_spref <- TRUE
    cpp_code <- paste(cpp_code, get_func_d(fname$spref_func,
                                           fname$spref_XPtr,
                                           preference$spref))
  }
  else if (typeof(preference$spref) != "externalptr") {
    stop('Type of "spref" must be "externalptr" or "character".')
  }
  if (typeof(preference$tpref) == "character") {
    compile_tpref <- TRUE
    cpp_code <- paste(cpp_code, get_func_d(fname$tpref_func,
                                           fname$tpref_XPtr,
                                           preference$tpref))
  }
  else if (typeof(preference$tpref) != "externalptr") {
    stop('Type of "tpref" must be "externalptr" or "character".')
  }

  if (typeof(preference$pref) == "character") {
    compile_pref <- TRUE
    cpp_code <- paste(cpp_code, get_func_und(fname$pref_func,
                                             fname$pref_XPtr,
                                             preference$pref))
  }
  else if (typeof(preference$pref) != "externalptr") {
    stop('Type of "pref" must be "externalptr" or "character".')
  }
  if (compile_spref | compile_tpref | compile_pref) {
    cat("Compiling preference function(s)...\t")
    Rcpp::sourceCpp(code = cpp_code)
    cat("Done.\n")
  }
  preference$spref.pointer <- ifelse(compile_spref,
                                     get(fname$spref_XPtr, envir = .GlobalEnv)(),
                                     preference$spref)
  preference$tpref.pointer <- ifelse(compile_tpref,
                                     get(fname$tpref_XPtr)(),
                                     preference$tpref)
  preference$pref.pointer <- ifelse(compile_pref,
                                    get(fname$pref_XPtr)(),
                                    preference$pref)
  # test functions
  test_pref_func_directed(preference$spref.pointer, 1, 1)
  test_pref_func_directed(preference$tpref.pointer, 1, 1)
  test_pref_func_undirected(preference$pref.pointer, 1)
  preference
}