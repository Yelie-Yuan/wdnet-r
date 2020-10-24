##
## wdnet: Weighted directed network
## Copyright (C) 2020  Yelie Yuan, Panpan Zhang, and Jun Yan
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

##' wenet: Weighted and Directed Networks
##'
##' This package provides functions to conduct network analysis
##' \itemize{
##' \item B-splines
##' \item M-splines
##' \item I-splines
##' \item convex splines (C-splines)
##' \item generalized Bernstein polynomials
##' \item their integrals (except C-splines) and derivatives of given order by
##' close-form recursive formulas
##' }
##'
##' In addition to the R interface, it also provides a C++ header-only library
##' integrated with \pkg{Rcpp}, which allows construction of spline basis matrix
##' directly in C++ with the help of \pkg{Rcpp} and \pkg{RcppArmadillo}.  So it
##' can also be treated as one of the \pkg{Rcpp*} packages.  A toy example
##' package that uses the C++ interface is available at
##' <https://github.com/wenjie2wang/example-pkg-Rcpp-splines2>.
##'
##' @docType package
##' @name wdnet
##' @useDynLib wdnet, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @importFrom Rcpp evalCpp

NULL
