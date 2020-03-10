# Copyright (C) {2020} {PB, LSTH, GB, PL, MAS}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

##############################################################################
## Data

#' @title Create xml lines of data
#'
#' @description
#' Create a line with continuous data
#'
#' @param data a p x ntaxa matrix of data, with species as colnames
#' @param trait the name of the trait
#' 
#' @return A line of text
##
get_all_data_lines <- function(data, trait){
  fun <- function(x){
    get_data_line(x[2:length(x)], x[1], trait)
  }
  data <- unname(rbind(colnames(data), data))
  data_lines <- apply(data, 2, fun)
  return(data_lines)
}

#' @title Create xml line of data
#'
#' @description
#' Create a line with continuous data
#'
#' @param data a vector of continuous data
#' @param species the name of the species
#' @param trait the name of the trait
#' 
#' @return A line of text
##
get_data_line <- function(data, species, trait){
  data_line <- paste0("\t\t<taxon id=\"", species, "\"> <attr name=\"", trait, "\"> ")
  data_line <- paste0(data_line, paste(data, collapse=" "))
  data_line <- paste0(data_line, " </attr> </taxon>")
  return(data_line)
}

#' @title Create xml lines of dated tips
#'
#' @description
#' Create all lines with dated tips
#'
#' @param dates a date frame wntaxa vector of dates, with species as names
#' 
#' @return Several lines of text
##
get_all_dated_tips <- function(dates){
  fun <- function(x){
    get_dated_tip_line(x[2], x[1])
  }
  data <- unname(rbind(names(dates), dates))
  data_lines <- apply(data, 2, fun)
  return(data_lines)
}

#' @title Create xml line of dated taxa
#'
#' @description
#' Create a line with a dated taxa
#'
#' @param species the name of the species
#' @param date the date
#' 
#' @return A line of text
##
get_dated_tip_line <- function(date, species){
  data_line <- paste0("\t\t<taxon id=\"", species, "\"> <date value=\"", date, "\"")
  data_line <- paste0(data_line, " direction=\"forwards\" units=\"years\"/>")
  data_line <- paste0(data_line, " </taxon>")
  return(data_line)
}

##############################################################################
## Parameters

#' @title Create xml line of matrix
#'
#' @description
#' Create a line with matrix
#'
#' @param value a vector of values
#' @param id id of the parameter
#' @param lower a vector with lower values for the parameter
#' @param indent number of identation needed
#' 
#' @return A line of text
##
get_parameter_line <- function(value, id = NULL, lower = NULL, indent = 4, ref = FALSE){
  if (ref) { # ref to a param only
    if (is.null(id)) stop("Please precise the idref to the parameter.")
    return(paste0(paste(rep("\t", indent), collapse = ""), "<parameter idref=\"", id, "\"/>"))
  }
  param_line <- paste0(paste(rep("\t", indent), collapse = ""), "<parameter")
  if (!is.null(id)) param_line <- paste0(param_line, " id=\"", id, "\"")
  param_line <- paste0(param_line, " value=\"", paste(value, collapse=" "), "\"")
  if (!is.null(lower)) {
    if (length(value) != length(lower)) stop("lower should have the same dimension as value.")
    param_line <- paste0(param_line, " lower=\"", paste(lower, collapse=" "), "\"")
  }
  param_line <- paste0(param_line, "/>")
  return(param_line)
}

#' @title Create an xml block with several parameter declarations
#'
#' @description
#' Create an xml block with several parameter declarations
#'
#' @param values a vector of values
#' @param ids id of the parameters
#' @param ref should the parameters be ref to existing parameters ?
#' 
#' @return A bloc of text
##
get_parameter_lines <- function(values, ids = NULL, ref = FALSE){
  if (is.null(ids)){
    return(sapply(values, get_parameter_line))
  } else {
    fun <- function(x){
      get_parameter_line(x[2:length(x)], x[1], ref = ref)
    }
    values <- unname(rbind(ids, values))
    return(apply(values, 2, fun))
  }
}

#' @title Create xml blox of matrix
#'
#' @description
#' Create a blox with matrix
#'
#' @param values a matrix of values
#' @param ids id of the parameters
#' 
#' @return A bloc of text
##
get_matrix_bloc <- function(values, ids = NULL){
  if (is.null(ids)){
    return(apply(values, 2, get_parameter_line))
  } else {
    fun <- function(x){
      get_parameter_line(x[2:length(x)], x[1])
    }
    values <- unname(rbind(ids, values))
    return(apply(values, 2, fun))
  }
}

get_diagonal_matrix <- function(name, value = NULL, lower = NULL) {
  return(c(
    paste0("\t<diagonalMatrix id=\"", name, ".matrix\">"),
    get_parameter_line(value = value, id = paste0(name, ".values"), lower = lower, indent = 2),
    "\t</diagonalMatrix>"
  ))
}

get_compound_diagonal_matrix <- function(name, value = NULL, lower = NULL) {
  return(c(
    paste0("\t<diagonalMatrix id=\"", name, ".matrix\">"),
    paste0("\t\t<compoundParameter id=\"", name, ".values\">"),
    sapply(1:length(value), function(z) get_parameter_line(value = value[z], id = paste0(name, ".values.", z), lower = lower[z], indent = 3)),
    "\t\t</compoundParameter>",
    "\t</diagonalMatrix>"
  ))
}

get_compound_symmetric_matrix <- function(name, valueDiag = NULL, valueOffDiag = NULL, lowerDiag = NULL, indent = 1) {
  indent_text <- paste0(rep("\t", indent), collapse = "")
  return(c(
    paste0(indent_text, "<compoundSymmetricMatrix id=\"", name, ".matrix\" asCorrelation=\"true\" isCholesky=\"true\">"),
    paste0(indent_text, "\t<diagonal>"),
    get_parameter_line(value = valueDiag, id = paste0(name, ".diagonal"), lower = lowerDiag, indent = indent + 2),
    paste0(indent_text, "\t</diagonal>"),
    paste0(indent_text, "\t<offDiagonal>"),
    get_parameter_line(value = valueOffDiag, id = paste0(name, ".offDiagonal"), indent = indent + 2),
    paste0(indent_text, "\t</offDiagonal>"),
    paste0(indent_text, "</compoundSymmetricMatrix>")
  ))
}

#' @title Create xml bloc of rate
#'
#' @description
#' Create a block
#'
#' @param value a vector of values
#' @param id id of the parameter
#' 
#' @return A line of text
##
get_rate_bloc <- function(value, id = NULL, ref = FALSE){
  param_line <- get_parameter_line(value, id, ref = ref)
  header <- c("\t\t<strictClockBranchRates>", "\t\t\t<rate>")
  footer <- c("\t\t\t</rate>", "\t\t</strictClockBranchRates>")
  return(c(header, param_line, footer))
}

#' @title Create xml blox of matrix
#'
#' @description
#' Create a blox with matrix
#'
#' @param values a vector of values
#' @param ids id of the parameters
#' 
#' @return A bloc of text
##
get_rate_blocs <- function(values, ids = NULL, ref = FALSE){
  if (is.null(ids)){
    return(sapply(values, get_rate_bloc))
  } else {
    fun <- function(x){
      get_rate_bloc(x[2:length(x)], x[1], ref)
    }
    values <- unname(rbind(ids, values))
    return(apply(values, 2, fun))
  }
}

#' @title Create xml bloc of random walk
#'
#' @description
#' Create a block
#'
#' @param param_id id of the parameter
#' @param windowSize window size of the random walk
#' @param weight weight of the operator
#' 
#' @return A line of text
##
get_operator_random_walk_bloc <- function(param_id, windowSize = 0.1, weight = 1, logBoundary = FALSE, LKJ = FALSE, dimLKJ = NULL){
  header <- c("\t\t<randomWalkOperator")
  boundaryCondition <- ""
  if (logBoundary) boundaryCondition <- "boundaryCondition=\"log\""
  header <- paste0(header, " windowSize=\"", windowSize, "\" weight=\"", weight, "\" ", boundaryCondition, ">")
  param <- paste0("\t\t\t<parameter idref=\"", param_id, "\"/>")
  footer <- c("\t\t</randomWalkOperator>")
  if (LKJ) {
    if (is.null(dimLKJ)) stop("The dimension of the LKJ transform must be specified.")
    param <- c(param, paste0("\t\t\t<LKJTransform dimension=\"", dimLKJ, "\"/>"))
  }
  return(c(header, param, footer))
}

get_operator_random_walk_blocs <- function(param_ids, windowSize = 0.1, weight = 1, logBoundary = FALSE){
  return(as.vector(sapply(param_ids, get_operator_random_walk_bloc, windowSize = windowSize, weight = weight, logBoundary = logBoundary)))
}

get_operator_random_walk_spherical_bloc <- function(param_id, dim, windowSize = 0.1, weight = 1){
  header <- c("\t\t<randomWalkOperator")
  header <- paste0(header, " windowSize=\"", windowSize, "\" weight=\"", weight, "\">")
  param <- paste0("\t\t\t<parameter idref=\"", param_id, "\"/>")
  trans <- paste0("\t\t\t<sphericalTransform dim=\"", dim, "\"/>")
  footer <- c("\t\t</randomWalkOperator>")
  return(c(header, param, trans, footer))
}

get_operator_random_walk_spherical_blocs <- function(param_ids, dim, windowSize = 0.1, weight = 1){
  # return(as.vector(sapply(param_ids, get_operator_random_walk_spherical_bloc, dim = dim, windowSize = windowSize, weight = weight)))
  header <- c("\t\t<randomWalkOperator")
  header <- paste0(header, " windowSize=\"", windowSize, "\" weight=\"", weight, "\">")
  param <- paste0("\t\t\t<matrixParameter idref=\"", param_ids, "\"/>")
  trans <- paste0("\t\t\t<sphericalTransform dim=\"", dim, "\"/>")
  footer <- c("\t\t</randomWalkOperator>")
  return(c(header, param, trans, footer))
}

#' @title Create xml bloc of normal prior
#'
#' @description
#' Create a block
#'
#' @param param_id id of the parameter
#' @param mean 
#' @param stdev
#' 
#' @return A line of text
##
get_prior_normal_bloc <- function(param_id, mean = 0, stdev = 5, half = FALSE){
  normal <- "normalPrior"
  if (half) normal <- "halfNormalPrior"
  header <- paste0("\t\t\t\t<", normal)
  header <- paste0(header, " id=\"", param_id, ".prior\"")
  header <- paste0(header, " mean=\"", mean, "\" stdev=\"", stdev, "\">")
  param <- paste0("\t\t\t\t\t<parameter idref=\"", param_id, "\"/>")
  footer <- paste0("\t\t\t\t</", normal, ">")
  return(c(header, param, footer))
}

get_prior_normal_blocs <- function(param_ids,  mean = 0, stdev = 5, half = FALSE){
  return(as.vector(sapply(param_ids, get_prior_normal_bloc, mean = mean, stdev = stdev, half = half)))
}

#' @title Create xml bloc of log
#'
#' @description
#' Create a block
#'
#' @param param_id id of the parameter
#' 
#' @return A line of text
##
get_log_bloc <- function(param_id, mean = 0, stdev = 5){
  return(paste0("\t\t\t\t<parameter idref=\"", param_id, "\"/>"))
}

get_log_blocs <- function(param_ids,  mean = 0, stdev = 5){
  return(as.vector(sapply(param_ids, get_log_bloc)))
}

shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

circulating <- function(x) {
  C <- matrix(NA, length(x), length(x))
  C[, 1] <- x
  for (i in 2:length(x)) {
    C[, i] <- shifter(x, i-1)
  }
  return(C)
}

#' @title Create xml bloc of attenuation matrix
#'
#' @description
#' Create a block
#'
#' @param method one of "full", "symmetric" of "diagonal"
#' 
#' @return A line of text
##
get_attenuation_matrix <- function(method) {
  if (method == "full") {
    ## Init
    header <- c("\t<compoundEigenMatrix id=\"attenuation.matrix\">",
                "\t\t<eigenValues>")
    line_values <- paste0("\t\t\t\t<parameter id=\"attenuation.values\" value=\"", paste0(1:p_base/10, collapse = " "), "\" />")
    middle <- c("\t\t</eigenValues>",
                "\t\t<eigenVectors>",
                "\t\t\t<matrixParameter id=\"attenuation.vectors\">")
    attVecs <- circulating(1:p_base / 20)
    ids <- paste0("attVec.col", 1:p_base)
    lines_vectors <- get_matrix_bloc(attVecs[-nrow(attVecs), ], ids)
    footer <- c("\t\t\t</matrixParameter>",
                "\t\t</eigenVectors>",
                "\t</compoundEigenMatrix>")
    return(c(header, line_values, middle, lines_vectors, footer))
  }
  if (method == "diagonal") {
    att_mat <- c("\t\t<diagonalMatrix id=\"attenuation.matrix\">",
              get_parameter_line(rep(1, p_base), "attenuation.values", lower = rep(0, p_base)),
              "\t\t</diagonalMatrix>")
    attVecs <- diag(rep(1, p_base))
    ids <- paste0("attVec.col", 1:p_base)
    lines_vectors <- get_matrix_bloc(attVecs[-nrow(attVecs), ], ids)
    att_vecs <- c("\t\t<matrixParameter id=\"attenuation.vectors\">",
                  lines_vectors,
             "\t\t</matrixParameter>")
    return(c(att_mat, att_vecs))
  }
  if (method == "symmetric") {
    ## Init
    header <- c("\t<compoundSymmetricMatrix id=\"attenuation.matrix\" asCorrelation=\"true\" isCholesky=\"true\">",
                "\t\t<diagonal>")
    line_values <- get_parameter_line(rep(1, p_base), "attenuation.values", lower = rep(0, p_base))
    middle <- c("\t\t</diagonal>",
                "\t\t<offDiagonal>")
    attVecs <- circulating(1:p_base / 20)
    ids <- paste0("attVec.col", 1:p_base)
    lines_vectors <- paste0("\t\t\t<parameter id=\"attenuation.vectors\" value=\"", paste0(rep(0, p_base * (p_base - 1) / 2), collapse = " "), "\" />")
    footer <- c("\t\t</offDiagonal>",
                "\t</compoundSymmetricMatrix>")
    return(c(header, line_values, middle, lines_vectors, footer))
  }
}

#' @title Create xml bloc of operators for attenuation matrix
#'
#' @description
#' Create a block
#'
#' @param method one of "full", "symmetric" of "diagonal"
#' 
#' @return A line of text
##
get_operators_attenuation <- function(method, windowSize = 0.1, weightVal = 1, weightVec = 1) {
  if (method == "full") {
    op_vecs <- get_operator_random_walk_spherical_blocs("attenuation.vectors", p_base-1, windowSize = windowSize, weight = weightVec)
    op_vals <- c(paste0("\t\t<randomWalkOperator windowSize=\"", windowSize, "\" weight=\"", weightVal, "\">"),
                 "\t\t\t<parameter idref=\"attenuation.values\"/>",
                 paste0("\t\t\t<positiveOrderedTransform dim=\"", p_base, "\"/>"),
                 # "\t\t\t<transform type=\"positiveOrdered\"/>",
                 "\t\t</randomWalkOperator>")
    return(c(op_vecs, op_vals))
  }
  if (method == "diagonal") {
    op_vals <- c(paste0("\t\t<randomWalkOperator windowSize=\"", windowSize, "\" weight=\"", weightVal, "\">"),
                 "\t\t\t<parameter idref=\"attenuation.values\"/>",
                 "\t\t</randomWalkOperator>")
    return(op_vals)
  }
  if (method == "symmetric") {
    op_vecs <- c(paste0("\t\t<randomWalkOperator windowSize=\"", windowSize, "\" weight=\"", weightVec, "\">"),
                 "\t\t\t<parameter idref=\"attenuation.vectors\"/>",
                 paste0("\t\t\t<LKJTransform dimension=\"", p_base, "\"/>"),
                 "\t\t</randomWalkOperator>")
    op_vals <- c(paste0("\t\t<randomWalkOperator windowSize=\"", windowSize, "\" weight=\"", weightVal, "\">"),
                 "\t\t\t<parameter idref=\"attenuation.values\"/>",
                 "\t\t</randomWalkOperator>")
    return(c(op_vecs, op_vals))
  }
}

#' @title Create xml bloc of priors for attenuation matrix
#'
#' @description
#' Create a block
#'
#' @param method one of "full", "symmetric" of "diagonal"
#' 
#' @return A line of text
##
get_prior_attenuation <- function(method) {
  if (method == "full") {
    ## Prior att vals
    header <- c("\t<multivariateNormalPrior id=\"prior.attVals\">",
                "\t\t<meanParameter>")
    line_values <- paste0("\t\t\t\t<parameter id=\"prior.attVals.mean\" value=\"", paste0(rep(0.0, p_base), collapse = " "), "\" />")
    middle <- c("\t\t</meanParameter>",
                "\t\t<precisionParameter>",
                "\t\t\t<matrixParameter id=\"prior.attVals.precision\">")
    prec <- diag(rep(1, p_base))
    lines_vectors <- get_matrix_bloc(prec)
    footer <- c("\t\t\t</matrixParameter>",
                "\t\t</precisionParameter>",
                "\t\t<data>",
                "\t\t\t<parameter idref=\"attenuation.values\"/>",
                "\t\t</data>",
                paste0("\t\t<positiveOrderedTransform dim=\"", p_base, "\"/>"),
                # "\t\t<transform type=\"positiveOrdered\"/>",
                "\t</multivariateNormalPrior>")
    prior_val <- c(header, line_values, middle, lines_vectors, footer)
    ## Prior att vecs           
    prior_vec <- c("\t<sphericalBetaPrior id=\"prior.attVecs\" shapeParameter=\"1.0\">",
                   "\t\t<data>",
                   "\t\t\t<matrixParameter idref=\"attenuation.vectors\"/>",
                   "\t\t</data>",    
                   "\t</sphericalBetaPrior>")
    ## Res
    return(c(prior_val, "", prior_vec))
  }
  if (method == "diagonal") {
    return(c("\t\t\t\t<gammaPrior id=\"prior.attVals\" shape=\"1.0\" scale=\"3.0\" offset=\"0.0\">",
             "\t\t\t\t\t<parameter idref=\"attenuation.values\"/>",
             "\t\t\t\t</gammaPrior>"))
  }
  if (method == "symmetric") {
    prior_val <- c("\t\t\t\t<gammaPrior id=\"prior.attVals\" shape=\"1.0\" scale=\"3.0\" offset=\"0.0\">",
                   "\t\t\t\t\t<parameter idref=\"attenuation.values\"/>",
                   "\t\t\t\t</gammaPrior>")
    prior_vec <- c(paste0("\t\t\t\t<LKJCorrelationPrior id=\"prior.attVecs\" shapeParameter=\"1.0\" dimension=\"", p_base, "\">"),
                   "\t\t\t\t\t<data>",
                   "\t\t\t\t\t\t<parameter idref=\"attenuation.vectors\"/>",
                   "\t\t\t\t\t</data>",
                   "\t\t\t\t</LKJCorrelationPrior>")
    return(c(prior_val, prior_vec))
  }
}