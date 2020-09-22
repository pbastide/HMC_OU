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

library(tidyverse)

read_log <- function(file, burning){
  ## Read in file
  dat <- read_tsv(file, skip = 3)

  ## Burning
  length_chain <- max(dat$state)
  burning_length <- floor(length_chain * burning)

  ## data
  return(dat %>% filter(state >= burning_length))
}

get_statistics <- function(dat, stats){
  return(dat %>% select(-state) %>% summarise_all(funs_(stats)))
}

get_median <- function(dat){
  return(dat %>% select(-state) %>% summarise_all(list(med = median)))
}

get_mean <- function(dat){
  return(dat %>% select(-state) %>% summarise_all(list(mean = mean)))
}

get_var <- function(dat){
  return(dat %>% select(-state) %>% summarise_all(list(var = var)))
}

get_covar <- function(dat){
  return(var(as.matrix(dat[, -1])))
}

get_HPDI95_up <- function(dat){
  if (!requireNamespace("HDInterval", quietly = TRUE)) {
    stop("Package \"HDInterval\" is needed to compute the HPD interval.",
         call. = FALSE)
  }
  fun <- function(vec){
    res <- HDInterval::hdi(vec)
    return(unname(res[2]))
  }
  return(dat %>% select(-state) %>% summarise_all(funs(Q95 = fun)))
}

get_HPDI95_bottom <- function(dat){
  if (!requireNamespace("HDInterval", quietly = TRUE)) {
    stop("Package \"HDInterval\" is needed to compute the HPD interval.",
         call. = FALSE)
  }
  fun <- function(vec){
    res <- HDInterval::hdi(vec)
    return(unname(res[1]))
  }
  return(dat %>% select(-state) %>% summarise_all(funs(Q05 = fun)))
}

get_mean_min_max <- function(dat){
  fun_min <- function(vec){
    res <- HDInterval::hdi(vec)
    return(unname(res[1]))
  }
  fun_max <- function(vec){
    res <- HDInterval::hdi(vec)
    return(unname(res[2]))
  }
  return(dat %>% select(-state) %>% summarise_all(funs(mean = mean, min = fun_min, max = fun_max)))
}

get_symmetric_matrix <- function(diag, offdiag) {
  mat <- diag(diag)
  n <- length(diag)
  k <- 1
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      mat[i, j] <- offdiag[k]
      mat[j, i] <- offdiag[k]
      k <- k + 1
    }
  }
  return(mat)
}

get_sphere_matrix <- function(vec, dim) {
  res_mat <- matrix(NA, dim, dim)
  for (i in 1:dim) {
    res_mat[, i] <- get_sphere_vector(vec[((dim - 1) * (i - 1) + 1):((dim - 1) * i)])
  }
  return(res_mat)
}

get_sphere_vector <- function(vec) {
  return(unlist(c(vec, sqrt(1 - sum(vec^2)))))
}

normalize <- function(mat) {
  apply(mat, 2, norm_vec)
}

norm_vec <- function(vec) {
  vec <- vec / sum(vec^2)
  if (vec[length(vec)] < 0) vec <- -vec
  return(vec)
}

extract_tip_traits <- function(tree, name) {
  trait_Edges <- sapply(tree$annotations, function(x) unlist(x[[name]]))
  if (is.null(dim(trait_Edges))){
    dim(trait_Edges) <- c(length(trait_Edges), 1)
  } else {
    trait_Edges <- t(trait_Edges)
  }
  trait_Tips <- trait_Edges[match(1:length(tree$tip.label), tree$edge[, 2]), , drop=F]
  rownames(trait_Tips) <- tree$tip.label
  return(trait_Tips)
}


##########################
## Heritability
##########################
##
#' @title Compute Mean Height heritability for BM
#'
#' @description
#' Compute the heritability using the mean height of the tree t_h,
#' using the BM diffusion variance sigma_a and the sampling
#' variance sigma_e, independently for each trait q:
#' h_q = sigma_a * t_h / [ sigma_e + sigma_a * t_h ]
#'
#' @param df a log dataframe containing the "precision.matrix" and the
#' "sampling.precision" for each state of the MCMC.
#' @param meanHeight the mean height of the tree
#' @param dimTrait the dimension of the trait
#'
#' @return A matrix of heritability for each trait and for each state of the MCMC.
##
compute_heritability_bm <- function(df, meanHeight, dimTrait) {
  ## Trait tree variance
  colPrecMat <- grep("precision.matrix", colnames(df))
  if (length(colPrecMat) == 0) colPrecMat <- grep("precisionMatrix", colnames(df))
  traitTreeVar <- apply(df[, colPrecMat], 1, function(z) diag(solve(matrix(z, dimTrait, dimTrait))))
  ## Sampling variance
  samplingVar <- apply(df[, grep("sampling.precision", colnames(df))], 1, function(z) 1 / z)
  ## Both
  bothVars <- rbind(traitTreeVar, samplingVar)
  ## Function
  meanHeight_heritability <- function(both) {
    return(sapply(1:dimTrait,
                  function(i) return(both[i] * meanHeight / (both[i + dimTrait] + both[i] * meanHeight))))
  }
  h_meanHeight <- t(apply(bothVars, 2, meanHeight_heritability))
  return(h_meanHeight)
}

##
#' @title Compute Mean Height heritability for OU
#'
#' @description
#' Compute the heritability using the mean height of the tree t_h,
#' using the OU diffusion variance sigma_a, the sampling
#' variance sigma_e and the selection strenght alpha,
#' independently for each trait q:
#' h_q = sigma_a * t_h' / [ sigma_e + sigma_a * t_h' ]
#' with t_h' the actualized mean height:
#' t_h' = (1 - exp(- 2 * alpha * t_h)) / (2 * alpha)
#'
#' @param df a log dataframe containing the "precisionMatrix", the
#' "sampling.precision" and the "attenuation.values" for each state of the MCMC.
#' @param meanHeight the mean height of the tree
#' @param dimTrait the dimension of the trait
#'
#' @return A matrix of heritability for each trait and for each state of the MCMC.
##
compute_heritability_ou <- function(df, meanHeight, dimTrait) {
  get_ou_variance <- function(t, r, alpha) {
    return(r / (2 * alpha) * (1 - exp(- 2 * alpha * t)))
  }
  ## Trait tree variance
  colPrecMat <- grep("precision.matrix", colnames(df))
  if (length(colPrecMat) == 0) colPrecMat <- grep("precisionMatrix", colnames(df))
  traitTreeVar <- apply(df[, colPrecMat], 1, function(z) diag(solve(matrix(z, dimTrait, dimTrait))))
  ## Sampling variance
  samplingVar <- apply(df[, grep("sampling.precision", colnames(df))], 1, function(z) 1 / z)
  ## Selection Strength
  selectionStrengths <- t(df[, grep("attenuation.values[0-9]", colnames(df))])
  ## All
  allPars <- rbind(traitTreeVar, samplingVar, selectionStrengths)
  ## Function
  meanHeight_heritability <- function(all) {
    return(sapply(1:dimTrait,
                  function(i) return(get_ou_variance(meanHeight, all[i], all[i + 2 * dimTrait]) / (all[i + dimTrait] + get_ou_variance(meanHeight, all[i], all[i + 2 * dimTrait])))))
  }
  return(t(apply(allPars, 2, meanHeight_heritability)))
}

##
#' @title Compute and save variance from precision matrix
#'
#' @param df a log dataframe containing the "precisionMatrix", the
#' "sampling.precision" and the "attenuation.values" for each state of the MCMC.
#' @param dimTrait the dimension of the trait
#'
#' @return Modified df
##
add_variance <- function(df, dimTrait) {
  # precision
  colPrecMat <- grep("precision.matrix", colnames(df))
  if (length(colPrecMat) == 0) colPrecMat <- grep("precisionMatrix", colnames(df))
  traitTreeVar <- apply(df[, colPrecMat], 1, function(z) solve(matrix(z, dimTrait, dimTrait)))
  traitTreeCor <- apply(traitTreeVar, 2, function(z) cov2cor(matrix(z, dimTrait, dimTrait)))
  # sampling
  samplingVar <- apply(df[, grep("sampling.precision", colnames(df))], 1, function(z) 1 / z)
  # name
  traitName <- strsplit(colnames(df[, colPrecMat[1]]), ".", fixed = TRUE)[[1]][1]
  # save
  for (i in 1:dimTrait) {
    for (j in i:dimTrait) {
      if (i == j) df[[paste0(traitName, ".sampling.variance", i - 1, i - 1)]] <- samplingVar[i, ]
      if (i == j) df[[paste0(traitName, ".variance.matrix", i - 1, j - 1)]] <- traitTreeVar[(i - 1) * dimTrait + j, ]
      if (i != j) df[[paste0(traitName, ".correlation.matrix", i - 1, j - 1)]] <- traitTreeCor[(i -1) * dimTrait + j, ]
    }
  }
  return(df)
}


##
#' @title Compute multivariate heritability
#'
#' @description
#' Compute eritability using tip variances.
#'
#' @param varTips a p x p x ntaxa array of variances for each tip
#' @param noise a p x p matrix of variance noise
#'
#' @return A matrix of heritability.
##
compute_heritability_multi <- function(traitTreeVar, samplingVar) {
  ## Both
  bothVars <- rbind(traitTreeVar, samplingVar)
  ## Function
  meanHeight_heritability <- function(both) {
    treeVar <- matrix(both[1:(dimTrait * dimTrait)], dimTrait, dimTrait)
    samplingVar <- matrix(both[(dimTrait * dimTrait + 1): (2 * dimTrait * dimTrait)], dimTrait, dimTrait)
    # return(treeVar %*% solve(treeVar + samplingVar))
    return(cor_extended(treeVar, samplingVar))
  }
  h_meanHeight <- t(apply(bothVars, 2, meanHeight_heritability))
  return(h_meanHeight)
}

cor_extended <- function(A, B) {
  C <- A
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      C[i, j] <- C[i, j] / sqrt((A[i, i] + B[i, i]) * (A[j, j] + B[j, j]))
    }
  }
  return(C)
}

get_sampling_variance <- function(vv, dimTrait, ntaxa, times_shared, treeScaled){
  samplingVar <- as.vector(solve(matrix(vv[grep("sampling.precision.matrix", names(vv))], ncol = dimTrait)))
  if (treeScaled) samplingVar <- mean(diag(times_shared)[1:ntaxa]) * samplingVar
  return(samplingVar)
}

##
#' @title Compute Multivariate heritability for BM
#'
#' @description
#'
#' @param df a log dataframe containing the "precision.matrix" and the
#' "sampling.precision" for each state of the MCMC.
#' @param meanHeight the mean height of the tree
#' @param dimTrait the dimension of the trait
#'
#' @return A matrix of heritability for each trait and for each state of the MCMC.
##
compute_heritability_multi_bm <- function(df, times_shared, dimTrait, ntaxa, sampleSizeRoot, treeScaled = FALSE) {
  get_bm_variance <- function(p, variance, sampleSizeRoot) {
    if (is.finite(sampleSizeRoot)){
      paramsBEAST <- params_BM(p = p,
                               variance = variance,
                               random = TRUE,
                               exp.root = rep(0, dimTrait),
                               var.root = variance / sampleSizeRoot)
    } else {
      paramsBEAST <- params_BM(p = p,
                               variance = variance,
                               random = FALSE,
                               value.root = rep(0, dimTrait))
    }
    V <-PhylogeneticEM:::compute_variance_block_diagonal.BM(times_shared, paramsBEAST, ntaxa)
  }
  sum_vars <- function(z) {
    vv <- get_bm_variance(dimTrait, matrix(z, dimTrait, dimTrait), sampleSizeRoot)
    apply(vv, c(1, 2), mean)
  }
  ## Trait tree variance
  colPrecMat <- grep("varianceMatrix", colnames(df))
  if (length(colPrecMat) == 0) colPrecMat <- grep("variance.matrix", colnames(df))
  traitTreeVar <- apply(df[, colPrecMat], 1, function(z) matrix(z, dimTrait, dimTrait))
  traitTreeVar <- apply(traitTreeVar, 2, sum_vars)
  ## Sampling variance
  samplingVar <- apply(df, 1, get_sampling_variance, dimTrait = dimTrait,ntaxa = ntaxa, times_shared = times_shared, treeScaled = treeScaled)
  ## heritability
  return(compute_heritability_multi(traitTreeVar, samplingVar))
}

##
#' @title Compute Multivariate heritability for BM, empirical definition from Hassler et al 2019
#'
#' @description
#'
#' @param df a log dataframe containing the "precision.matrix" and the
#' "sampling.precision" for each state of the MCMC.
#' @param meanHeight the mean height of the tree
#' @param dimTrait the dimension of the trait
#'
#' @return A matrix of heritability for each trait and for each state of the MCMC.
##
compute_heritability_multi_bm_empirical <- function(df, times_shared, dimTrait, ntaxa, sampleSizeRoot, treeScaled = FALSE) {
  ## Trait tree variance
  colPrecMat <- grep("varianceMatrix", colnames(df))
  if (length(colPrecMat) == 0) colPrecMat <- grep("variance.matrix", colnames(df))
  traitTreeVar <- apply(df[, colPrecMat], 1, function(z) matrix(z, dimTrait, dimTrait))
  traitTreeVar <- (sum(diag(times_shared)) / ntaxa - sum(times_shared) / (ntaxa^2)) * traitTreeVar
  ## Sampling variance
  samplingVar <- apply(df, 1, get_sampling_variance, dimTrait = dimTrait, ntaxa = ntaxa, times_shared = times_shared, treeScaled = treeScaled)
  samplingVar <- (ntaxa - 1) / ntaxa * samplingVar
  ## heritability
  return(compute_heritability_multi(traitTreeVar, samplingVar))
}

##
#' @title Compute Multivariate Heritability
#'
#' @description
#'
#' @param dimTrait the dimension of the trait
#' @param variance the variance of the OU
#' @param selectionStrength the selection strength of the OU
#' @param sampleSizeRoot the root sample size (Inf for fixed root)
#' @param times_shared
#' @param ntaxa
#'
#' @return A matrix of heritability.
##
compute_heritability_multi_ou_model <- function(dimTrait, variance, selectionStrength, sampleSizeRoot, samplingVar, times_shared, ntaxa) {
  get_ou_variance <- function(dimTrait, variance, selectionStrength, sampleSizeRoot) {
    if (anyNA(selectionStrength)) selectionStrength <- matrix(0, dimTrait, dimTrait)
    if (is.finite(sampleSizeRoot)){
      paramsBEAST <- params_OU(p = dimTrait, variance = variance,
                               selection.strength = selectionStrength,
                               optimal.value = rep(0, dimTrait),
                               random = TRUE, stationary.root = FALSE,
                               exp.root = rep(0, dimTrait), var.root = variance / sampleSizeRoot)
    } else {
      paramsBEAST <- params_OU(p = dimTrait, variance = variance,
                               selection.strength = selectionStrength,
                               optimal.value = rep(0, dimTrait),
                               random = FALSE, stationary.root = FALSE,
                               value.root = rep(0, dimTrait))
    }
    V <- PhylogeneticEM:::compute_variance_block_diagonal.OU(times_shared, paramsBEAST, ntaxa)
  }
  vv <- get_ou_variance(dimTrait, variance, selectionStrength, sampleSizeRoot)
  treeVar <- apply(vv, c(1, 2), mean)
  # return(treeVar %*% solve(treeVar + samplingVar))
  cor_extended(treeVar, samplingVar)
}

##
#' @title Compute Multivariate Heritability
#'
#' @description
#'
#' @param df a log dataframe containing the "precisionMatrix", the
#' "sampling.precision" and the "attenuation.values" for each state of the MCMC.
#' @param meanHeight the mean height of the tree
#' @param dimTrait the dimension of the trait
#'
#' @return A matrix of heritability for each trait and for each state of the MCMC.
##
compute_heritability_multi_ou <- function(df, times_shared, dimTrait, ntaxa, sampleSizeRoot, treeScaled = FALSE) {
  sum_vars <- function(z) {
    var <- matrix(as.numeric(z[colPrecMat]), dimTrait, dimTrait)
    att <- matrix(as.numeric(z[colAttMat]), dimTrait, dimTrait)
    samplingVar <- matrix(get_sampling_variance(z, dimTrait, ntaxa = ntaxa, times_shared = times_shared, treeScaled = treeScaled), dimTrait, dimTrait)
    compute_heritability_multi_ou_model(dimTrait, var, att, sampleSizeRoot, samplingVar, times_shared, ntaxa)
  }
  ## Trait tree variance
  colPrecMat <- grep("variance.matrix", colnames(df))
  if (length(colPrecMat) == 0) colPrecMat <- grep("varianceMatrix", colnames(df))
  colAttMat <- grep("attenuation.matrix", colnames(df))
  if (length(colAttMat) == 0) colAttMat <- grep("attenuationMatrix", colnames(df))
  return(t(apply(df, 1, sum_vars)))
  ## Sampling variance
  # samplingVar <- apply(df[, grep("sampling.precision", colnames(df))], 1, function(z) 1 / z)
  ## heritability
  # return(compute_heritability_multi(traitTreeVar, samplingVar))
}

##
#' @title Compute Multivariate Mean Correlation
#'
#' @description
#' Computes the mean correlation over all the tips, as expected by the model.
#'
#' @param dimTrait the dimension of the trait
#' @param variance the variance of the OU
#' @param selectionStrength the selection strength of the OU
#' @param sampleSizeRoot the root sample size (Inf for fixed root)
#' @param times_shared
#' @param ntaxa
#'
#' @return A matrix of heritability.
##
compute_mean_correlation_multi_ou_model <- function(dimTrait, variance, selectionStrength, sampleSizeRoot, samplingVar, times_shared, ntaxa) {
  get_ou_variance <- function(dimTrait, variance, selectionStrength, sampleSizeRoot) {
    if (anyNA(selectionStrength)) selectionStrength <- matrix(0, dimTrait, dimTrait)
    if (is.finite(sampleSizeRoot)){
      paramsBEAST <- params_OU(p = dimTrait, variance = variance,
                               selection.strength = selectionStrength,
                               optimal.value = rep(0, dimTrait),
                               random = TRUE, stationary.root = FALSE,
                               exp.root = rep(0, dimTrait), var.root = variance / sampleSizeRoot)
    } else {
      paramsBEAST <- params_OU(p = dimTrait, variance = variance,
                               selection.strength = selectionStrength,
                               optimal.value = rep(0, dimTrait),
                               random = FALSE, stationary.root = FALSE,
                               value.root = rep(0, dimTrait))
    }
    V <- PhylogeneticEM:::compute_variance_block_diagonal.OU(times_shared, paramsBEAST, ntaxa)
  }
  vv <- get_ou_variance(dimTrait, variance, selectionStrength, sampleSizeRoot)
  corMats <- apply(vv, 3, function(x) cor_sum(x, samplingVar))
  corMean <- apply(corMats, 1, mean)
  return(matrix(corMean, dimTrait, dimTrait))
  # cov2cor(treeVar + samplingVar)
}

cor_sum <- function(A, B) {
  C <- A
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      C[i, j] <- (A[i, j] + B[i, j]) / sqrt((A[i, i] + B[i, i]) * (A[j, j] + B[j, j]))
    }
  }
  return(C)
}

##
#' @title Compute Multivariate Mean Tip Correlation
#'
#' @description
#' Computes the mean correlation over all the tips, as expected by the model.
#'
#' @param df a log dataframe containing the "precisionMatrix", the
#' "sampling.precision" and the "attenuation.values" for each state of the MCMC.
#' @param meanHeight the mean height of the tree
#' @param dimTrait the dimension of the trait
#'
#' @return A matrix of heritability for each trait and for each state of the MCMC.
##
compute_mean_correlation_multi_ou <- function(df, times_shared, dimTrait, ntaxa, sampleSizeRoot, treeScaled = FALSE) {
  sum_vars <- function(z) {
    var <- matrix(as.numeric(z[colPrecMat]), dimTrait, dimTrait)
    att <- matrix(as.numeric(z[colAttMat]), dimTrait, dimTrait)
    samplingVar <- matrix(get_sampling_variance(z, dimTrait, ntaxa = ntaxa, times_shared = times_shared, treeScaled = treeScaled), dimTrait, dimTrait)
    compute_mean_correlation_multi_ou_model(dimTrait, var, att, sampleSizeRoot, samplingVar, times_shared, ntaxa)
  }
  ## Trait tree variance
  colPrecMat <- grep("variance.matrix", colnames(df))
  if (length(colPrecMat) == 0) colPrecMat <- grep("varianceMatrix", colnames(df))
  colAttMat <- grep("attenuation.matrix", colnames(df))
  if (length(colAttMat) == 0) colAttMat <- grep("attenuationMatrix", colnames(df))
  return(t(apply(df, 1, sum_vars)))
  ## Sampling variance
  # samplingVar <- apply(df[, grep("sampling.precision", colnames(df))], 1, function(z) 1 / z)
  ## heritability
  # return(compute_heritability_multi(traitTreeVar, samplingVar))
}

##
#' @title Sample trees from the posterior
#'
#' @description Take a trees log file, and write a new one with Ntree
#' sampled trees at random without replacement.
#'
#' @param tree_log a trees log file.
#' @param burning Percentage of trees to burn (default to 0.1)
#' @param Ntrees the number of trees to be sampled
#'
#' @return The name of the new log file with Ntree sampled trees.
##
sample_trees <- function(tree_log, Ntrees, burning = 0.1) {
  ## File Name
  name_new_file <- sub(".trees", paste0("_sample_", Ntrees, ".trees"), tree_log)
  if (file.exists(name_new_file)) {
    warning("A sample file already exists, cannot override it. Not doing anything.")
    return(name_new_file)
  }
  ## Read
  beast_trees <- readLines(tree_log)
  ## States
  states <- grep("STATE_", beast_trees)
  not_states <- !(1:length(beast_trees) %in% states)
  ## Burning
  Nburned <- ceiling(burning * length(states))
  burned_trees <- c(rep(TRUE, Nburned), rep(FALSE, length(states) - Nburned))
  ## Sampling
  sampled_states <- sample(states[!burned_trees], Ntrees)
  not_sampled_states <- !(1:length(beast_trees) %in% sampled_states) & !not_states
  beast_trees <- beast_trees[!not_sampled_states]
  ## Writing
  write(beast_trees, name_new_file)
  return(name_new_file)
}

##
#' @title Sub-sample trees
#'
#' @description Take a trees log file, and write a new one with trees
#' sampled every lagSample
#'
#' @param tree_log a trees log file.
#' @param lag the number of trees to be sampled
#'
#' @return The name of the new log file with sub sampled trees.
##
sub_sample_trees <- function(tree_log, lagSample = 1) {
  ## File Name
  name_new_file <- sub(".trees", paste0("_subsample_", lagSample, ".trees"), tree_log)
  if (file.exists(name_new_file)) {
    warning("A sample file already exists, cannot override it. Not doing anything.")
    return(name_new_file)
  }
  ## Read
  beast_trees <- readLines(tree_log)
  ## Sub-Sampling
  sub_sample_states <- sapply(beast_trees, is_multiple_lag, lagSample = lagSample)
  beast_trees <- beast_trees[sub_sample_states]
  ## Writing
  write(beast_trees, name_new_file)
  return(name_new_file)
}

is_multiple_lag <- function(tree_log_line, lagSample) {
  if (!grepl("STATE_", tree_log_line)) return(TRUE);
  tree_log_line <- sub("\\s\\[.*$", "", tree_log_line)
  tree_log_line <- sub("tree STATE_", "", tree_log_line)
  return(as.numeric(tree_log_line) %% lagSample == 0)
}

