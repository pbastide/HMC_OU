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

rm(list = ls())

library(here)

source(here("R_Utility_Files/write_xml_utility.R"))
source(here("R_Utility_Files/read_log_functions.R"))
source(here("R_Utility_Files/utility.R"))

library(PhylogeneticEM)

###############################################################################
## House keeping

## Directory
directory <- "exampleDataSets/BEEHIVE/simulations"

## Folder to keep the generated files
datestamp_day <- format(Sys.time(), "%Y-%m-%d")
new_folder <- paste0("xml_files_", datestamp_day)
dir.create(file.path(directory, new_folder))
## Folder for results
results_folder <- paste0("results_", datestamp_day)
dir.create(file.path(directory, results_folder))

## template xml file
fileName <- "blanquartMSM_sim"
template_file <- paste0(fileName, "_template_hmc.xml")

###############################################################################
## Parameters

## Tree MSM
treeMSM <- read.tree(file = "exampleDataSets/BEEHIVE/fixedTree/blanquartMSM.tree")
tree <- treeMSM
# rescale tree
tree$edge.length <- tree$edge.length / max(diag(vcv(tree))) # 53.0911

ntaxa <- length(tree$tip.label)

times_shared <- compute_times_ca(tree)

## Real data
# load(file = "exampleDataSets/BEEHIVE/fixedTree/data_blanquart/MSM_trait_beast.RData")

## Trait Name
traitName <- "traitSim"
dimTrait <- 2

## Parameters
alpha <- c(log(2)/0.3, 0)
opt <- c(1, -0.1)
sigma2 <- c(2 * alpha[1] * 0.1, 0.01)
r <- -0.9
var_e <- c(sigma2[1] / 2 / alpha[1], sigma2[2])
factor_noise <- c(0.1, 0.25, 0.5, 0.75, 1, 2, 3)

sampleSizeRoot <- Inf

nsim <- 1:2

sigma_mat <- diag(sigma2)
sigma_mat[1, 2] <- sigma_mat[2, 1] <- r * sqrt(prod(sigma2))

###############################################################################
## Heritability

compute_heritability <- function(ff) {
  var_e_f <- ff * var_e
  herrTrue <- compute_heritability_multi_ou_model(dimTrait, sigma_mat, diag(alpha), sampleSizeRoot, diag(var_e_f), times_shared, ntaxa)
  herrTrue <- data.frame(noise = ff, h1 = herrTrue[1, 1], h2 = herrTrue[2, 2], states = 1, method = "true")
  return(herrTrue)
}

herrTrue <- do.call(rbind, lapply(factor_noise, compute_heritability))
save(tree, ntaxa, times_shared,
     traitName, dimTrait,
     alpha, opt, sigma2, r, var_e, factor_noise, sampleSizeRoot, nsim, sigma_mat,
     herrTrue,
     file = here(directory, paste0(datestamp_day, "_true_parameters.RData")))


###############################################################################
## Simulations
set.seed(20110721)

for (ff in factor_noise) {

  var_e_f <- ff * var_e

  herrTrue <- compute_heritability_multi_ou_model(dimTrait, sigma_mat, diag(alpha), sampleSizeRoot, diag(var_e_f), times_shared, ntaxa)
  herrTrue <- data.frame(h1 = herrTrue[1, 1], h2 = herrTrue[2, 2], states = 1, method = "true")

  ## Sim with PhylogeneticEM
  params_sim <- params_OU(p = dimTrait,
                          variance = sigma_mat,
                          selection.strength = diag(alpha),
                          random = FALSE, stationary.root = FALSE,
                          value.root = opt,
                          optimal.value = opt)
  params_sim$process <- "OUBM"

  for (nrep in nsim) {
    sims <- simul_process(params_sim, phylo = tree)
    data_sim <- t(extract(sims, where = "tips", what = "states"))
    data_sim <- add_errors(data_sim, diag(var_e_f))

    colnames(data_sim) <- c("GSVL", "CD4")

    ###############################################################################
    ## Write to file

    ## Template
    xml_file <- file.path(directory, template_file)
    xml_file <- readLines(xml_file)

    ## Log Name
    xml_file <- sub("fileName=\"", paste0("fileName=\"", results_folder, "/"), xml_file)
    xml_file <- sub("_template", paste0("_bm_noise_", ff * 100, "_rep_", nrep), xml_file)

    ## Data
    for (taxon in tree$tip.label) {
      value <- data_sim[taxon, ]
      xml_file <- insert_below(xml_file,
                               paste0("\t\t\t<attr name=\"", traitName, "\">", paste(value, collapse = " "), "</attr>"),
                               grep(paste0("<taxon id=\"", taxon, "\">"), xml_file) + 1)
    }

    ## Tree
    posTree <- grep("<newick id=\"startingTree\">", xml_file)
    xml_file <- replace_line(xml_file,
                             write.tree(phy = tree),
                             posTree + 1)
    xml_file <- replace_line(xml_file,
                             "<newick id=\"startingTree\" usingDates=\"false\" usingHeights=\"true\">",
                             posTree)

    ## Write
    write(xml_file, file.path(directory, new_folder, paste0(fileName, "_bm_noise_", ff * 100, "_rep_", nrep, ".xml")))

  }
}

## Inference xml for OU models
source(here("exampleDataSets/BEEHIVE/simulations/write_xml_hmc.R"))

## Inference xml for models with added observation errors
source(here("exampleDataSets/BEEHIVE/simulations/write_xml_factor_hmc.R"))
