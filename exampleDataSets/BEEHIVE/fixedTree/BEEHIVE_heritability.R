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

## Working directory should be "HMC_OU"
library(here)

source(here("R_Utility_Files/read_log_functions.R"))

library(PhylogeneticEM)

###############################################################################
## House keeping

datestamp_day <- "2020-01-30" ## Change date here

## Directory
directory <- "exampleDataSets/BEEHIVE/fixedTree"
file_results <- paste0("exampleDataSets/BEEHIVE/fixedTree/results_", datestamp_day, "/blanquartMSM_no_miss_var_2_sfr_hmc_")

## Folder for results
results_folder <- paste0("results_", datestamp_day)

## Parameters
MSMtree <- read.tree(file = here("exampleDataSets/BEEHIVE/fixedTree/blanquartMSM.tree"))
MSMtree$edge.length <- MSMtree$edge.length / 53.0911
dimTrait <- 2
ntaxa <- length(MSMtree$tip.label)
sampleSizeRoot <- Inf

times_shared <- compute_times_ca(MSMtree)

traitName <- "trait2"


###############################################################################
## Heritability

files_res <- list.files(here(directory, results_folder))
files_res <- files_res[grep("\\.log", files_res)]
files_res <- files_res[grepl("_factor", files_res)]
files_res <- files_res[!grepl("_console", files_res)]
files_res <- files_res[!grepl("_lambda", files_res)]
files_res <- files_res[!grepl("_heritability", files_res)]
for (ff in files_res) {
  herr_file <- here(directory, results_folder, sub("\\.log", "_heritability.log", ff))
  if (!file.exists(herr_file)) {
    cat(paste0(herr_file, "\n"))
    tree_scaled <- grepl("_scaled", ff)
    dat_log <- read_log(file = here(directory, results_folder, ff), burning = 0)
    herr <- compute_heritability_multi_ou(dat_log, times_shared, dimTrait, ntaxa, sampleSizeRoot, treeScaled = tree_scaled)
    colnames(herr) <- paste0(traitName, ".h", outer(1:dimTrait, 1:dimTrait, paste0))
    herr <- data.frame(herr)
    herr$state <- dat_log$state
    write.csv(herr, file = herr_file)
  }
}