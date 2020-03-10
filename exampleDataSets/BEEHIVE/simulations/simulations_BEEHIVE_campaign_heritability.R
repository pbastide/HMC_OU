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

source(here("R_Utility_Files/read_log_functions.R"))

library(PhylogeneticEM)

###############################################################################
## House keeping

datestamp_day <- "2019-11-29" ## Chage date here

## Directory
directory <- "exampleDataSets/BEEHIVE/simulations"
file_results <- paste0("exampleDataSets/BEEHIVE/fixedTree/results_", datestamp_day, "/blanquartMSM_var_2_sfr_")

## Folder for results
results_folder <- paste0("results_", datestamp_day)

###############################################################################
## Parameters

load(here(directory, paste0(datestamp_day, "_true_parameters.RData")))


###############################################################################
## Heritability

files_res <- list.files(file.path(directory, results_folder))
files_res <- files_res[grep("[0-9]*_hmc\\.log", files_res)]
files_res <- files_res[!grepl("_heritability", files_res)]
for (ff in files_res) {
  herr_file <- file.path(directory, results_folder, sub("\\.log", "_heritability.log", ff))
  if (!file.exists(herr_file)) {
    cat(paste0(herr_file, "\n"))
    dat_log <- read_log(file = file.path(directory, results_folder, ff), burning = 0)
    herr <- compute_heritability_multi_ou(dat_log, times_shared, dimTrait, ntaxa, sampleSizeRoot)
    colnames(herr) <- paste0(traitName, ".h", outer(1:dimTrait, 1:dimTrait, paste0))
    herr <- data.frame(herr)
    herr$state <- dat_log$state
    write.csv(herr, file = herr_file)
  }
}