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

## Working directory should be "HMC_OU"
library(here)
source(here("R_Utility_Files/read_log_functions.R"))

###############################################################################
## Parameters

## Directory
directory <- "exampleDataSets/TIMING"

datestamp_day <- "2020-09-03"

results_folder <- paste0("results_", datestamp_day)

## Utility functions
## Function to get all esses
get_esses <- function(fileName, iter) {
  ## Read log
  burnin <- 0.1
  if (grepl("Schnitzler2017_hmc", fileName))  burnin <- 0.0002
  if (grepl("blanquartMSM_hmc", fileName))  burnin <- 0.01
  dat_log <- read_log(here(directory, results_folder, paste0(fileName, "_", iter, ".log")), burning = burnin)
  time_step <- dat_log$state[2] - dat_log$state[1]
  ## Keep parameters only
  dat_log <- dat_log[, 5:ncol(dat_log)]
  ## get esses
  esses <- apply(dat_log, 2, function(trace) tracerer::calc_ess(trace, time_step))
  return(esses)
}

get_time_minutes <- function(fileName, iter) {
  console_file <- readLines(here(directory, results_folder, paste0(fileName, "_", iter, "_console.log")))
  timing <- console_file[grep("Time Chain:", console_file) + 1]
  timing <- strsplit(timing, " ")
  time <- as.double(timing[[1]][1])
  unit <- timing[[1]][2]
  if (unit == "minutes") return(time)
  if (unit == "hours") return(time * 60)
  return(NA)
}

get_esses_per_minutes <- function(fileName, iter) {
  esses <- get_esses(fileName, iter)
  time <- get_time_minutes(fileName, iter)
  return(esses / time)
}

res <- matrix(ncol = 4, nrow = 0)
res_all <- NULL
## Template file
for (fileName in c("Schnitzler2017", "blanquartMSM")) {
  for (model in c("OU_Cor", "BM_No")) {
    # fileName <- "blanquartMSM"
    # fileName <- "Schnitzler2017"
    
    mcmc_esses <- as.data.frame(t(sapply(1:10, function(iter) get_esses_per_minutes(paste0(fileName, "_", model), iter))))
    hmc_esses <- as.data.frame(t(sapply(1:10, function(iter) get_esses_per_minutes(paste0(fileName, "_hmc_", model), iter))))
    
    mcmc_mean_esses <- colMeans(mcmc_esses)
    hmc_mean_esses <- colMeans(hmc_esses)
    
    res_all[[paste0(fileName, "_", model)]] <- rbind(mcmc_mean_esses, hmc_mean_esses, hmc_mean_esses / mcmc_mean_esses)
    
    
    res_f <- data.frame(dataset = paste0(fileName, "_", model),
                        method = c("RW", "HMC", "Speed-up"),
                        median = c(median(mcmc_mean_esses),
                                   median(hmc_mean_esses),
                                   median(hmc_mean_esses) / median(mcmc_mean_esses)),
                        min = c(min(mcmc_mean_esses),
                                min(hmc_mean_esses),
                                min(hmc_mean_esses) / min(mcmc_mean_esses)))
    res <- rbind(res, res_f)
  }
}

## Summary result (median and mean)
res

## Complete results (all parameters)
res_all

## Median of all speed-ups
median(do.call(c, sapply(res_all, function(x) x[3, ])))
