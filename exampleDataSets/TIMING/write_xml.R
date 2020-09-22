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

###############################################################################
## Parameters

## Directory
directory <- "exampleDataSets/TIMING"
directory_template <- "xml_templates"

datestamp_day <- format(Sys.time(), "%Y-%m-%d")

new_folder <- paste0("xml_files_", datestamp_day)
dir.create(here(directory, new_folder))

results_folder <- paste0("results_", datestamp_day)
dir.create(here(directory, results_folder))

## Template file
fileNames <- c("blanquartMSM", "blanquartMSM_hmc",
               "Schnitzler2017", "Schnitzler2017_hmc")

for (fileName in fileNames) {
  for (model in c("OU_Cor", "OU_Ind", "BM_No")) {
    for (iter in 1:10) {
      ###############################################################################
      ## template xml file
      fileNameModel <- paste0(fileName, "_", model)
      template_file <- paste0(fileNameModel, ".xml")

      xml_file <- here(directory, directory_template, template_file)
      xml_file <- readLines(xml_file)

      ###############################################################################
      ## Log names
      xml_file <- sub("fileName=\"", paste0("fileName=\"", results_folder, "/"), xml_file)
      xml_file <- sub(paste0(fileName, ".log\""), paste0(fileNameModel, "_", iter, ".log\""), xml_file)


      ###############################################################################
      ## Write the file
      write(xml_file, here(directory, new_folder, paste0(fileNameModel, "_", iter, ".xml")))
    }
  }
}
