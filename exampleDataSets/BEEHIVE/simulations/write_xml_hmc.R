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
## Utilities function
source(here("R_Utility_Files", "utilities_xml_functions.R"))
source(here("R_Utility_Files", "write_xml_utility.R"))

###############################################################################
## Parameters

## Directory
datestamp_day <- format(Sys.time(), "%Y-%m-%d")
directory <- paste0("exampleDataSets/BEEHIVE/simulations/xml_files_", datestamp_day)
directory_xml <- directory

## Template file
fileName <- "blanquartMSM_sim"
GSS <- TRUE
if (grepl("psss", fileName)) GSS <- FALSE
precvar <- "variance"
nrep <- NULL
factor_noise <- NULL
for (nrep in 1:2) {
  for (factor_noise in c(0.1, 0.25, 0.5, 0.75, 1, 2, 3)) {
    bm_name <- paste0("_bm")
    if (!is.null(nrep)) bm_name <- paste0(bm_name, "_noise_", factor_noise * 100, "_rep_", nrep)
    template_file <- paste0(fileName, bm_name, ".xml")

    ## Trait Name
    traitName <- "traitSim"
    dimTrait <- 2

    traitNameOU <- traitName
    traitNameDrift <- traitName
    dimTraitOU <- dimTrait
    dimTraitDrift <- dimTrait

    ###############################################################################
    ## Model

    for (model_name in c(
      "bm",
      "ou_diagonal",
      "oubm"
    )) {

      model_name <- paste0(model_name, "_noise_", factor_noise * 100, "_rep_", nrep)

      fileNameModel <- paste0(fileName, "_", model_name, ".xml")

      ###############################################################################
      ## template xml file
      xml_file <- here(directory_xml, template_file)
      xml_file <- readLines(xml_file)

      ###############################################################################
      ## Log names
      xml_file <- sub(bm_name, paste0("_", model_name), xml_file)


      ###############################################################################
      ## OU
      if (grepl("ou", model_name)) {

        traitName <- traitNameOU
        dimTrait <- dimTraitOU

        dimOU <- dimTrait
        masked <- ""
        canBeNeg <- grepl("acdc", model_name)
        lowerOU <- rep(0, dimTrait)
        if (canBeNeg) lowerOU <- NULL

        ## Declaration Attenuation
        posToWrite <- grep(paste0("<traitDataLikelihood id=\"", traitName), xml_file)

        if (grepl("ou_diagonal", model_name)) {
          att_mat <- get_diagonal_matrix(name = paste0(traitName, ".attenuation"),
                                         value = rep(1, dimTrait), lower = lowerOU)
          att_mat <- c(att_mat, "")
        } else if (grepl("oubm", model_name)) {
          att_mat <- get_compound_diagonal_matrix(name = paste0(traitName, ".attenuation"),
                                                  value = c(1, 0), lower = lowerOU)
          att_mat <- c(att_mat,
                       c(paste0("\t<maskedParameter id=\"", traitName, ".attenuation.values.masked", "\">"),
                         paste0("\t\t<parameter idref=\"", traitName, ".attenuation.values", "\"/>"),
                         "\t\t<mask>",
                         "\t\t\t<parameter value=\"1.0 0.0\"/>",
                         "\t\t</mask>",
                         "\t</maskedParameter>"))
          dimOU <- dimTrait - 1
          masked <- ".masked"
        }

        xml_file <- insert_above(xml_file, att_mat, posToWrite)

        ## Model Likelihood
        posLL <- grep(paste0("<traitDataLikelihood id=\"", traitName), xml_file)
        posToWrite <- grep("</conjugateRootPrior>", xml_file)
        posToWrite <- min(posToWrite[posToWrite > posLL])

        optVals <- as.vector(get_rate_blocs(rep(0, dimTrait), ids = paste0(traitName, ".opt.", 1:dimTrait), ref = TRUE))

        block <- c(
          paste0("\t\t<optimalTraits id=\"", paste0(traitName, ".opt"), "\">"),
          optVals,
          "\t\t</optimalTraits>",
          "\t\t<strengthOfSelectionMatrix>",
          paste0("\t\t\t<matrixParameter idref=\"", paste0(traitName, ".attenuation.matrix"), "\"/>"),
          "\t\t</strengthOfSelectionMatrix>"
        )

        xml_file <- insert_below(xml_file, block, posToWrite)

        ## Operators
        # Gradient likelihood attenuation
        posG <- grep(paste0("<diffusionGradient id=\"", traitName), xml_file)
        posToWrite <- grep("<meanGradient parameter=", xml_file)
        posToWrite <- min(posToWrite[posToWrite > posG])
        xml_file <- insert_above(xml_file,
                                 c(paste0("\t\t<attenuationGradient parameter=\"diagonal\" traitName=\"", traitName, "\">"),
                                   paste0("\t\t\t<traitDataLikelihood idref=\"", traitName, ".traitLikelihood\"/>"),
                                   paste0("\t\t\t<compoundSymmetricMatrix idref=\"", traitName, ".attenuation.matrix\"/>"),
                                   "\t\t</attenuationGradient>"),
                                 posToWrite)
        # Gradient likelihood optimum
        posToWrite <- grep("<meanGradient parameter=\"root\"", xml_file)
        posToWrite <- min(posToWrite[posToWrite > posG])
        xml_file <- replace_line(xml_file,
                                 paste0("\t\t<meanGradient parameter=\"both\" traitName=\"", traitName, "\">"),
                                 posToWrite)
        # Gradient Prior
        xml_file <- insert_below(xml_file,
                                 c("\t\t<gradient>",
                                   paste0("\t\t\t<distributionLikelihood idref=\"", traitName, ".attenuation.values.prior\"/>"),
                                   paste0("\t\t\t<parameter idref=\"", traitName, ".attenuation.values\"/>"),
                                   "\t\t</gradient>"),
                                 grep(paste0("<LKJCorrelationPrior idref=\"", traitName), xml_file)[1] + 1)
        # compound parameter
        xml_file <- insert_above(xml_file,
                                 paste0("\t\t\t\t<parameter idref=\"", traitName, ".attenuation.values\"/>"),
                                 grep(paste0("<multivariateCompoundTransform id=\"", traitName), xml_file) - 2)
        # Mask
        if (grepl("oubm", model_name)) {
          xml_file <- insert_above(xml_file,
                                   c("\t\t\t<mask>",
                                     get_parameter_line(value = c(rep(1, dimTrait * (dimTrait + 1) / 2 + dimOU), 0, rep(1, dimTrait)),
                                                        id = paste0(traitName, ".mask.parameter")),
                                     "\t\t\t</mask>"),
                                   grep(paste0("<multivariateCompoundTransform id=\"", traitName), xml_file))
        }
        # Transform
        if (grepl("oubm", model_name)) {
          xml_file <- insert_below(xml_file,
                                   c(paste0("\t\t\t\t<transform type=\"log\" dim=\"", dimOU, "\"/>"),
                                     paste0("\t\t\t\t<transform type=\"none\" dim=\"", 1, "\"/>")),
                                   grep("<LKJTransform dimension=", xml_file))
        } else {
          xml_file <- insert_below(xml_file,
                                   c(paste0("\t\t\t\t<transform type=\"log\" dim=\"", dimOU, "\"/>")),
                                   grep("<LKJTransform dimension=", xml_file))
        }

        ## Priors
        # attenuation values
        posP <- grep(paste0("<normalPrior id=\"", traitName), xml_file)
        posToWrite <- grep("</normalPrior>", xml_file)
        posToWrite <- min(posToWrite[posToWrite > posP])
        xml_file <- insert_below(xml_file,
                                 gsub("attenuation.values\"",
                                      paste0("attenuation.values", masked, "\""),
                                      get_prior_normal_bloc(paste0(traitName, ".attenuation.values"), mean = 0, stdev = 7.07306, half = !canBeNeg)),
                                 posToWrite)
        xml_file <- insert_below(xml_file,
                                 paste0("\t\t\t\t<halfNormalPrior idref=\"", traitName, ".attenuation.values.prior\"/>"),
                                 grep(paste0("<LKJCorrelationPrior idref=\"", traitName), xml_file)[2])


        ## Log
        # attenuation matrix
        xml_file <- write_in_log("attenuation.matrix", 2)

        # attenuation values
        xml_file <- write_in_log("attenuation.values", 2)
        xml_file <- write_in_log_prior("attenuation.values", 2)

        ## GSS
        if (GSS) {
          # attenuation values
          posGSS <- grep(paste0("<logTransformedNormalReferencePrior id=\"", traitName, ".", precvar, ".diagonal.workingPrior\""), xml_file)
          posGSSGradient <- grep(paste0("<gradient id=\"", traitName, ".", precvar, ".diagonal.workingPrior.gradient\""), xml_file)
          if (grepl("oubm", model_name)) {
            block <- NULL
            for (i in 1:dimOU) {
              block <- c(block, gsub(paste0(traitName, ".", precvar, ".diagonal"), paste0(traitName, ".attenuation.values.", i), xml_file[posGSS + 0:2]))
            }
            block <- gsub("dimension=\"[0-9]\"", "dimension=\"1\"", block)
            block <- gsub("attenuation.values.1.workingPrior", "attenuation.values.workingPrior", block)
            block <- gsub(paste0("idref=\"", traitName, ".attenuation.values.1"), paste0("idref=\"", traitName, ".attenuation.values.masked"), block)
          } else {
            # working prior
            block <- gsub(paste0(traitName, ".", precvar, ".diagonal"), paste0(traitName, ".attenuation.values"), xml_file[posGSS + 0:2])
          }
          blockGradient <- gsub(paste0(traitName, ".", precvar, ".diagonal"), paste0(traitName, ".attenuation.values"), xml_file[posGSSGradient + 0:3])

          xml_file <- insert_above(xml_file, block,
                                   grep(paste0("<normalReferencePrior id=\"", traitName, ".meanParameter.workingPrior\""), xml_file))
          xml_file <- insert_above(xml_file, blockGradient,
                                   grep(paste0("<gradient id=\"", traitName, ".meanParameter.workingPrior.gradient\""), xml_file))

          # Mask
          if (grepl("oubm", model_name)) {
            xml_file <- insert_above(xml_file,
                                     c("\t\t\t<mask>",
                                       get_parameter_line(value = NA, id = paste0(traitName, ".mask.parameter"), ref = TRUE),
                                       "\t\t\t</mask>"),
                                     grep(paste0("<multivariateCompoundTransform idref=\"", traitName), xml_file))
          }
        }

      }

      ###############################################################################
      ## Write the file
      write(xml_file, here(directory, fileNameModel))
    }
  }
}
