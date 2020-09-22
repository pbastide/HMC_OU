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

#######
## Working directory should be "HMC_OU"
library(here)
## Utilities function
source(here("R_Utility_Files", "utilities_xml_functions.R"))
source(here("R_Utility_Files", "write_xml_utility.R"))

###############################################################################
## Parameters

traitName <- "traitSim"
dimTrait <- 2

datestamp_day <- format(Sys.time(), "%Y-%m-%d")
directory <- paste0("exampleDataSets/BEEHIVE/simulations/xml_files_", datestamp_day)

fileName <- "blanquartMSM_sim"

GSS <- TRUE
if (grepl("psss", fileName)) GSS <- FALSE
precvar <- "variance"
precvar_factor <- "variance"

samplingInit <- c(1, 0.01)

###############################################################################
## Models

for (nrep in 1:2) {
  for (factor_noise in c(0.1, 0.25, 0.5, 0.75, 1, 2, 3)) {
    for (model_name in c(
      "bm",
      "ou_diagonal",
      "oubm"
    )){

      model_name_file <- paste0(model_name, "_noise_", factor_noise * 100, "_rep_", nrep)

      for (scaled in c(FALSE, TRUE)) {
        for (correlated in c(FALSE, TRUE)) {
          scaled_string <- ""
          if (scaled) scaled_string <- "_scaled"

          correlated_string <- ""
          if (correlated) correlated_string <- "_cor"

          factor_name <- paste0(model_name, "_factor", correlated_string, scaled_string)
          factor_name <- paste0(factor_name, "_noise_", factor_noise * 100, "_rep_", nrep)

          ## Template file
          template_file <- paste0(fileName, "_", model_name_file, ".xml")
          fileNameModel <- paste0(fileName, "_", factor_name, ".xml")

          ## template xml file
          xml_file <- here(directory, template_file)
          xml_file <- readLines(xml_file)

          ###############################################################################
          ## Log names
          xml_file <- sub(paste0("_", model_name_file), paste0("_", factor_name), xml_file)

          ###############################################################################
          ## Node Trait paramter definitions
          posToWrite <- grep("</treeModel>", xml_file)

          nodeTraits <- c(
            paste0("\t\t<nodeTraits name=\"", traitName, ".tipTraits\" rootNode=\"false\" internalNodes=\"false\" leafNodes=\"true\" traitDimension=\"", dimTrait, "\" asMatrix=\"true\">"),# asMatrix=\"true\" initialValue=\"", traitValues, "\">"),
            paste0("\t\t\t<parameter id=\"", traitName, ".leafTraits\"/>"),
            "\t\t</nodeTraits>"
          )

          xml_file <- insert_above(xml_file, nodeTraits, posToWrite)

          ###############################################################################
          ## Factor Model

          if (precvar_factor == "precision") {
            block_sampling <- get_compound_symmetric_matrix(paste0(traitName, ".sampling.precision"),
                                                            valueDiag = samplingInit,
                                                            valueOffDiag = rep(0, dimTrait * (dimTrait - 1) / 2),
                                                            lowerDiag = rep(0, dimTrait),
                                                            indent = 3)
          } else if (precvar_factor == "variance") {
            block_sampling <- c(
              paste0("\t\t\t<cachedMatrixInverse id=\"", traitName, ".sampling.precision.matrix\">"),
              get_compound_symmetric_matrix(paste0(traitName, ".sampling.variance"),
                                            valueDiag = samplingInit,
                                            valueOffDiag = rep(0, dimTrait * (dimTrait - 1) / 2),
                                            lowerDiag = rep(0, dimTrait),
                                            indent = 4),
              "\t\t\t</cachedMatrixInverse>"
            )
          }

          scaleTrue <- ""
          if (scaled) scaleTrue <-"scaleByTipHeight=\"true\""

          block_factor <- c(
            paste0("\t<repeatedMeasuresModel id=\"", traitName, ".repeatedMeasures\" traitName=\"", traitName, "\" ", scaleTrue, ">"),
            "\t\t<treeModel idref=\"treeModel\"/>",
            "\t\t<traitParameter>",
            paste0("\t\t\t<parameter idref=\"", traitName, ".leafTraits\"/>"),
            "\t\t</traitParameter>",
            "\t\t<samplingPrecision>",
            block_sampling,
            "\t\t</samplingPrecision>",
            paste0("\t\t<multivariateDiffusionModel idref=\"", traitName, ".diffusionModel\"/>"),
            "\t</repeatedMeasuresModel>",
            ""
          )

          posToWrite <- grep("<traitDataLikelihood id=", xml_file)[1]
          xml_file <- insert_above(xml_file, block_factor, posToWrite)

          ## Factor trait parameter
          posToWrite <- grep("</traitDataLikelihood>", xml_file)
          for (pp in rev(posToWrite)) {
            xml_file <- insert_above(xml_file, paste0("\t\t<repeatedMeasuresModel idref=\"", traitName, ".repeatedMeasures\"/>"), pp)
          }

          posToDelete <- grep(paste0("<parameter id=\"leaf.", traitName, "\"/>"), xml_file)
          xml_file <- delete_lines(xml_file,
                                   posToDelete - 1,
                                   posToDelete + 1)

          ###############################################################################
          ## Sampling prior, operator and log

          ## Prior
          prior <- c(
            paste0("\t\t\t\t<halfTPrior id=\"", traitName, ".sampling.", precvar_factor, ".diagonal.prior\" df=\"1\" scale=\"2.5\">"),
            paste0("\t\t\t\t\t<parameter idref=\"", traitName,".sampling.", precvar_factor, ".diagonal\"/>"),
            "\t\t\t\t</halfTPrior>"
          )
          xml_file <- insert_above(xml_file, prior,
                                   grep("diffusionGradient id=", xml_file))

          xml_file <- insert_above(xml_file,
                                   paste0("\t\t\t\t<halfTPrior idref=\"", traitName, ".sampling.", precvar_factor, ".diagonal.prior\"/>"),
                                   grep("</prior>", xml_file))

          if (correlated) {
            prior <- c(
              paste0("\t\t\t\t<LKJCorrelationPrior id=\"", traitName, ".sampling.", precvar_factor, ".offDiagonal.prior\" shapeParameter=\"1.0\" dimension=\"2\">"),
              "\t\t\t\t\t<data>",
              paste0("\t\t\t\t\t<parameter idref=\"", traitName,".sampling.", precvar_factor, ".offDiagonal\"/>"),
              "\t\t\t\t\t</data>",
              "\t\t\t\t</LKJCorrelationPrior>"
            )
            xml_file <- insert_above(xml_file, prior,
                                     grep("diffusionGradient id=", xml_file))

            xml_file <- insert_above(xml_file,
                                     paste0("\t\t\t\t<LKJCorrelationPrior idref=\"", traitName, ".sampling.", precvar_factor, ".offDiagonal.prior\"/>"),
                                     grep("</prior>", xml_file))
          }

          ## Operator
          xml_file <- replace_line(xml_file,
                                   "\t<diffusionGradient id=\"diffusion.gradient\">",
                                   grep("<diffusionGradient id=", xml_file))
          # Gradient likelihood
          param_der <- "diagonal"
          if (correlated) param_der <- "both"
          xml_file <- insert_below(xml_file,
                                   c(paste0("\t<precisionGradient id=\"errorModel.gradient\"  parameter=\"", param_der, "\" traitName=\"", traitName, "\">"),
                                     paste0("\t\t<repeatedMeasuresModel idref=\"", traitName, ".repeatedMeasures\"/>"),
                                     paste0("\t\t<traitDataLikelihood idref=\"", traitName, ".traitLikelihood\"/>"),
                                     paste0("\t\t<compoundSymmetricMatrix idref=\"", traitName, ".sampling.precision.matrix\"/>"),
                                     "\t</precisionGradient>",
                                     "",
                                     paste0("\t<compoundGradient id=\"", traitName, ".traitLikelihood.gradient\">"),
                                     "\t\t<diffusionGradient idref=\"diffusion.gradient\"/>",
                                     "\t\t<precisionGradient idref=\"errorModel.gradient\"/>",
                                     "\t</compoundGradient>"
                                   ),
                                   grep("</diffusionGradient>", xml_file) + 1)
          # Gradient Prior
          if (grepl("drift", model_name)) {
            xml_file <- insert_below(xml_file,
                                     c("\t\t<gradient>",
                                       paste0("\t\t\t<halfTPrior idref=\"", traitName, ".sampling.", precvar_factor, ".diagonal.prior\"/>"),
                                       paste0("\t\t\t<parameter idref=\"", traitName, ".sampling.", precvar_factor, ".diagonal\"/>"),
                                       "\t\t</gradient>"),
                                     grep(paste0("<distributionLikelihood idref=\"", traitName, ".drift.prior\"/>"), xml_file)[1] + 2)
          } else {
            xml_file <- insert_below(xml_file,
                                     c("\t\t<gradient>",
                                       paste0("\t\t\t<halfTPrior idref=\"", traitName, ".sampling.", precvar_factor, ".diagonal.prior\"/>"),
                                       paste0("\t\t\t<parameter idref=\"", traitName, ".sampling.", precvar_factor, ".diagonal\"/>"),
                                       "\t\t</gradient>"),
                                     grep(paste0("<distributionLikelihood idref=\"", traitName, ".meanParameter.prior\"/>"), xml_file)[1] + 2)
          }
          if (correlated) {
            xml_file <- insert_below(xml_file,
                                     c("\t\t<gradient>",
                                       paste0("\t\t\t<LKJCorrelationPrior idref=\"", traitName, ".sampling.", precvar_factor, ".offDiagonal.prior\"/>"),
                                       paste0("\t\t\t<parameter idref=\"", traitName, ".sampling.", precvar_factor, ".offDiagonal\"/>"),
                                       "\t\t</gradient>"),
                                     grep(paste0("\t\t\t<halfTPrior idref=\"", traitName, ".sampling.", precvar_factor, ".diagonal.prior\"/>"), xml_file)[1] + 2)
          }
          if (grepl("oubm", model_name)) {
            # Mask
            posMask <- grep("<mask>", xml_file)[2] + 1
            xml_file <- replace_line(xml_file,
                                     gsub("\"/>", paste0(" ", paste(rep(1, dimTrait), collapse = " "), "\"/>"), xml_file[posMask]),
                                     posMask)
            # compound parameter
            xml_file <- insert_above(xml_file,
                                     paste0("\t\t\t\t<parameter idref=\"", traitName, ".sampling.", precvar_factor, ".diagonal\"/>"),
                                     grep("<multivariateCompoundTransform id=", xml_file) - 4)
            if (correlated) {
              # Mask
              posMask <- grep("<mask>", xml_file)[2] + 1
              xml_file <- replace_line(xml_file,
                                       gsub("\"/>", paste0(" ", paste(rep(1, dimTrait * (dimTrait - 1) / 2), collapse = " "), "\"/>"), xml_file[posMask]),
                                       posMask)
              # compound parameter
              xml_file <- insert_above(xml_file,
                                       paste0("\t\t\t\t<parameter idref=\"", traitName, ".sampling.", precvar_factor, ".offDiagonal\"/>"),
                                       grep("<multivariateCompoundTransform id=", xml_file) - 4)
            }
          } else {
            # compound parameter
            xml_file <- insert_above(xml_file,
                                     paste0("\t\t\t\t<parameter idref=\"", traitName, ".sampling.", precvar_factor, ".diagonal\"/>"),
                                     grep("<multivariateCompoundTransform id=", xml_file) - 1)
            if (correlated) {
              xml_file <- insert_above(xml_file,
                                       paste0("\t\t\t\t<parameter idref=\"", traitName, ".sampling.", precvar_factor, ".offDiagonal\"/>"),
                                       grep("<multivariateCompoundTransform id=", xml_file) - 1)
            }
          }
          # Transform
          xml_file <- insert_above(xml_file,
                                   c(paste0("\t\t\t\t<transform type=\"log\" dim=\"", dimTrait, "\"/>")),
                                   grep("</multivariateCompoundTransform>", xml_file))
          if (correlated) {
            xml_file <- insert_above(xml_file,
                                     c(paste0("\t\t\t\t<LKJTransform dimension=\"", dimTrait, "\"/>")),
                                     grep("</multivariateCompoundTransform>", xml_file))
          }

          ## Log
          posToWrite <- grep("</log>", xml_file)[2]
          xml_file <- insert_above(xml_file, paste0("\t\t\t<parameter idref=\"", traitName,".sampling.", precvar_factor, ".matrix\"/>"), posToWrite)
          xml_file <- insert_above(xml_file, paste0("\t\t\t<parameter idref=\"", traitName,".sampling.", precvar_factor, ".diagonal\"/>"), posToWrite)
          xml_file <- insert_above(xml_file, paste0("\t\t\t<parameter idref=\"", traitName,".sampling.", precvar_factor, ".offDiagonal\"/>"), posToWrite)
          xml_file <- insert_above(xml_file, paste0("\t\t\t<parameter idref=\"", traitName,".sampling.precision.matrix\"/>"), posToWrite)

          ########################################################################################
          ## GSS
          if (GSS) {
            # working prior
            posGSS <- grep(paste0("<logTransformedNormalReferencePrior id=\"", traitName, ".", precvar, ".diagonal.workingPrior\""), xml_file)
            block <- gsub(paste0(traitName, ".", precvar, ".diagonal"),
                          paste0(traitName, ".sampling.", precvar_factor, ".diagonal"),
                          xml_file[posGSS + 0:2])
            xml_file <- insert_below(xml_file, block,
                                     grep(paste0("<normalReferencePrior id=\"", traitName, ".meanParameter.workingPrior\""), xml_file) + 2)
            if (correlated) {
              posGSS <- grep(paste0("<normalReferencePrior id=\"", traitName, ".", precvar, ".offDiagonal.workingPrior\""), xml_file)
              block <- gsub(paste0(traitName, ".", precvar, ".offDiagonal"),
                            paste0(traitName, ".sampling.", precvar_factor, ".offDiagonal"),
                            xml_file[posGSS + 0:2])
              xml_file <- insert_below(xml_file, block,
                                       grep(paste0("<logTransformedNormalReferencePrior id=\"", traitName, ".sampling.", precvar_factor, ".diagonal.workingPrior\""), xml_file) + 2)
            }
            # working prior gradient
            posGSSGradient <- grep(paste0("<gradient id=\"", traitName, ".", precvar, ".diagonal.workingPrior.gradient\""), xml_file)
            blockGradient <- gsub(paste0(traitName, ".", precvar, ".diagonal"),
                                  paste0(traitName, ".sampling.", precvar_factor, ".diagonal"),
                                  xml_file[posGSSGradient + 0:3])
            if (grepl("drift", model_name)) {
              xml_file <- insert_below(xml_file, blockGradient,
                                       grep(paste0("<gradient id=\"", traitName, ".drift.workingPrior.gradient\""), xml_file) + 3)
            } else {
              xml_file <- insert_below(xml_file, blockGradient,
                                       grep(paste0("<gradient id=\"", traitName, ".meanParameter.workingPrior.gradient\""), xml_file) + 3)
            }
            if (correlated) {
              posGSSGradient <- grep(paste0("<gradient id=\"", traitName, ".", precvar, ".offDiagonal.workingPrior.gradient\""), xml_file)
              blockGradient <- gsub(paste0(traitName, ".", precvar, ".offDiagonal"),
                                    paste0(traitName, ".sampling.", precvar_factor, ".offDiagonal"),
                                    xml_file[posGSSGradient + 0:3])
              xml_file <- insert_below(xml_file, blockGradient,
                                       grep(paste0("<gradient id=\"", traitName, ".sampling.", precvar_factor, ".diagonal.workingPrior.gradient\""), xml_file) + 3)
            }

            ###############################################################################
            ## Actual writing

            write(xml_file, here(directory, fileNameModel))
          }
        }
      }
    }
  }
}