## XML writing utility functions

##############################################################################
## Insert and replace

insert_above <- function(xml_file, block, linePos) {
  return(c(xml_file[1:(linePos - 1)],
           block,
           xml_file[linePos:length(xml_file)]))
}

insert_below <- function(xml_file, block, linePos) {
  return(c(xml_file[1:linePos],
           block,
           xml_file[(linePos + 1):length(xml_file)]))
}

replace_line <- function(xml_file, block, linePos) {
  return(c(xml_file[1:(linePos - 1)],
           block,
           xml_file[(linePos + 1):length(xml_file)]))
}

delete_lines <- function(xml_file, from, to) {
  return(c(xml_file[1:(from - 1)],
           xml_file[(to + 1):length(xml_file)]))
}

comment_lines <- function(xml_file, from, to) {
  xml_file <- insert_above(xml_file, "<!--", from)
  xml_file <- insert_below(xml_file, "-->", to)
  return(xml_file)
}

##############################################################################
## Log

write_in_log <- function(nameParam, pos) {
  return(insert_above(xml_file,
                      paste0("\t\t\t<parameter idref=\"", traitName,".", nameParam, "\"/>"),
                      grep("</log>", xml_file)[pos]))
}

write_in_log_prior <- function(nameParam, pos) {
  return(insert_above(xml_file,
                      paste0("\t\t\t<prior idref=\"", traitName,".", nameParam, ".prior\"/>"),
                      grep("</log>", xml_file)[pos]))
}

write_in_log_workingPrior <- function(nameParam, pos) {
  return(insert_above(xml_file,
                      paste0("\t\t\t<prior idref=\"", traitName,".", nameParam, ".workingPrior\"/>"),
                      grep("</log>", xml_file)[pos]))
}

write_in_log_traitDataLikelihood <- function(nameParam, pos) {
  return(insert_above(xml_file,
                      paste0("\t\t\t<traitDataLikelihood idref=\"", traitName,".", nameParam, "\"/>"),
                      grep("</log>", xml_file)[pos]))
}