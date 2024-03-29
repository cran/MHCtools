#' GetHpltTable() function
#'
#' \code{\link{GetHpltTable}} uses the output files produced by the HpltFind()
#' function to produce a table with the mean proportion of incongruent sequences
#' for each nest. If the mean proportion of incongruent sequences is generally
#' low, but certain nests have many incongruent sequences, biological reasons
#' may be causing the mismatches, e.g. extra-pair fertilizations or
#' recombination events.
#'
#' If you publish data or results produced with MHCtools, please cite both of
#' the following references:
#' Roved, J. 2022. MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J., Hansson, B., Stervander, M., Hasselquist, D., & Westerdahl, H. 2022.
#' MHCtools - an R package for MHC high-throughput sequencing data: genotyping,
#' haplotype and supertype inference, and downstream genetic analyses in non-model
#' organisms. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13645
#'
#' @param filepath is a user defined path to the folder where the output files
#'   from the HpltFind() function have been saved.
#' @return  A table with the mean proportion of incongruent sequences for each
#'   nest.
#' @seealso \code{\link{HpltFind}}; \code{\link{GetHpltStats}}
#' @examples
#' filepath <- system.file("extdata/HpltFindOut/", package="MHCtools")
#' GetHpltTable(filepath)
#' @export

GetHpltTable <- function(filepath) {

  # Get the file names of the .Rds output generated by the HpltFind function
  file_names <- dir(filepath)

  # Sort the file names by  nest number
  file_names <- file_names[order(as.numeric(gsub("[^0-9]", "", file_names)))]

  mean_props <- numeric()

  # Extract the mean proportion of incongruent sequences for each nest
  for(i in 1:length(file_names)) {

    mean_props[i] <- readRDS(file.path(filepath, file_names[i]))$Mean_prop_incongr_seqs

  }

  # List the observed mean proportions of incongruent sequences for each nest
  Table_obs_means <- as.data.frame(cbind(file_names, mean_props))

  Table_obs_means

}
