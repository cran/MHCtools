#' CreateHpltOccTable() function
#'
#' \code{\link{CreateHpltOccTable}} is designed to create a haplotype-sequence
#' occurrence matrix from the set of R lists with putative haplotypes output
#' by the HpltFind() function. CreateHpltOccTable() assumes that data originated
#' from a diploid species.
#'
#' If you publish data or results produced with MHCtools, please cite both of
#' the following references:
#' Roved, J. 2022. MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J., Hansson, B., Stervander, M., Hasselquist, D., & Westerdahl, H. 2022.
#' MHCtools - an R package for MHC high-throughput sequencing data: genotyping,
#' haplotype and supertype inference, and downstream genetic analyses in non-model
#' organisms. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13645
#'
#' @param seq_table seq_table is a sequence table as output by the 'dada2'
#'   pipeline, which has samples in rows and nucleotide sequence variants in
#'   columns.
#' @param filepath is a user defined path to the folder where the output files
#'   from the HpltFind() function have been saved.
#' @param path_out is a user defined path to the folder where the output files
#'   will be saved.
#' @return  A binary (logical) occurrence matrix with the data set sequences
#'   (inherited from seq_table) in columns and the putative haplotypes inferred
#'   by the HpltFind() function in rows.
#' @seealso \code{\link{HpltFind}}; for more information about 'dada2' visit
#' <https://benjjneb.github.io/dada2/>
#' @examples
#' seq_table <- sequence_table
#' filepath <- system.file("extdata/HpltFindOut/", package="MHCtools")
#' path_out <- tempdir()
#' CreateHpltOccTable(seq_table, filepath, path_out)
#' @export

CreateHpltOccTable <- function(seq_table, filepath, path_out) {

  # Extract the sequence names from seq_table
  seq_names <- colnames(seq_table)

  # Get the file names of the .Rds output (created by HpltFind)
  file_names <- dir(filepath)

  # Sort the file names by  nest number
  file_names <- file_names[order(as.numeric(gsub("[^0-9]", "", file_names)))]

  # create a matrix
  hplt_occ_matrix <- matrix(data = NA, nrow = 4*length(file_names), ncol = length(seq_names))
  colnames(hplt_occ_matrix) <- seq_names
  rownames(hplt_occ_matrix) <- character(length=dim(hplt_occ_matrix)[1])

  # read in the .Rds files one by one and assign one line in the matrix to each haplotype
  for(i in 1:length(file_names)) {

    nest_i <- readRDS(file.path(filepath, file_names[i]))

    # loop over the haplotypes in nest_i
    for(j in 1:4) {

      # paste the haplotype name as rowname in the matrix together with the nest number
      rownames(hplt_occ_matrix)[(4*i-4+j)] <- paste0(names(nest_i$Putative_haplotypes[j]),"_nest_",i)

      # assign 1 to the columns where the colnames are found in the putative haplotype
      # assign 0 to all other columns
      hplt_occ_matrix[4*i-4+j,] <- as.numeric(seq_names %in% nest_i$Putative_haplotypes[[j]])

    }

  }

  # Save the table as .csv
  write.csv(hplt_occ_matrix, file=paste0(path_out,"/Hplt_occ_matrix_",c(format(Sys.Date(),"%Y%m%d")),".csv"))

}
