#' CreateFas() function
#'
#' \code{\link{CreateFas}} creates a FASTA file with all the sequences in a
#' 'dada2' sequence table.
#'
#' If you publish data or results produced with MHCtools, please cite both of
#' the following references:
#' Roved, J. (2022). MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J. (2024). MHCtools 1.5: Analysis of MHC sequencing data in R. In S.
#' Boegel (Ed.), HLA Typing: Methods and Protocols (2nd ed., pp. 275–295).
#' Humana Press. https://doi.org/https://doi.org/10.1007/978-1-0716-3874-3_18
#'
#' @param seq_table seq_table is a sequence table as output by the 'dada2'
#'   pipeline, which has samples in rows and nucleotide sequence variants in
#'   columns.
#' @param path_out is a user defined path to the folder where the output files
#'   will be saved.
#' @return A FASTA file with all the sequences in a 'dada2' sequence table. The
#'   sequences are named in the FASTA file by an index number corresponding to
#'   their column number in the sequence table.
#' @seealso \code{\link{CreateSamplesFas}}; for more information about 'dada2'
#'   visit <https://benjjneb.github.io/dada2/>
#' @examples
#' seq_table <- sequence_table_fas
#' path_out <- tempdir()
#' CreateFas(seq_table, path_out)
#' @export

CreateFas <- function(seq_table, path_out) {

  # Create empty file

  file.create(paste0(path_out, "/Sequences_", c(format(Sys.Date(),"%Y%m%d")), ".fas"))

  # for loop over all the sequences in the data set

  for (i in 1:length(colnames(seq_table))) {

    # Create sequence name line

    seq_name <- paste0(">Sequence_", formatC(i, width = nchar(length(colnames(seq_table))), format = "d", flag = "0"))
    # the formatC() expression creates index numbers of the sequences with zeroes
    # padded in front, so that all numbers have the same number of digits to prevent
    # RegEx pattern matching between e.g. "Sequence_1" and "Sequence_1X".

    # Extract nucleotide sequence from dada2 sequence table

    seq <- colnames(seq_table)[i]

    # Concatenate sequence name and nucleotide sequence in two lines

    lines <- paste(seq_name, seq, sep="\n")

    # Append sequence name line and nucleotide sequence to the file followed by
    # a line break

    cat(paste0(lines, "\n"), append = T, file = paste0(path_out, "/Sequences_", c(format(Sys.Date(),"%Y%m%d")), ".fas"))

  }

}



