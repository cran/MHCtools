#' CreateFas() function
#'
#' \code{\link{CreateFas}} creates a FASTA file with all the sequences in a
#' 'dada2' sequence table.
#'
#' If you publish data produced with MHCtools, please cite:
#' Roved, J. 2020. MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J., Hansson, B., Stervander, M., Hasselquist, D., Westerdahl, H. 2020.
#' Non-random association of MHC-I alleles in favor of high diversity haplotypes
#' in wild songbirds revealed by computer-assisted MHC haplotype inference using
#' the R package MHCtools. bioRxiv.
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

  file.create(paste(path_out, "/Sequences_", c(format(Sys.Date(),"%Y%m%d")), ".fas", sep=""))

  # for loop over all the sequences in the data set

  for (i in 1:length(colnames(seq_table))) {

    # Create sequence name line

    seq_name <- paste(">Sequence_", i, sep = "")

    # Extract nucleotide sequence from dada2 sequence table

    seq <- colnames(seq_table)[i]

    # Concatenate sequence name and nucleotide sequence in two lines

    lines <- paste(seq_name, seq, sep="\n")

    # Append sequence name line and nucleotide sequence to the file followed by
    # a line break

    cat(paste(lines, "\n", sep = ""), append = T, file = paste(path_out, "/Sequences_", c(format(Sys.Date(),"%Y%m%d")), ".fas", sep=""))

  }

}



