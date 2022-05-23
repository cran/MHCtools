#' CreateSamplesFas() function
#'
#' \code{\link{CreateSamplesFas}} creates a set of FASTA files with the
#' sequences present in each sample in a 'dada2' sequence table.
#'
#' If you publish data or results produced with MHCtools, please cite both of
#' the following references:
#' Roved, J. 2022. MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J., Hansson, B., Stervander, M., Hasselquist, D., & Westerdahl, H. 2022.
#' MHCtools – an R package for MHC high‐throughput sequencing data: genotyping,
#' haplotype and supertype inference, and downstream genetic analyses in non‐model
#' organisms. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13645
#'
#' @param seq_table seq_table is a sequence table as output by the 'dada2'
#'   pipeline, which has samples in rows and nucleotide sequence variants in
#'   columns.
#' @param path_out is a user defined path to the folder where the output files
#'   will be saved.
#' @return  A set of FASTA files with the sequences present in each sample in
#'   the sequence table.The sequences are named in the FASTA files by an index
#'   number corresponding to their column number in the sequence table, thus
#'   identical sequences will have identical sample names in all the FASTA
#'   files.
#' @seealso \code{\link{CreateFas}}; for more information about 'dada2' visit
#'   <https://benjjneb.github.io/dada2/>
#' @examples
#' seq_table <- sequence_table_fas
#' path_out <- tempdir()
#' CreateSamplesFas(seq_table, path_out)
#' @export

CreateSamplesFas <- function(seq_table, path_out) {

  # The individual fasta files will be saved in a subfolder in the output path
  # called "Sample_fastas".

  dir.create(paste(path_out, "/Sample_fastas", sep=""))

  # Extract the sample names to a new vector

  sample_names <- rownames(seq_table)

  # for loop over all the samples in the data set

  for (i in 1:length(sample_names)) {

    # Create an empty file

    file.create(paste(path_out, "/Sample_fastas/", sample_names[i], "_", c(format(Sys.Date(),"%Y%m%d")), ".fas", sep=""))

    # Create a vector Sample_seqs that will contain the names of the nucleotide
    # sequences that are found in each sample

    Sample_seqs <- factor()

    # Fetch column numbers for the sequences in sample i in the dada2 sequence
    # table

    z <- which(seq_table[i,] > 0)

    # Put the nucleotide sequences found in sample i into Sample_seqs

    Sample_seqs <- colnames(seq_table[,z])

    ## Add the sequences in Sample_seqs to the fasta file for sample i

    # for loop over all the sequences in the data set

    for (j in 1:length(Sample_seqs)) {

      # Create sequence name line

      seq_name <- paste(">Sequence_", z[j], sep = "")

      # Extract nucleotide sequence from Sample_seqs

      seq <- Sample_seqs[j]

      # Concatenate sequence name and nucleotide sequence in two lines

      lines <- paste(seq_name, seq, sep="\n")

      # Append sequence name line and nucleotide sequence to the file followed
      # by a line break

      cat(paste(lines, "\n", sep = ""), append = T, file = paste(path_out, "/Sample_fastas/", sample_names[i], "_", c(format(Sys.Date(),"%Y%m%d")), ".fas", sep=""))

    }

  }

}

