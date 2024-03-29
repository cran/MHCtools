#' ReplMatch() function
#'
#' In amplicon filtering it is sometimes valuable to compare technical
#' replicates in order to estimate the accuracy of a genotyping experiment. This
#' may be done both to optimize filtering settings and to estimate repeatability
#' to report in a publication. \code{\link{ReplMatch}} is designed to
#' automatically compare technical replicates in an amplicon filtering data set
#' and report the proportion of mismatches. The functions GetReplTable() and
#' GetReplStats() are designed to evaluate the output files.
#'
#' Note: ReplMatch() will throw a warning if all samples in a replicate set have
#' 0 sequences. In that case, the mean_props for that replicate set and the
#' repeatability for the data set will be NaN, and ReplMatch() will report which
#' replicate set is problematic and suggest to remove it from the repl_table. If
#' removing replicate sets, beware that the replicate sets in repl_table must be
#' numbered consecutively beginning at 1.
#'
#' If you publish data or results produced with MHCtools, please cite both of
#' the following references:
#' Roved, J. 2022. MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J., Hansson, B., Stervander, M., Hasselquist, D., & Westerdahl, H. 2022.
#' MHCtools - an R package for MHC high-throughput sequencing data: genotyping,
#' haplotype and supertype inference, and downstream genetic analyses in non-model
#' organisms. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13645
#'
#' @param repl_table is a table containing the sample names of technical
#'   replicates in the data set. This table should be organized so that the
#'   individual names are in the first column (Sample_ID), and the index number
#'   of the replicate set is in the second column (Replic_set). Replicate sets
#'   may contain more than two replicates, but sets must be numbered consecutively
#'   beginning at 1.
#' @param seq_table seq_table is a sequence table as output by the 'dada2'
#'   pipeline, which has samples in rows and nucleotide sequence variants in
#'   columns.
#' @param path_out is a user defined path to the folder where the output files
#'   will be saved.
#' @return  A set of R lists containing for each replicate set the observed
#'   sequence variants, the names of the sequences that were incongruent in the
#'   replicates, and the mean proportion of incongruent sequences (if 100%
#'   matches are expected between the replicates, this is equivalent of an error
#'   rate in the sequencing process). The sequences are named in the output by
#'   an index number corresponding to their column number in the sequence table,
#'   thus identical sequences will have identical sample names in all the output
#'   files. These files can be reopened in R e.g. using the readRDS() function
#'   in the base package.
#' @seealso \code{\link{GetReplTable}}; \code{\link{GetReplStats}}; for more
#'   information about 'dada2' visit <https://benjjneb.github.io/dada2/>
#' @examples
#' repl_table <- replicates_table
#' seq_table <- sequence_table_repl
#' path_out <- tempdir()
#' ReplMatch(repl_table, seq_table, path_out)
#' @export

ReplMatch <- function(repl_table, seq_table, path_out) {

  # The dada2 sequence table does not use sequence names, but identifies
  # sequence variants by their nuceotide sequence. Here I create a vector for
  # naming the sequences by their column number in the seq_table

  seq_names <- vector("character", length=length(colnames(seq_table)))

  seq_names <- paste0("Sequence_", formatC(seq(1:length(colnames(seq_table))), width = nchar(length(colnames(seq_table))), format = "d", flag = "0"))
  # the formatC() expression creates index numbers of the sequences with zeroes
  # padded in front, so that all numbers have the same number of digits. This is
  # necessary to avoid RegEx pattern matching, e.g. grepl matching "Sequence_1"
  # and any "Sequence_1X" or "Sequence_1XX".

  colnames(seq_table) <- seq_names

  # Define the number of replicate sets in the data set
  No_repl <- max(repl_table$Replic_set)

  # Loop over the rpelicate sets
  for (i in 1:No_repl) {

    Repl_samples <- factor()

    Repl_samples <- repl_table[repl_table[,"Replic_set"]==i,"Sample_ID"]

    # Create a list that will contain vectors for each replicate in the replicate set
    Repl_seqs <- vector("list", length(Repl_samples))


    # Loop over the samples in the replicate set
    for (j in 1:length(Repl_samples)) {

      # Create vectors Repl_seqs[[j]] which for each replicate will contain the
      # names of the sequences that are found in each sample
      Repl_seqs[[j]] <- factor()

      # Fetch column numbers for the sequences in replicate j in the sequence table
      z <- which(seq_table[paste(Repl_samples[j]),] > 0)

      # Put the names of the seqs in replicate j into Repl_seqs[[j]]
      Repl_seqs[[j]] <- seq_names[z]

    }


    # Create a list that will contain vectors for the (names of) sequences that
    # are found in either replicate but missing in at least one other
    Inc_seqs <- vector("list", length(Repl_samples))

    # Create a list that will contain vectors for the (names of) sequences that
    # are matched in all replicates
    Match_seqs <- vector("list", length(Repl_samples))

    # Loop over the samples in the replicate set
    for (j in 1:length(Repl_samples)) {

      Inc_seqs[[j]] <- factor()

      Match_seqs[[j]] <- factor()

      # Loop over the sequences in each sample
      for (Seq in Repl_seqs[[j]]) {

		# Match the sequences in replicate j against the sequences in the other replicates
        match_global <- grepl(Seq, Repl_seqs, perl = TRUE)

        if (sum(match_global) == length(match_global)) {

          Match_seqs[[j]] <- append(Match_seqs[[j]], Seq)

        } else {

          Inc_seqs[[j]] <- append(Inc_seqs[[j]], Seq)

        }

      }

    }


    ### Now the proportions of sequences that are incongruent in each replicate
    ### can be calculated, i.e. the sequences that are present in one replicate
    ### but missing in at least one other

    PrInc <- vector("numeric", length(Repl_samples))

    # Loop over the samples in the replicate set
    for (j in 1:length(Repl_samples)) {

      if(length(Repl_seqs[[j]]) > 0) PrInc[j] <- length(Inc_seqs[[j]])/length(Repl_seqs[[j]]) else PrInc[j] <- NA

    }

    ### If all samples in a replicate set has 0 sequences, i.e. if the length of
    ### all elements in Repl_seqs is 0, throw a warning that the mean_props for
    ### that replicate set and the repeatability for the data set will be NaN.

    if(sum(as.numeric(summary(Repl_seqs)[,1])) == 0) {warning(paste0("All samples in replicate set ", i, " have 0 sequences; mean_props will be NaN for this replicate set; global repeatability will be NaN. Suggestion: Remove this replicate set from repl_table."))}

    ### Finally, output a .Rds file with the mean proportion of incongruent
    ### sequences, the observed sequence variants, the matched sequence
    ### variants, and the incongruent sequence variants for each replicate set

    Mean_PrInc <- mean(PrInc,na.rm=TRUE)

    names(Repl_seqs) <- Repl_samples

    names(Inc_seqs) <- Repl_samples

    Output <- list(Mean_prop_incongr_seqs=c(Mean_PrInc), Matching_seqs=c(Match_seqs[[1]]), Observed_seqs=c(Repl_seqs), Incongruent_seqs=c(Inc_seqs))

    saveRDS(Output, file=paste0(path_out,"/Repl_set_", i, "_", c(format(Sys.Date(),"%Y%m%d")), ".Rds"))

  }

}
