#' PapaDiv() function
#'
#' \code{\link{PapaDiv}} calculates the joint major histocompatibility complex
#' (MHC) diversity in parent pairs, taking into account alleles that are shared
#' between the parents. The joint diversity in parent pairs is often of interest
#' in studies of mate choice, fitness, and heritability.
#'
#' The PapaDiv() function outputs a set of R lists containing for the joint
#' diversity of each parent pair, the proportion of sequences that are shared
#' between the parents, the diversity of each of the parents, the observed sequence
#' variants in each parent, the matched sequence variants, and the incongruent
#' sequence variants in each parent.
#'
#' In addition, PapaDiv() produces a summary table with the names of the parents in
#' a pair, their respective MHC diversities, and the joint parent pair diversity.
#'
#' If you publish data or results produced with MHCtools, please cite both of
#' the following references:
#' Roved, J. 2022. MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J., Hansson, B., Stervander, M., Hasselquist, D., & Westerdahl, H. 2022.
#' MHCtools - an R package for MHC high-throughput sequencing data: genotyping,
#' haplotype and supertype inference, and downstream genetic analyses in non-model
#' organisms. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13645
#'
#' @param parents_table is a table containing the sample names of the parents in
#'   each nest. This table should be organized so that each row represents one
#'   nest, with the individual names of the mothers in the first column
#'   (Mother), and the individual names of the fathers in the second column
#'   (Father).
#' @param seq_table seq_table is a sequence table as output by the 'dada2'
#'   pipeline, which has samples in rows and nucleotide sequence variants in
#'   columns.
#' @param path_out is a user defined path to the folder where the output files
#'   will be saved.
#' @return  a set of R lists containing for the joint diversity of each parent
#'   pair, the proportion of sequences that are shared between the parents, the
#'   diversity of each of the parents, the observed sequence variants in each
#'   parent, the matched sequence variants, and the incongruent sequence
#'   variants in each parent. The sequences are named in the output by an index
#'   number corresponding to their column number in the sequence table, thus
#'   identical sequences will have identical sample names in all the output
#'   files. These files are saved in a sub folder in the output path called
#'   Parent_pairs (created by PapaDiv()) and can be reopened in R e.g. using
#'   the readRDS() function in the base package. For downstream data analysis,
#'   the PapaDiv() function also produces a summary table with the names of the
#'   parents in a pair, their respective MHC diversities, and the joint parent
#'   pair diversity. This table is saved as a .csv file in the output path.
#' @seealso For more information about 'dada2' visit
#'   <https://benjjneb.github.io/dada2/>
#' @examples
#' parents_table <- parents_table
#' seq_table <- sequence_table
#' path_out <- tempdir()
#' PapaDiv(parents_table, seq_table, path_out)
#' @importFrom "utils" "write.csv"
#' @export

PapaDiv <- function(parents_table, seq_table, path_out) {

  # The .Rds files will be saved in a sub folder in the output path called
  # "Parent_pairs".

  dir.create(paste0(path_out, "/Parent_pairs"))

  # The dada2 sequence table does not use sequences names, but identifies
  # sequence variants by their nuceotide sequence. Here I create a vector for
  # naming the sequences by their column number in the seq_table

  seq_names <- vector("character", length=length(colnames(seq_table)))

  seq_names <- paste0("Sequence_", formatC(seq(1:length(colnames(seq_table))), width = nchar(length(colnames(seq_table))), format = "d", flag = "0"))
  # the formatC() expression creates index numbers of the sequences with zeroes
  # padded in front, so that all numbers have the same number of digits. This is
  # necessary to avoid RegEx pattern matching, e.g. grepl matching "Sequence_1"
  # and any "Sequence_1X" or "Sequence_1XX".

  colnames(seq_table) <- seq_names

  ### Create two vectors that will contain the joint diversity of the parent
  ### pair and the proportion of sequences that are shared between the parents
  ### in the pair

  Joint_div <- vector("numeric", length=length(parents_table$Mother))

  Pr_match <- vector("numeric", length=length(parents_table$Mother))

  ### Create to vectors that will contain the diversity of each of the parents

  Mother_div <- vector("numeric", length=length(parents_table$Mother))

  Father_div <- vector("numeric", length=length(parents_table$Mother))


  for (i in 1:length(parents_table$Mother)) {

    # Create a list that will contain vectors for each parent in the parent pair

    Parent_seqs <- vector("list", length=2)

    for (j in 1:2) {

      # Create vectors Parent_seqs[[j]] which for each parent will contain the
      # names of the sequences that are found in each sample

      Parent_seqs[[j]] <- factor()

      # Fetch column numbers for the sequences in parent j in the sequence table

      z <- which(seq_table[paste(parents_table[i,j]),] > 0)

      # Put the names of the seqs in parent j into Parent_seqs[[j]]

      Parent_seqs[[j]] <- seq_names[z]

    }


    # Create a list that will contain vectors for the (names of) sequences that
    # are found in either parent but missing in the other

    Inc_seqs <- vector("list", length=2)

    # Create a list that will contain vectors for the (names of) sequences that
    # are matched in the parents

    Match_seqs <- vector("list", length=2)


    for (j in 1:2) {

      Inc_seqs[[j]] <- factor()

      Match_seqs[[j]] <- factor()

      # Match the sequences in parent j against the sequences in the other parent

      for (Seq in Parent_seqs[[j]]) {

        match_global <- grepl(Seq, Parent_seqs, perl = TRUE)

        if (sum(match_global) == length(match_global)) {

          Match_seqs[[j]] <- append(Match_seqs[[j]], Seq)

        } else {

          Inc_seqs[[j]] <- append(Inc_seqs[[j]], Seq)

        }

      }

    }


    ### Now the joint diversity of the parent pair and the proportion of
    ### sequences that are shared between the parents can be calculated

    Joint_div[i] <- length(Inc_seqs[[1]])+length(Inc_seqs[[2]])+length(Match_seqs[[1]])

    Pr_match[i] <- length(Match_seqs[[1]])/(length(Inc_seqs[[1]])+length(Inc_seqs[[2]])+length(Match_seqs[[1]]))

    ### Assign the diversity of each of the parents

    Mother_div[i] <- length(Parent_seqs[[1]])

    Father_div[i] <- length(Parent_seqs[[2]])


    ### Output a .Rds file with the joint diversity of the parent pair, the
    ### proportion of sequences that are shared between the parents, the
    ### diversity of each of the parents,  the observed sequence variants, the
    ### matched sequence variants, and the incongruent sequence variants for
    ### each parent pair

    names(Parent_seqs) <- c(paste(parents_table[i,1]), paste(parents_table[i,2]))

    names(Inc_seqs) <- c(paste(parents_table[i,1]), paste(parents_table[i,2]))

    Output <- list(Joint_pair_diversity=c(Joint_div[i]), Prop_matching_seqs=c(Pr_match[i]), Mother_diversity=c(Mother_div[i]), Father_diversity=c(Father_div[i]), Observed_seqs=c(Parent_seqs), Matching_seqs=c(Match_seqs[[1]]), Incongruent_seqs=c(Inc_seqs))

    saveRDS(Output, file=paste0(path_out, "/Parent_pairs/Parent_pair_", i, "_", c(format(Sys.Date(),"%Y%m%d")), ".Rds"))

  }

  ### Finally, output a .csv table that summarizes the joint diversity of the
  ### parent pairs, the diversity of each of the parents, and the proportion of
  ### sequences that are shared between the parents

  Summary_table <- as.data.frame(cbind(paste(parents_table[,1]), paste(parents_table[,2]), Mother_div, Father_div, Joint_div, Pr_match))

  rownames(Summary_table) <- seq(1:length(parents_table[,1]))
  colnames(Summary_table) <- c("Mother", "Father", "Mother_diversity", "Father_diversity", "Joint_pair_diversity", "Prop_matching_seqs")

  write.csv(Summary_table,file=paste0(path_out,"/Parent_pair_diversity_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))

}
