#' HpltFind() function
#'
#' \code{\link{HpltFind}} is designed to automatically infer major
#' histocompatibility complex (MHC) haplotypes from the genotypes of parents and
#' offspring in families (defined as nests) in non-model species, where MHC
#' sequence variants cannot be identified as belonging to individual loci. The
#' functions GetHpltTable() and GetHpltStats() are designed to evaluate the
#' output files.
#'
#' If you publish data or results produced with MHCtools, please cite both of
#' the following references:
#' Roved, J. 2022. MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J., Hansson, B., Stervander, M., Hasselquist, D., & Westerdahl, H. 2022.
#' MHCtools – an R package for MHC high‐throughput sequencing data: genotyping,
#' haplotype and supertype inference, and downstream genetic analyses in non‐model
#' organisms. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13645
#'
#' @param nest_table is a table containing the sample names of parents and
#'   offspring in each nest. This table should be organized so that the
#'   individual names are in the first column (Sample_ID), and the nest number
#'   is in the second column (Nest). For each nest, the first two rows should be
#'   the parents, followed immediately by the offspring in the subsequent rows,
#'   and then followed by the next nest, and so on. It is assumed that nests are
#'   numbered consecutively beginning at 1.
#' @param seq_table seq_table is a sequence table as output by the 'dada2'
#'   pipeline, which has samples in rows and nucleotide sequence variants in
#'   columns.
#' @param path_out is a user defined path to the folder where the output files
#'   will be saved.
#' @return  A set of R lists containing for each nest the putative haplotypes,
#'   the names of sequences that could not be resolved with certainty in each
#'   parent, the names of the sequences that were incongruent in the genotypes
#'   of the nest, and the mean proportion of incongruent sequences (which is a
#'   measure of the haplotype inference success and largely influenced by the
#'   exactness of the genotyping experiment). The sequences are named in the
#'   output by an index number corresponding to their column number in the
#'   sequence table, thus identical sequences will have identical sample names
#'   in all the output files. These files can be reopened in R e.g. using the
#'   readRDS() function in the base package.
#' @seealso \code{\link{GetHpltTable}}; \code{\link{GetHpltStats}}; for more
#'   information about 'dada2' visit <https://benjjneb.github.io/dada2/>
#' @examples
#' nest_table <- nest_table
#' seq_table <- sequence_table
#' path_out <- tempdir()
#' HpltFind(nest_table, seq_table, path_out)
#' @export

HpltFind <- function(nest_table, seq_table, path_out) {

  # The dada2 sequence table does not use sequences names, but identifies
  # sequence variants by their nuceotide sequence. Here I create a vector for
  # naming the sequences by their column number in the seq_table

  seq_names <- vector("character", length=length(colnames(seq_table)))

  seq_names <- paste("Sequence_", seq(1:length(colnames(seq_table))), sep = "")

  colnames(seq_table) <- seq_names

  # Define the number of nests in the data set

  No_nests <- max(nest_table$Nest)

  for (i in 1:No_nests) {

    Nest_samples <- factor()

    Nest_samples <- nest_table[nest_table[,"Nest"]==i,"Sample_ID"]

    No_Chicks <- length(Nest_samples)-2


    # Put the names of the parents into a vector called Parent_names

    Parent_names <- c(paste(Nest_samples[1]), paste(Nest_samples[2]))

    ##  Note: In the nest table, I list the female followed by the male parent,
    ##  so Parent_names[1] is always the female and Parent_names[2] is always
    ##  the male


    # Put the names of the chicks into a vector called Chick_names

    Chick_names <- 0

    for (j in 3:length(Nest_samples)) {

      Chick_names[(j-2)] <- paste(Nest_samples[j])

    }


    # Create two lists that will contain vectors for each chick with the names
    # of the sequences that are found in each parent

    CP1 <- vector("list", No_Chicks)
    CP2 <- vector("list", No_Chicks)

    # Create a list that will contain vectors for each chick with the names of
    # the sequences that are not found in either parent

    Cinc <- vector("list", No_Chicks)


    for (j in 1:No_Chicks) {

      # Create vectors CP1[[j]] and CP2[[j]] which for each chick will contain
      # the names of the sequences that are found in each parent

      CP1[[j]] <- factor()
      CP2[[j]] <- factor()

      # Create a vector Cinc[[j]] which for each chick will contain the names of
      # the sequences that are not found in either parent

      Cinc[[j]] <- factor()

      # Fetch row numbers for the sequences in chick j in the occurrence matrix

      z <- which(seq_table[Chick_names[j],] > 0)

      # Enter the names of the sequences in chick j into a vector called Cseqs

      Cseqs <- seq_names[z]

      for (Seq in Cseqs) {

        # If Seq has > 0 reads in parent 1, append Seq to CP1[[j]]

        if (seq_table[Parent_names[1],Seq] > 0) {

          CP1[[j]] <- append(CP1[[j]], Seq)

        }

        # If Seq has > 0 reads in parent 2, append Seq to CP2[[j]]

        if (seq_table[Parent_names[2],Seq] > 0) {

          CP2[[j]] <- append(CP2[[j]], Seq)

        }

        # If Seq has 0 reads in both parents, append Seq to Cinc[[j]]

        if (seq_table[Parent_names[1],Seq] == 0 & seq_table[Parent_names[2],Seq] == 0) {

          Cinc[[j]] <- append(Cinc[[j]], Seq)

        }

      }

    }


    # Create vectors P1A , P1B, P2A, P2B, which will hold the names of the
    # sequences in the putative haplotypes

    # Note: P1 is the female and P2 the male parent

    P1A <- factor()
    P1B <- factor()
    P2A <- factor()
    P2B <- factor()

    # Create vectors P1Ainc , P1Binc, P2Ainc, P2Binc, which will hold the names
    # of sequences that are incongruent in the putative haplotypes

    P1Ainc <- factor()
    P1Binc <- factor()
    P2Ainc <- factor()
    P2Binc <- factor()

    # Assign the sequences from the first chick to the A haplotype of the parent
    # in which they were matched

    P1A <- CP1[[1]]
    P2A <- CP2[[1]]


    if (No_Chicks > 1) {

      for (j in 2:No_Chicks) {

        # Logical vectors with the congruence between CP1[[j]] and P1A

        x <- CP1[[j]] %in% P1A
        y <- P1A %in% CP1[[j]]

        # If proportion of matches in comparisons between CP1[[j]] and P1A > 0.8,
        # do

        if ((length(which(x==TRUE))+length(which(y==TRUE)))/(length(x)+length(y)) > 0.8) {

          for (Seq in CP1[[j]]) {

            # If Seq is not present in P1A, do

            if (Seq %in% P1A == FALSE) {

              # Append Seq to P1A
              P1A <- append(P1A, Seq)

              # Append Seq to P1Ainc
              P1Ainc <- append(P1Ainc, Seq)

            }

          }

        } else {

          # If no sequences were yet assigned to the P1B haplotype, it gets the
          # sequences from CP1[[j]]

          if (length(P1B)==0) {

            P1B <- CP1[[j]]

          } else {

            for (Seq in CP1[[j]]) {

              # If Seq is not present in P1B, do

              if (Seq %in% P1B == FALSE) {

                # Append Seq to P1B
                P1B <- append(P1B, Seq)

                # Append Seq to P1Binc
                P1Binc <- append(P1Binc, Seq)

              }

            }

          }

        }

      }



      for (j in 2:No_Chicks) {

        # Logical vectors with the congruence between CP2[[j]] and P2A

        x <- CP2[[j]] %in% P2A
        y <- P2A %in% CP2[[j]]

        # If proportion of matches in comparisons between CP2[[j]] and P2A > 0.8,
        # do

        if ((length(which(x==TRUE))+length(which(y==TRUE)))/(length(x)+length(y)) > 0.8) {

          for (Seq in CP2[[j]]) {

            # If Seq is not present in P2A, do

            if (Seq %in% P2A == FALSE) {

              # Append Seq to P2A
              P2A <- append(P2A, Seq)

              # Append Seq to P2Ainc
              P2Ainc <- append(P2Ainc, Seq)

            }

          }

        } else {

          # If no sequences were yet assigned to the P2B haplotype, it gets the
          # sequences from CP2[[j]]

          if (length(P2B)==0) {

            P2B <- CP2[[j]]

          } else {

            for (Seq in CP2[[j]]) {

              # If Seq is not present in P2B, do

              if (Seq %in% P2B == FALSE) {

                # Append Seq to P2B
                P2B <- append(P2B, Seq)

                # Append Seq to P2Binc
                P2Binc <- append(P2Binc, Seq)

              }

            }

          }

        }

      }

    }

    ### Now the proportions of sequences that are incongruent in the chicks can
    ### be calculated (i.e. sequences that were found in a chick, but were
    ### absent either from both parents or from other chicks that shared the
    ### same haplotype)

    if(length(P1A) > 0) PrInc_P1A <- length(P1Ainc)/length(P1A) else PrInc_P1A <- NA
    if(length(P1B) > 0) PrInc_P1B <- length(P1Binc)/length(P1B) else PrInc_P1B <- NA
    if(length(P2A) > 0) PrInc_P2A <- length(P2Ainc)/length(P2A) else PrInc_P2A <- NA
    if(length(P2B) > 0) PrInc_P2B <- length(P2Binc)/length(P2B) else PrInc_P2B <- NA

    PrInc_C <- vector("numeric", No_Chicks)

    for (j in 1:No_Chicks) {

      if(length(c(CP1[[j]],CP2[[j]],Cinc[[j]])) > 0) PrInc_C[j] <- length(Cinc[[j]])/(length(c(CP1[[j]],CP2[[j]],Cinc[[j]]))) else PrInc_C[j] <- NA

    }


    ### Fetch the names of the sequences in parent 1 and 2 and match them
    ### against P1A and P1B, P2A and P2B, respectively, and add sequences that
    ### didn’t match against either the A or the B haplotype to vectors in the
    ### list Pinc


    # Create a list that will contain vectors for the (names of) sequences that
    # are found in each parent

    Pseqs <- vector("list", 2)


    # Create a list that will contain vectors for the (names of) sequences that
    # are found in either parent but not in any of the chicks

    Pinc <- vector("list", 2)

    # Create a list that will contain vectors for the (names of) sequences that
    # are found in both haplotypes of a parent. These sequences cannot be
    # assigned to either haplotype with 100% certainty, and will therefore be
    # listed as unresolved.

    Punrs <- vector("list", 2)


    # Create two lists of each two vectors that allow looping over the haplotype
    # names

    A_hplts <- vector("list", 2)
    A_hplts[[1]] <- P1A
    A_hplts[[2]] <- P2A
    B_hplts <- vector("list", 2)
    B_hplts[[1]] <- P1B
    B_hplts[[2]] <- P2B


    for (j in 1:2) {

      Pseqs[[j]] <- factor()

      Pinc[[j]] <- factor()

      Punrs[[j]] <- factor()

      # Fetch row numbers for the sequences in parent j in the occurrence matrix

      z <- which(seq_table[Parent_names[j],] > 0)

      # Put the names of the seqs in parent j into Pseqs[[j]]

      Pseqs[[j]] <- seq_names[z]

      for (Seq in Pseqs[[j]]) {

        # If Seq is not present in either the A or B haplotype in parent j, do

        if (Seq %in% A_hplts[[j]] == FALSE & Seq %in% B_hplts[[j]] == FALSE) {

          # Append Seq to Pinc[[j]]

          Pinc[[j]] <- append(Pinc[[j]], Seq)

        }

        # If Seq is present in both the A and B haplotypes in parent j, do

        if (Seq %in% A_hplts[[j]] == TRUE & Seq %in% B_hplts[[j]] == TRUE) {

          # Append Seq to Punrs[[j]]

          Punrs[[j]] <- append(Punrs[[j]], Seq)

        }

      }

    }


    ### Now the proportions of sequences that are incongruent in the parents can
    ### be calculated, i.e. sequences that are present in the parents but not
    ### found in any chicks

    PrInc_P <- vector("numeric", 2)

    for (j in 1:2) {

      if(length(c(P1A,P1B,Pinc[[j]])) > 0) PrInc_P[j] <- length(Pinc[[j]])/(length(c(P1A,P1B,Pinc[[j]]))) else PrInc_P[j] <- NA

    }



    ### Finally, output a .Rds file with the mean proportion of incongruent
    ### sequences, the putative haplotypes, and lists of the unresolved
    ### sequences and incongruent sequences

    PrInc <- mean(c(PrInc_P1A, PrInc_P1B, PrInc_P2A, PrInc_P2B, PrInc_C[seq(1:No_Chicks)], PrInc_P[1], PrInc_P[2]),na.rm=TRUE)

    Put_haplotypes <- list(P1A, P1B, P2A, P2B)

    names(Put_haplotypes) <- c(paste(Parent_names[1], "A", sep=""), paste(Parent_names[1], "B", sep=""), paste(Parent_names[2], "A", sep=""), paste(Parent_names[2], "B", sep=""))

    names(Punrs) <- Parent_names

    names(Cinc) <- Chick_names

    names(Pinc) <- Parent_names

    Inc_Seqs <- list(P1Ainc, P1Binc, P2Ainc, P2Binc, Cinc, Pinc)

    names(Inc_Seqs) <- c(paste(Parent_names[1], "A", sep=""), paste(Parent_names[1], "B", sep=""), paste(Parent_names[2], "A", sep=""), paste(Parent_names[2], "B", sep=""),"Chicks","Parents")

    Output <- list(Mean_prop_incongr_seqs=c(PrInc), Putative_haplotypes=c(Put_haplotypes), Unresolved_seqs=c(Punrs), Incongruent_seqs=c(Inc_Seqs))

    saveRDS(Output, file=paste(path_out, "/Haplotypes_nest", i, "_", c(format(Sys.Date(),"%Y%m%d")), ".Rds", sep=""))

  }

}
