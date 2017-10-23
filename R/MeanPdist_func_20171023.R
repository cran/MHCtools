#' MeanPdist() function
#'
#' \code{\link{MeanPdist}} calculates the mean p-distance from pairwise
#' comparisons of the sequences in each sample in a 'dada2' sequence table.
#' Note: Sequences are required to be of equal length.
#'
#' @param seq_table seq_table is a sequence table as output by the 'dada2'
#'   pipeline, which has samples in rows and nucleotide sequence variants in
#'   columns.
#' @param codon_pos optional, a vector of codon positions to include in
#'   p-distance calculations, if this argument is omitted, p-distance
#'   calculations are made using all codons.
#' @return A table with the mean p-distance for each sample.
#' @seealso For more information about 'dada2'visit
#'   <https://benjjneb.github.io/dada2>
#' @examples
#' seq_table <- sequence_table_fas
#' codon_pos <- c(1,2,3,4,5,6,7,8)
#' MeanPdist(seq_table, codon_pos)
#' @importFrom "utils" "combn"
#' @export

MeanPdist <- function(seq_table,codon_pos=NULL) {

  # Extract the sample names to a new vector

  sample_names <- rownames(seq_table)

  # Create a vector mean_Pdist

  mean_Pdist <- vector("numeric", length=length(sample_names))

  # Create a vector seq_list containing all the sequences in seq_table,  using
  # the strsplit() function to split the nucleotides in each sequence into
  # separate elements

  seqs <- colnames(seq_table)

  seq_list <- list()

  for (i in 1:length(seqs)) {

    seq_list[i] <- strsplit(seqs[i],"")

  }

  # for loop over all the samples in the data set

  for (i in 1:length(sample_names)) {

    # Fetch column numbers for the sequences in sample i in the dada2 sequence
    # table

    z <- which(seq_table[i,] > 0)

    # Create a vector lengths with the lengths of the sequences in sample i

    lengths <- vector()

    for(j in 1:length(z))  {

      lengths[j] <- length(seq_list[[j]])

    }

    # Throw a warning if sequences in sample i are of different lengths

    if(max(lengths) != min(lengths)) {

      stop("Pairwise comparisons not meaningful for sequences of different length")

    } else {

      # Create a vector pd

      pd <- vector("numeric", length=length(z))

      # Generate a list of all pairwise combinations of the elements in z

      pwc <- combn(sort(z),2,simplify=T)

      # For each combination in pwc

      for(j in 1:length(pwc[1,]))  {

        if(is.null(codon_pos)) {

          # Calculate the Pdist in each pairwise comparison of the sequences
          # in seq_list, using the numbers from pwc as indices to extract the
          # sequences from the list

          pd[j] <- length(which(seq_list[[pwc[1,j]]] != seq_list[[pwc[2,j]]]))/length(seq_list[[pwc[1,j]]])

        } else {

          # Throw a warning if selected codons exceed sequence length

          if(max(codon_pos) > max(lengths)) {

            stop("Selected codons exceed sequence length")

          } else {

            # Calculate the Pdist in each pairwise comparison of the sequences
            # in seq_list, using the numbers from pwc as indices to extract the
            # sequences from the list, and using the vector codon_pos to define
            # which codons to compare

            pd[j] <- length(which(seq_list[[pwc[1,j]]][codon_pos] != seq_list[[pwc[2,j]]][codon_pos]))/length(seq_list[[pwc[1,j]]][codon_pos])

          }

        }

      }

    }

    # Calculate the mean of the pd vector and add this to the mean_Pdist vector

    mean_Pdist[i] <- mean(pd)

  }

  # Output a table with sample names and mean Pdistances

  tab <- as.data.frame(mean_Pdist)

  rownames(tab) <- sample_names

  tab

}
