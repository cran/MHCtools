#' CalcPdist() function
#'
#' \code{\link{CalcPdist}} calculates p-distances from pairwise sequence
#' comparisons and mean p-distances for each sample in a 'dada2' sequence
#' table.
#'
#' @param seq_file seq_file is a sequence table as output by the 'dada2'
#'   pipeline, which has samples in rows and nucleotide sequence variants in
#'   columns. Optionally, a fasta file can be supplied as input in the format
#'   rendered by e.g. read.fasta() from the package 'seqinr'.
#' @param path_out is a user defined path to the folder where the output files
#'   will be saved.
#' @param aa_pdist optional, a logical (TRUE/FALSE) that determines whether
#'   nucleotide sequences should be translated to amino acid sequences before
#'   p-distance calculation, default is NULL/FALSE.
#' @param codon_pos optional, a vector of codon positions to include in
#'   p-distance calculations, if this argument is omitted, p-distance
#'   calculations are made using all codons.
#' @param input_fasta optional, a logical (TRUE/FALSE) that indicates whether
#'   the input file is a fasta file (TRUE) or a dada2 sequence table
#'   (NULL/FALSE), default is NULL/FALSE.
#' @return If a sequence table is given as input file, the function returns a
#'   table with the mean p-distance for each sample. Additionally, the function
#'   produces a matrix with p-distance values of all pairwise sequence
#'   comparisons. This table is saved as a .csv file in the output path. If a
#'   fasta file is used as input, only the p-distance matrix will be produced.
#' @seealso For more information about 'dada2'visit
#'   <https://benjjneb.github.io/dada2>
#' @examples
#' seq_file <- sequence_table_fas
#' path_out <- tempdir()
#' CalcPdist(seq_file, path_out, aa_pdist=NULL, codon_pos=c(1,2,3,4,5,6,7,8), input_fasta=NULL)
#' @importFrom "utils" "combn"
#' @export

CalcPdist <- function(seq_file,path_out,aa_pdist=NULL,codon_pos=NULL,input_fasta=NULL) {

  ###### Sequence table as input file ######

  if(is.null(input_fasta) || isFALSE(input_fasta)) {

    # Extract the sample names to a new vector

    sample_names <- rownames(seq_file)

    # Create a vector mean_Pdist

    mean_Pdist <- vector("numeric", length=length(sample_names))

    # create a matrix that will contain the pairwise p-distance between all the sequences in the sequence table

    pdist_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=length(colnames(seq_file))))
    rownames(pdist_matrix) <- colnames(seq_file)
    colnames(pdist_matrix) <- colnames(seq_file)

    # Create a vector seq_list containing all the sequences in seq_file, using
    # the strsplit() function to split the nucleotides in each sequence into
    # separate elements

    seqs <- colnames(seq_file)

    seq_list <- list()

    for (i in 1:length(seqs)) {

      seq_list[i] <- strsplit(seqs[i],"")

    }

    # Create a vector with the lengths of the sequences

    lengths <- vector()

    for(j in 1:length(seq_list))  {

      lengths[j] <- length(seq_list[[j]])

    }

    # Throw a warning if sequences in the sequence table are of different lengths

    if(max(lengths) != min(lengths)) {

      stop("Pairwise comparisons not meaningful for sequences of different length")

    }

    ### Create a vector seq_list_aa containing all the sequences in seq_file translated to amino acids

    if(isTRUE(aa_pdist)) {

      seq_list_aa <- list(length=length(seq_list))

      # Define the genetic code

      Phe <- c("TTT","TTC")
      Leu <- c("TTA","TTG","CTT","CTC","CTA","CTG")
      Ile <- c("ATT","ATC","ATA")
      Met <- c("ATG")
      Val <- c("GTT","GTC","GTA","GTG")
      Ser <- c("TCT","TCC","TCA","TCG","AGT","AGC")
      Pro <- c("CCT","CCC","CCA","CCG")
      Thr <- c("ACT","ACC","ACA","ACG")
      Ala <- c("GCT","GCC","GCA","GCG")
      Tyr <- c("TAT","TAC")
      His <- c("CAT","CAC")
      Gln <- c("CAA","CAG")
      Asn <- c("AAT","AAC")
      Lys <- c("AAA","AAG")
      Asp <- c("GAT","GAC")
      Glu <- c("GAA","GAG")
      Cys <- c("TGT","TGC")
      Trp <- c("TGG")
      Arg <- c("CGT","CGC","CGA","CGG","AGA","AGG")
      Gly <- c("GGT","GGC","GGA","GGG")
      Stop <- c("TAA","TAG","TGA")

      # translate the nucleotide sequences to amino acid sequences

      if(readline(prompt="Are the sequences aligned in open reading frame 5' -> 3'? y/n: ") == "y") {

        for(i in 1:length(seq_list)) {

          # fill in NaNs in seq_list_aa
          seq_list_aa[[i]] <- NaN*seq(length(seq_list[[i]])/3)

          if((length(seq_list[[i]])/3) %% 1 == 0) {

            for(j in 1:(length(seq_list[[i]])/3)) {

              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Phe) {seq_list_aa[[i]][j] <- "F"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Leu) {seq_list_aa[[i]][j] <- "L"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Ile) {seq_list_aa[[i]][j] <- "I"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Met) {seq_list_aa[[i]][j] <- "M"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Val) {seq_list_aa[[i]][j] <- "V"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Ser) {seq_list_aa[[i]][j] <- "S"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Pro) {seq_list_aa[[i]][j] <- "P"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Thr) {seq_list_aa[[i]][j] <- "T"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Ala) {seq_list_aa[[i]][j] <- "A"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Tyr) {seq_list_aa[[i]][j] <- "Y"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% His) {seq_list_aa[[i]][j] <- "H"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Gln) {seq_list_aa[[i]][j] <- "Q"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Asn) {seq_list_aa[[i]][j] <- "N"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Lys) {seq_list_aa[[i]][j] <- "K"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Asp) {seq_list_aa[[i]][j] <- "D"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Glu) {seq_list_aa[[i]][j] <- "E"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Cys) {seq_list_aa[[i]][j] <- "C"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Trp) {seq_list_aa[[i]][j] <- "W"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Arg) {seq_list_aa[[i]][j] <- "R"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Gly) {seq_list_aa[[i]][j] <- "G"}
              if(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="") %in% Stop) {seq_list_aa[[i]][j] <- "X"; warning("Stop codon detected in sequence.")}

            }

          } else {

            stop("Sequences should be aligned in open reading frame 5' -> 3'")

          }

        }

      } else {

        stop("Sequences should be aligned in open reading frame 5' -> 3'")

      }

    }

    ### calculate pairwise p-distances between all sequences in the
    ### sequence table, and save the values in the matrix
    ### (upper right matrix, values rounded to five digits)

    for (i in 1:(length(colnames(seq_file))-1)) {

      for (j in (i+1):length(colnames(seq_file))) {

        if(is.null(aa_pdist) || isFALSE(aa_pdist)) {

          # nucleotide p-distance calculation

          if(is.null(codon_pos)) {

            # Calculate the nucleotide Pdist in each pairwise comparison of the sequences
            # in seq_list

            pdist_matrix[i,j] <- round(length(which(seq_list[[i]] != seq_list[[j]]))/length(seq_list[[i]]), digits=5)

          } else {

            # Throw a warning if selected codons exceed sequence length

            if(max(codon_pos) > max(lengths)) {

              stop("Selected codons exceed sequence length")

            } else {

              # Calculate the nucleotide Pdist in each pairwise comparison of the sequences
              # in seq_list, using the vector codon_pos to define
              # which codons to compare

              pdist_matrix[i,j] <- round(length(which(seq_list[[i]][codon_pos] != seq_list[[j]][codon_pos]))/length(seq_list[[i]][codon_pos]), digits=5)

            }

          }

        }

        if(isTRUE(aa_pdist)) {

          # amino acid p-distance calculation

          if(is.null(codon_pos)) {

            # Calculate the amino acid Pdist in each pairwise comparison of the sequences
            # in seq_list_aa

            pdist_matrix[i,j] <- round(length(which(seq_list_aa[[i]] != seq_list_aa[[j]]))/length(seq_list_aa[[i]]), digits=5)

          } else {

            # Throw a warning if selected codons exceed amino acid sequence length

            if(max(codon_pos) > length(seq_list_aa[[i]])) {

              stop("Selected codons exceed amino acid sequence length")

            } else {

              # Calculate the amino acid Pdist in each pairwise comparison of the sequences
              # in seq_list_aa, using the vector codon_pos to define
              # which codons to compare

              pdist_matrix[i,j] <- round(length(which(seq_list_aa[[i]][codon_pos] != seq_list_aa[[j]][codon_pos]))/length(seq_list_aa[[i]][codon_pos]), digits=5)

            }

          }

        }

      }

    }

    ### Calculate mean p-distances for each sample in the data set

    for (i in 1:length(sample_names)) {

      # Fetch column numbers for the sequences in sample i in the dada2 sequence
      # table

      z <- which(seq_file[i,] > 0)

      # Create a vector pd

      pd <- vector("numeric", length=length(z))

      # Generate a list of all pairwise combinations of the elements in z

      pwc <- combn(sort(z),2,simplify=T)

      # For each combination in pwc

      for(j in 1:length(pwc[1,]))  {

        # Extract the Pdist in each pairwise comparison of the sequences
        # in seq_list from the pdist_matrix, using the numbers from pwc as
        # indices to extract the values from the matrix

        pd[j] <- pdist_matrix[pwc[1,j],pwc[2,j]]

      }

      # Calculate the mean of the pd vector and add this to the mean_Pdist vector

      mean_Pdist[i] <- mean(pd)

    }

    # Output a table with sample names and mean Pdistances

    tab <- as.data.frame(mean_Pdist)

    rownames(tab) <- sample_names

    print(tab)

  }

  ###### Fasta as input file ######

  if(isTRUE(input_fasta)) {

    # fasta files are accepted in the list format rendered by
    # the read.fasta() function from the package 'seqinr'

    # create a matrix that will contain the pairwise p-distance between all the sequences in the fasta file

    pdist_matrix <- as.data.frame(matrix(nrow=length(seq_file),ncol=length(seq_file)))
    rownames(pdist_matrix) <- names(seq_file)
    colnames(pdist_matrix) <- names(seq_file)

    # Create a vector with the lengths of the sequences

    lengths <- vector()

    for(j in 1:length(seq_file))  {

      lengths[j] <- length(seq_file[[j]])

    }

    # Throw a warning if sequences in the fasta file are of different lengths

    if(max(lengths) != min(lengths)) {

      stop("Pairwise comparisons not meaningful for sequences of different length")

    }

    ### Create a vector seq_list_aa containing all the sequences in the fasta file translated to amino acids

    if(isTRUE(aa_pdist)) {

      seq_list_aa <- list(length=length(seq_file))

      # Define the genetic code

      Phe <- c("TTT","TTC")
      Leu <- c("TTA","TTG","CTT","CTC","CTA","CTG")
      Ile <- c("ATT","ATC","ATA")
      Met <- c("ATG")
      Val <- c("GTT","GTC","GTA","GTG")
      Ser <- c("TCT","TCC","TCA","TCG","AGT","AGC")
      Pro <- c("CCT","CCC","CCA","CCG")
      Thr <- c("ACT","ACC","ACA","ACG")
      Ala <- c("GCT","GCC","GCA","GCG")
      Tyr <- c("TAT","TAC")
      His <- c("CAT","CAC")
      Gln <- c("CAA","CAG")
      Asn <- c("AAT","AAC")
      Lys <- c("AAA","AAG")
      Asp <- c("GAT","GAC")
      Glu <- c("GAA","GAG")
      Cys <- c("TGT","TGC")
      Trp <- c("TGG")
      Arg <- c("CGT","CGC","CGA","CGG","AGA","AGG")
      Gly <- c("GGT","GGC","GGA","GGG")
      Stop <- c("TAA","TAG","TGA")

      # translate the nucleotide sequences to amino acid sequences

      if(readline(prompt="Are the sequences aligned in open reading frame 5' -> 3'? y/n: ") == "y") {

        for(i in 1:length(seq_file)) {

          # fill in NaNs in seq_list_aa
          seq_list_aa[[i]] <- NaN*seq(length(seq_file[[i]])/3)

          if((length(seq_file[[i]])/3) %% 1 == 0) {

            for(j in 1:(length(seq_file[[i]])/3)) {

              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Phe) {seq_list_aa[[i]][j] <- "F"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Leu) {seq_list_aa[[i]][j] <- "L"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Ile) {seq_list_aa[[i]][j] <- "I"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Met) {seq_list_aa[[i]][j] <- "M"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Val) {seq_list_aa[[i]][j] <- "V"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Ser) {seq_list_aa[[i]][j] <- "S"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Pro) {seq_list_aa[[i]][j] <- "P"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Thr) {seq_list_aa[[i]][j] <- "T"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Ala) {seq_list_aa[[i]][j] <- "A"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Tyr) {seq_list_aa[[i]][j] <- "Y"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% His) {seq_list_aa[[i]][j] <- "H"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Gln) {seq_list_aa[[i]][j] <- "Q"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Asn) {seq_list_aa[[i]][j] <- "N"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Lys) {seq_list_aa[[i]][j] <- "K"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Asp) {seq_list_aa[[i]][j] <- "D"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Glu) {seq_list_aa[[i]][j] <- "E"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Cys) {seq_list_aa[[i]][j] <- "C"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Trp) {seq_list_aa[[i]][j] <- "W"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Arg) {seq_list_aa[[i]][j] <- "R"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Gly) {seq_list_aa[[i]][j] <- "G"}
              if(toupper(paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Stop) {seq_list_aa[[i]][j] <- "X"; warning("Stop codon detected in sequence.")}

            }

          } else {

            stop("Sequences should be aligned in open reading frame 5' -> 3'")

          }

        }

      } else {

        stop("Sequences should be aligned in open reading frame 5' -> 3'")

      }

    }

    ### calculate pairwise p-distances between all sequences in the
    ### sequence table, and save the values in the matrix
    ### (upper right matrix, values rounded to five digits)

    for (i in 1:(length(seq_file)-1)) {

      for (j in (i+1):length(seq_file)) {

        if(is.null(aa_pdist) || isFALSE(aa_pdist)) {

          # nucleotide p-distance calculation

          if(is.null(codon_pos)) {

            # Calculate the nucleotide Pdist in each pairwise comparison of the sequences
            # in seq_file

            pdist_matrix[i,j] <- round(length(which(seq_file[[i]] != seq_file[[j]]))/length(seq_file[[i]]), digits=5)

          } else {

            # Throw a warning if selected codons exceed sequence length

            if(max(codon_pos) > max(lengths)) {

              stop("Selected codons exceed sequence length")

            } else {

              # Calculate the nucleotide Pdist in each pairwise comparison of the sequences
              # in seq_file, using the vector codon_pos to define
              # which codons to compare

              pdist_matrix[i,j] <- round(length(which(seq_file[[i]][codon_pos] != seq_file[[j]][codon_pos]))/length(seq_file[[i]][codon_pos]), digits=5)

            }

          }

        }

        if(isTRUE(aa_pdist)) {

          # amino acid p-distance calculation

          if(is.null(codon_pos)) {

            # Calculate the amino acid Pdist in each pairwise comparison of the sequences
            # in seq_list_aa

            pdist_matrix[i,j] <- round(length(which(seq_list_aa[[i]] != seq_list_aa[[j]]))/length(seq_list_aa[[i]]), digits=5)

          } else {

            # Throw a warning if selected codons exceed amino acid sequence length

            if(max(codon_pos) > length(seq_list_aa[[i]])) {

              stop("Selected codons exceed amino acid sequence length")

            } else {

              # Calculate the amino acid Pdist in each pairwise comparison of the sequences
              # in seq_list_aa, using the vector codon_pos to define
              # which codons to compare

              pdist_matrix[i,j] <- round(length(which(seq_list_aa[[i]][codon_pos] != seq_list_aa[[j]][codon_pos]))/length(seq_list_aa[[i]][codon_pos]), digits=5)

            }

          }

        }

      }

    }

  }

  # Export the p-distance matrix as a .csv file
  write.csv(pdist_matrix,file=paste(path_out,"/Pdist_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))

}
