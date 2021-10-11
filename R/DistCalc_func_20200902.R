#' DistCalc() function
#'
#' \code{\link{DistCalc}} calculates Grantham distances, Sandberg distances,
#' or p-distances from pairwise comparisons of aligned sequences.
#'
#' The DistCalc() function takes a fasta file or a 'dada2'-style sequence
#' occurrence table (with aligned sequences as column names and samples in
#' rows) as input and produces a matrix with pairwise distances for all
#' sequences in the data set. If calculation of Sandberg distances is specified,
#' the function additionally outputs five tables with physico-chemical
#' z-descriptor values (based on Sandberg et al. 1998) for each amino acid
#' position in all sequences in the data set. These tables may be useful for
#' further downstream analyses, such as estimation of MHC supertypes. If a
#' sequence occurrence table is provided as input, the DistCalc() function
#' furthermore produces a table with the mean distances from all pairwise
#' comparisons of the sequences in each sample in the data set.
#'
#' Grantham distances and Sandberg distances are calculated as described in
#' Pierini & Lenz 2018. The Grantham distances produced by DistCalc() are
#' simply the mean Grantham distances (Grantham 1974) between all amino acid
#' codons in sequence pairs. When calculating Sandberg distances, DistCalc()
#' first computes Euclidian distances between all amino acid pairs based on
#' the five physico-chemical z-descriptors defined in Sandberg et al. 1998.
#' Sandberg distances are then calculated as the mean Euclidian distances
#' between all amino acid codons in sequence pairs. P-distances calculated
#' by DistCalc() are simply the proportion of varying codons between pairs
#' of sequences.
#'
#' The DistCalc() function includes an option for the user to specify which
#' codons to compare, which is useful e.g. if conducting the analysis only
#' on codon positions involved in specific functions, such as peptide binding
#' of an MHC molecule. It also accepts calculating amino acid distances
#' directly from protein-coding DNA sequences using the standard genetic
#' code.
#'
#' The DistCalc() function accepts the following characters in the sequences:
#' Nucleotide sequences: A,T,G,C
#' Amino acid sequences: A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V
#'
#' It accepts gaps defined by '-'. Nucleotide triplets containing gaps are
#' translated to 'X', if amino acid distances are calculated directly from DNA
#' nucleotide sequences. Please note that '-' or 'X' are treated as unique
#' characters in p-distance calculations. The function will not accept 'X' or
#' gaps in Grantham or Sandberg distance calculations. If you wish to exclude
#' codons with 'X' or gaps from distance calculations, please use the
#' codon_pos option to specify which codons to compare.
#'
#' If you publish data produced with MHCtools, please cite:
#' Roved, J. 2020. MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J., Hansson, B., Stervander, M., Hasselquist, D., Westerdahl, H. 2020.
#' Non-random association of MHC-I alleles in favor of high diversity haplotypes
#' in wild songbirds revealed by computer-assisted MHC haplotype inference using
#' the R package MHCtools. bioRxiv.
#'
#' If you calculated Grantham or Sandberg distances, please additionally cite:
#' Pierini, F., Lenz, T.L. 2018. Divergent allele advantage at human MHC genes:
#' Signatures of past and ongoing selection. Mol. Biol. Evol. 35, 2145–2158.
#'
#' ...and either of the following references:
#' Grantham R. 1974. Amino acid difference formula to help explain protein
#' evolution. Science 185:862–864.
#' Sandberg M, Eriksson L, Jonsson J, Sjostrom M, Wold S. 1998. New chemical
#' descriptors relevant for the design of biologically active peptides. A
#' multivariate characterization of 87 amino acids. JMed Chem. 41(14):2481–2491.
#'
#' @param seq_file is a sequence occurrence table as output by the 'dada2'
#'   pipeline, which has samples in rows and nucleotide sequence variants in
#'   columns. Optionally, a fasta file can be supplied as input in the format
#'   rendered by read.fasta() from the package 'seqinr'.
#' @param path_out is a user defined path to the folder where the output files
#'   will be saved.
#' @param input_fasta optional, a logical (TRUE/FALSE) that indicates whether
#'   the input file is a fasta file (TRUE) or a 'dada2'-style sequence table
#'   (NULL/FALSE). The default is NULL/FALSE.
#' @param input_seq defines the type of sequences in seq_file. It may take the
#'   values 'nucl' or 'aa'.
#' @param aa_dist is optional, a logical (TRUE/FALSE) that determines whether
#'   nucleotide sequences should be translated to amino acid sequences before
#'   distance calculation, default is NULL/FALSE. Note that aa_dist must be set
#'   to TRUE, if Grantham or Sandberg distances are calculated from an alignment
#'   of nucleotide sequences.
#' @param codon_pos is optional, a vector of comma separated integers specifying
#'   which codon positions to include in distance calculations. If omitted,
#'   distance calculations are made using all codons.
#' @param dist_type is used to specify which kind of distances that are
#'   calculated. It takes the values 'G' for Grantham distances, 'S'  for
#'   Sandberg distances, or 'P' for p-distances. The argument is optional with
#'   'G' as default setting.
#' @return The function returns a matrix with distances from all pairwise sequence
#'   comparisons, where n is the number of sequences. If a sequence occurrence
#'   table is given as input file, the function additionally returns a table with
#'   the mean distance for each sample in the data set. If a sequence occurrence
#'   table is given as input file, the sequences are named in the output matrix by
#'   an index number that corresponds to their column number in the input file. If
#'   calculation of Sandberg distances is specified, the function additionally
#'   outputs five tables with physico-chemical z-descriptor values for each amino
#'   acid position in all sequences in the data set. All output tables are saved
#'   as .csv files in the output path.
#' @seealso For more information about 'dada2', visit
#'   <https://benjjneb.github.io/dada2/>
#' @examples
#' seq_file <- sequence_table_fas
#' path_out <- tempdir()
#' DistCalc(seq_file, path_out, input_fasta=NULL, input_seq="nucl", aa_dist=NULL,
#' codon_pos=c(1,2,3,4,5,6,7,8), dist_type="P")
#' @importFrom "utils" "combn"
#' @importFrom "stats" "dist"
#' @export

DistCalc <- function(seq_file, path_out, input_fasta=NULL, input_seq="aa", aa_dist=NULL, codon_pos=NULL, dist_type="G") {

  ##########################################
  ###### Sequence table as input file ######
  ##########################################

  if(is.null(input_fasta) || isFALSE(input_fasta)) {

    # The dada2 sequence table does not use sequences names, but identifies
    # sequence variants by their nuceotide sequence. Here I create a vector for
    # naming the sequences by their column number in the seq_table

    seq_names <- vector("character", length=length(colnames(seq_file)))

    seq_names <- paste("Sequence_", seq(1:length(colnames(seq_file))), sep = "")

    # Extract the sample names to a new vector

    sample_names <- rownames(seq_file)

    # Create a vector mean_dist that will get the mean distances per individual

    mean_dist <- vector("numeric", length=length(sample_names))

    # create a matrix that will contain the pairwise distances between all the sequences in the sequence table

    dist_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=length(colnames(seq_file))))
    rownames(dist_matrix) <- seq_names
    colnames(dist_matrix) <- seq_names

    # Create a list seq_list containing all the sequences in seq_file, using
    # the strsplit() function to split the nucleotides in each sequence into
    # separate elements

    seqs <- colnames(seq_file)

    seq_list <- list()

    for (i in 1:length(seqs)) {

      seq_list[i] <- strsplit(seqs[i],"")

    }

    # Throw a warning if sequences in the sequence table are of different lengths

    if(max(lengths(seq_list)) != min(lengths(seq_list))) {

      stop("Pairwise comparisons not meaningful for sequences of different length.")

    }

    # Throw a warning if sequences in the sequence table contain non-standard characters

    if(input_seq == "aa") {

      for(i in 1:length(seq_list)) {

        for(j in 1:length(unique(seq_list[[i]]))) {

          if(!(unique(seq_list[[i]])[j] %in% c("-","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"))) {

            stop("Sequences contain non-standard characters. Accepted characters are -ARNDCQEGHILKMFPSTWYV.")

          }

        }

      }

    }

    if(input_seq == "nucl") {

      for(i in 1:length(seq_list)) {

        for(j in 1:length(unique(seq_list[[i]]))) {

          if(!(toupper(unique(seq_list[[i]])[j]) %in% c("-","A","C","G","T"))) {

            stop("Sequences contain non-standard characters. Accepted characters are -ATGC.")

          }

        }

      }

    }

    ### Create a vector seq_list_aa containing all the sequences in seq_file translated to amino acids

    if(isTRUE(aa_dist)) {

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
              if(list(c("-")) %in% paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)])) {seq_list_aa[[i]][j] <- "X"; warning("Gap detected in sequence.")}

            }

          } else {

            stop("Sequences should be aligned in open reading frame 5' -> 3'.")

          }

        }

      } else {

        stop("Sequences should be aligned in open reading frame 5' -> 3'.")

      }

    }

    ### calculate pairwise distances between all sequences in the
    ### sequence table, and save the values in the dist_matrix
    ### (upper right matrix, values rounded to five digits)

    ##############################
    ### p-distance calculation ###
    ##############################

    if(dist_type == "P") {

      ## For input files with nucleotide sequences ##

      if(input_seq == "nucl" && (is.null(aa_dist) || isFALSE(aa_dist))) {

        # nucleotide p-distance calculation

        if(is.null(codon_pos)) {

          # Calculate the nucleotide Pdist in each pairwise comparison of the sequences
          # in seq_list

          for (i in 1:(length(colnames(seq_file))-1)) {

            for (j in (i+1):length(colnames(seq_file))) {

              dist_matrix[i,j] <- round(length(which(seq_list[[i]] != seq_list[[j]]))/length(seq_list[[i]]), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed sequence length

          if(max(codon_pos) > max(lengths(seq_list))) {

            stop("Selected codons exceed sequence length.")

          } else {

            # Calculate the nucleotide Pdist in each pairwise comparison of the sequences
            # in seq_list, using the vector codon_pos to define which codons to compare

            for (i in 1:(length(colnames(seq_file))-1)) {

              for (j in (i+1):length(colnames(seq_file))) {

                dist_matrix[i,j] <- round(length(which(seq_list[[i]][codon_pos] != seq_list[[j]][codon_pos]))/length(seq_list[[i]][codon_pos]), digits=5)

              }

            }

          }

        }

      }

      if(input_seq == "nucl" && isTRUE(aa_dist)) {

        # amino acid p-distance calculation

        if(is.null(codon_pos)) {

          # Calculate the amino acid Pdist in each pairwise comparison of the sequences
          # in seq_list_aa

          for (i in 1:(length(colnames(seq_file))-1)) {

            for (j in (i+1):length(colnames(seq_file))) {

              dist_matrix[i,j] <- round(length(which(seq_list_aa[[i]] != seq_list_aa[[j]]))/length(seq_list_aa[[i]]), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_list_aa))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the amino acid Pdist in each pairwise comparison of the sequences
            # in seq_list_aa, using the vector codon_pos to define  which codons to
            # compare

            for (i in 1:(length(colnames(seq_file))-1)) {

              for (j in (i+1):length(colnames(seq_file))) {

                dist_matrix[i,j] <- round(length(which(seq_list_aa[[i]][codon_pos] != seq_list_aa[[j]][codon_pos]))/length(seq_list_aa[[i]][codon_pos]), digits=5)

              }

            }

          }

        }

      }

      ## For input files with amino acid sequences ##

      if(input_seq == "aa") {

        # amino acid p-distance calculation

        if(is.null(codon_pos)) {

          # Calculate the amino acid Pdist in each pairwise comparison of the sequences
          # in seq_list

          for (i in 1:(length(colnames(seq_file))-1)) {

            for (j in (i+1):length(colnames(seq_file))) {

              dist_matrix[i,j] <- round(length(which(seq_list[[i]] != seq_list[[j]]))/length(seq_list[[i]]), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_list_aa))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the amino acid Pdist in each pairwise comparison of the sequences
            # in seq_list, using the vector codon_pos to define which codons to compare

            for (i in 1:(length(colnames(seq_file))-1)) {

              for (j in (i+1):length(colnames(seq_file))) {

                dist_matrix[i,j] <- round(length(which(seq_list[[i]][codon_pos] != seq_list[[j]][codon_pos]))/length(seq_list[[i]][codon_pos]), digits=5)

              }

            }

          }

        }

      }

    }

    #####################################
    ### Grantham distance calculation ###
    #####################################

    if(dist_type == "G") {

      # Define a Grantham distance matrix

      g_dist_matrix <- rbind.data.frame(

        c(  0, 112, 111, 126, 195,  91, 107,  60,  86,  94,  96, 106,  84, 113,  27,  99,  58, 148, 112,  64),
        c(112,   0,  86,  96, 180,  43,  54, 125,  29,  97, 102,  26,  91,  97, 103, 110,  71, 101,  77,  96),
        c(111,  86,   0,  23, 139,  46,  42,  80,  68, 149, 153,  94, 142, 158,  91,  46,  65, 174, 143, 133),
        c(126,  96,  23,   0, 154,  61,  45,  94,  81, 168, 172, 101, 160, 177, 108,  65,  85, 181, 160, 152),
        c(195, 180, 139, 154,   0, 154, 170, 159, 174, 198, 198, 202, 196, 205, 169, 112, 149, 215, 194, 192),
        c( 91,  43,  46,  61, 154,   0,  29,  87,  24, 109, 113,  53, 101, 116,  76,  68,  42, 130,  99,  96),
        c(107,  54,  42,  45, 170,  29,   0,  98,  40, 134, 138,  56, 126, 140,  93,  80,  65, 152, 122, 121),
        c( 60, 125,  80,  94, 159,  87,  98,   0,  98, 135, 138, 127, 127, 153,  42,  56,  59, 184, 147, 109),
        c( 86,  29,  68,  81, 174,  24,  40,  98,   0,  94,  99,  32,  87, 100,  77,  89,  47, 115,  83,  84),
        c( 94,  97, 149, 168, 198, 109, 134, 135,  94,   0,   5, 102,  10,  21,  95, 142,  89,  61,  33,  29),
        c( 96, 102, 153, 172, 198, 113, 138, 138,  99,   5,   0, 107,  15,  22,  98, 145,  92,  61,  36,  32),
        c(106,  26,  94, 101, 202,  53,  56, 127,  32, 102, 107,   0,  95, 102, 103, 121,  78, 110,  85,  97),
        c( 84,  91, 142, 160, 196, 101, 126, 127,  87,  10,  15,  95,   0,  28,  87, 135,  81,  67,  36,  21),
        c(113,  97, 158, 177, 205, 116, 140, 153, 100,  21,  22, 102,  28,   0, 114, 155, 103,  40,  22,  50),
        c( 27, 103,  91, 108, 169,  76,  93,  42,  77,  95,  98, 103,  87, 114,   0,  74,  38, 147, 110,  68),
        c( 99, 110,  46,  65, 112,  68,  80,  56,  89, 142, 145, 121, 135, 155,  74,   0,  58, 177, 144, 124),
        c( 58,  71,  65,  85, 149,  42,  65,  59,  47,  89,  92,  78,  81, 103,  38,  58,   0, 128,  92,  69),
        c(148, 101, 174, 181, 215, 130, 152, 184, 115,  61,  61, 110,  67,  40, 147, 177, 128,   0,  37,  88),
        c(112,  77, 143, 160, 194,  99, 122, 147,  83,  33,  36,  85,  36,  22, 110, 144,  92,  37,   0,  55),
        c( 64,  96, 133, 152, 192,  96, 121, 109,  84,  29,  32,  97,  21,  50,  68, 124,  69,  88,  55,   0)

      )

      colnames(g_dist_matrix) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
      rownames(g_dist_matrix) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

      ## For input files with nucleotide sequences ##

      # Print error message if amino acid calculation is not specified

      if(input_seq == "nucl" && (is.null(aa_dist) || isFALSE(aa_dist))) {

        stop("Grantham distances can only be calculated for amino acid sequences. Please load amino acid sequences or specify aa_dist = T.")

      }

      if(input_seq == "nucl" && isTRUE(aa_dist)) {

        # Grantham distance calculation

        if(is.null(codon_pos)) {

          # Calculate the Grantham distances in each pairwise comparison of
          # the sequences in seq_list_aa

          for (i in 1:(length(colnames(seq_file))-1)) {

            for (j in (i+1):length(colnames(seq_file))) {

              g_dist_ij <- vector("numeric", length=length(seq_list_aa[[i]]))

              for(k in 1:length(seq_list_aa[[i]])) {

                g_dist_ij[k] <- g_dist_matrix[which(rownames(g_dist_matrix)==seq_list_aa[[i]][k]),which(colnames(g_dist_matrix)==seq_list_aa[[j]][k])]

              }

              dist_matrix[i,j] <- round(mean(g_dist_ij), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_list_aa))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the Grantham distance in each pairwise
            # comparison of the sequences in seq_list_aa, using
            # the vector codon_pos to define which codons to compare

            for (i in 1:(length(colnames(seq_file))-1)) {

              for (j in (i+1):length(colnames(seq_file))) {

                g_dist_ij <- vector("numeric", length=length(codon_pos))

                for(k in 1:length(codon_pos)) {

                  g_dist_ij[k] <- g_dist_matrix[which(rownames(g_dist_matrix)==seq_list_aa[[i]][codon_pos[k]]),which(colnames(g_dist_matrix)==seq_list_aa[[j]][codon_pos[k]])]

                }

                dist_matrix[i,j] <- round(mean(g_dist_ij), digits=5)

              }

            }

          }

        }

      }

      ## For input files with amino acid sequences ##

      if(input_seq == "aa") {

        # Grantham distance calculation

        if(is.null(codon_pos)) {

          # Calculate the Grantham distances in each pairwise comparison of
          # the sequences in seq_list

          for (i in 1:(length(colnames(seq_file))-1)) {

            for (j in (i+1):length(colnames(seq_file))) {

              g_dist_ij <- vector("numeric", length=length(seq_list[[i]]))

              for(k in 1:length(seq_list[[i]])) {

                g_dist_ij[k] <- g_dist_matrix[which(rownames(g_dist_matrix)==seq_list[[i]][k]),which(colnames(g_dist_matrix)==seq_list[[j]][k])]

              }

              dist_matrix[i,j] <- round(mean(g_dist_ij), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_list))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the Grantham distance in each pairwise
            # comparison of the sequences in seq_list, using
            # the vector codon_pos to define which codons to compare

            for (i in 1:(length(colnames(seq_file))-1)) {

              for (j in (i+1):length(colnames(seq_file))) {

                g_dist_ij <- vector("numeric", length=length(codon_pos))

                for(k in 1:length(codon_pos)) {

                  g_dist_ij[k] <- g_dist_matrix[which(rownames(g_dist_matrix)==seq_list[[i]][codon_pos[k]]),which(colnames(g_dist_matrix)==seq_list[[j]][codon_pos[k]])]

                }

                dist_matrix[i,j] <- round(mean(g_dist_ij), digits=5)

              }

            }

          }

        }

      }

    }

    #####################################
    ### Sandberg distance calculation ###
    #####################################

    if(dist_type == "S") {

      # Calculate Sandberg distances as the Euclidian distances between all 20 amino
      # acids based on the five physicochemical z-descriptors described in Sandberg
      # et al. (1998)

      # Define a matrix with z-descriptor values for all 20 amino acids

      z_desc_matrix <- rbind.data.frame(

        c(0.24, -2.32, 0.60, -0.14, 1.30),
        c(3.52, 2.50, -3.50, 1.99, -0.17),
        c(3.05, 1.62, 1.04, -1.15, 1.61),
        c(3.98, 0.93, 1.93, -2.46, 0.75),
        c(0.84, -1.67, 3.71, 0.18, -2.65),
        c(1.75, 0.50, -1.44, -1.34, 0.66),
        c(3.11, 0.26, -0.11, -3.04, -0.25),
        c(2.05, -4.06, 0.36, -0.82, -0.38),
        c(2.47, 1.95, 0.26, 3.90, 0.09),
        c(-3.89, -1.73, -1.71, -0.84, 0.26),
        c(-4.28, -1.30, -1.49, -0.72, 0.84),
        c(2.29, 0.89, -2.49, 1.49, 0.31),
        c(-2.85, -0.22, 0.47, 1.94, -0.98),
        c(-4.22, 1.94, 1.06, 0.54, -0.62),
        c(-1.66, 0.27, 1.84, 0.70, 2.00),
        c(2.39, -1.07, 1.15, -1.39, 0.67),
        c(0.75, -2.18, -1.12, -1.46, -0.40),
        c(-4.36, 3.94, 0.59, 3.44, -1.59),
        c(-2.54, 2.44, 0.43, 0.04, -1.47),
        c(-2.59, -2.64, -1.54, -0.85, -0.02)

      )

      colnames(z_desc_matrix) <- c("z1", "z2", "z3", "z4", "z5")
      rownames(z_desc_matrix) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

      # Calculate Euclidian distances between 20 amino acids based on the
      # z-descriptor values

      s_dist_matrix <- as.matrix(dist(z_desc_matrix, method="euclidian", diag=T, upper=T))

      ## For input files with nucleotide sequences ##

      # Print error message if amino acid calculation is not specified

      if(input_seq == "nucl" && (is.null(aa_dist) || isFALSE(aa_dist))) {

        # Print error message

        stop("Sandberg distances can only be calculated for amino acid sequences. Please load amino acid sequences or specify aa_dist = T.")

      }

      if(input_seq == "nucl" && isTRUE(aa_dist)) {

        # Create five matrices for outputting z-descriptor values
        # for all sequences in the data set

        z1_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=max(lengths(seq_list_aa))))
        rownames(z1_matrix) <- seq_names
        colnames(z1_matrix) <- as.character(c(1:max(lengths(seq_list_aa))))
        z2_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=max(lengths(seq_list_aa))))
        rownames(z2_matrix) <- seq_names
        colnames(z2_matrix) <- as.character(c(1:max(lengths(seq_list_aa))))
        z3_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=max(lengths(seq_list_aa))))
        rownames(z3_matrix) <- seq_names
        colnames(z3_matrix) <- as.character(c(1:max(lengths(seq_list_aa))))
        z4_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=max(lengths(seq_list_aa))))
        rownames(z4_matrix) <- seq_names
        colnames(z4_matrix) <- as.character(c(1:max(lengths(seq_list_aa))))
        z5_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=max(lengths(seq_list_aa))))
        rownames(z5_matrix) <- seq_names
        colnames(z5_matrix) <- as.character(c(1:max(lengths(seq_list_aa))))

        # Sandberg distance calculation

        if(is.null(codon_pos)) {

          # Calculate the Sandberg distances in each pairwise comparison of
          # the sequences in seq_list_aa

          for (i in 1:(length(colnames(seq_file))-1)) {

            for (j in (i+1):length(colnames(seq_file))) {

              s_dist_ij <- vector("numeric", length=length(seq_list_aa[[i]]))

              for(k in 1:length(seq_list_aa[[i]])) {

                s_dist_ij[k] <- s_dist_matrix[which(rownames(s_dist_matrix)==seq_list_aa[[i]][k]),which(colnames(s_dist_matrix)==seq_list_aa[[j]][k])]

              }

              dist_matrix[i,j] <- round(mean(s_dist_ij), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_list_aa))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the Sandberg distance in each pairwise
            # comparison of the sequences in seq_list_aa, using
            # the vector codon_pos to define which codons to compare

            for (i in 1:(length(colnames(seq_file))-1)) {

              for (j in (i+1):length(colnames(seq_file))) {

                s_dist_ij <- vector("numeric", length=length(codon_pos))

                for(k in 1:length(codon_pos)) {

                  s_dist_ij[k] <- s_dist_matrix[which(rownames(s_dist_matrix)==seq_list_aa[[i]][codon_pos[k]]),which(colnames(s_dist_matrix)==seq_list_aa[[j]][codon_pos[k]])]

                }

                dist_matrix[i,j] <- round(mean(s_dist_ij), digits=5)

              }

            }

          }

        }

        # add z-descriptor values to the z1-5 matrices

        for(i in 1:length(colnames(seq_file))) {

          for(k in 1:max(lengths(seq_list_aa))) {

            z1_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list_aa[[i]][k]),1]
            z2_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list_aa[[i]][k]),2]
            z3_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list_aa[[i]][k]),3]
            z4_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list_aa[[i]][k]),4]
            z5_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list_aa[[i]][k]),5]

          }

        }

      }

      ## For input files with amino acid sequences ##

      if(input_seq == "aa") {

        # Create five matrices for outputting z-descriptor values
        # for all sequences in the data set

        z1_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=max(lengths(seq_list))))
        rownames(z1_matrix) <- seq_names
        colnames(z1_matrix) <- as.character(c(1:max(lengths(seq_list))))
        z2_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=max(lengths(seq_list))))
        rownames(z2_matrix) <- seq_names
        colnames(z2_matrix) <- as.character(c(1:max(lengths(seq_list))))
        z3_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=max(lengths(seq_list))))
        rownames(z3_matrix) <- seq_names
        colnames(z3_matrix) <- as.character(c(1:max(lengths(seq_list))))
        z4_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=max(lengths(seq_list))))
        rownames(z4_matrix) <- seq_names
        colnames(z4_matrix) <- as.character(c(1:max(lengths(seq_list))))
        z5_matrix <- as.data.frame(matrix(nrow=length(colnames(seq_file)),ncol=max(lengths(seq_list))))
        rownames(z5_matrix) <- seq_names
        colnames(z5_matrix) <- as.character(c(1:max(lengths(seq_list))))

        # Sandberg distance calculation

        if(is.null(codon_pos)) {

          # Calculate the Sandberg distances in each pairwise comparison of
          # the sequences in seq_list

          for (i in 1:(length(colnames(seq_file))-1)) {

            for (j in (i+1):length(colnames(seq_file))) {

              s_dist_ij <- vector("numeric", length=length(seq_list[[i]]))

              for(k in 1:length(seq_list[[i]])) {

                s_dist_ij[k] <- s_dist_matrix[which(rownames(s_dist_matrix)==seq_list[[i]][k]),which(colnames(s_dist_matrix)==seq_list[[j]][k])]

              }

              dist_matrix[i,j] <- round(mean(s_dist_ij), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_list))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the Sandberg distance in each pairwise
            # comparison of the sequences in seq_list, using
            # the vector codon_pos to define which codons to compare

            for (i in 1:(length(colnames(seq_file))-1)) {

              for (j in (i+1):length(colnames(seq_file))) {

                s_dist_ij <- vector("numeric", length=length(codon_pos))

                for(k in 1:length(codon_pos)) {

                  s_dist_ij[k] <- s_dist_matrix[which(rownames(s_dist_matrix)==seq_list[[i]][codon_pos[k]]),which(colnames(s_dist_matrix)==seq_list[[j]][codon_pos[k]])]

                }

                dist_matrix[i,j] <- round(mean(s_dist_ij), digits=5)

              }

            }

          }

        }

        # add z-descriptor values to the z1-5 matrices

        for(i in 1:length(colnames(seq_file))) {

          for(k in 1:max(lengths(seq_list))) {

            z1_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list[[i]][k]),1]
            z2_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list[[i]][k]),2]
            z3_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list[[i]][k]),3]
            z4_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list[[i]][k]),4]
            z5_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list[[i]][k]),5]

          }

        }

      }

      # Export the z1-5 matrices as .csv files
      write.csv(z1_matrix,file=paste(path_out,"/z1_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))
      write.csv(z2_matrix,file=paste(path_out,"/z2_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))
      write.csv(z3_matrix,file=paste(path_out,"/z3_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))
      write.csv(z4_matrix,file=paste(path_out,"/z4_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))
      write.csv(z5_matrix,file=paste(path_out,"/z5_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))

    }

    ### Calculate mean distances for each sample in the data set

    for (i in 1:length(sample_names)) {

      # Fetch column numbers for the sequences in sample i in the sequence
      # table

      z <- which(seq_file[i,] > 0)

      # Create a vector md

      md <- vector("numeric", length=length(z))

      # Generate a list of all pairwise combinations of the elements in z

      pwc <- combn(sort(z),2,simplify=T)

      # For each combination in pwc

      for(j in 1:length(pwc[1,]))  {

        # Extract the distance in each pairwise comparison of the sequences
        # in seq_list from the dist_matrix, using the numbers from pwc as
        # indices to extract the values from the matrix

        md[j] <- dist_matrix[pwc[1,j],pwc[2,j]]

      }

      # Calculate the mean of the md vector and add this to the mean_dist vector

      mean_dist[i] <- mean(md)

    }

    # Export a table with sample names and mean distances as a .csv file

    mean_dist_tab <- as.data.frame(mean_dist)

    rownames(mean_dist_tab) <- sample_names

    write.csv(mean_dist_tab,file=paste(path_out,"/mean_dist_table_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))

  }

  #################################
  ###### Fasta file as input ######
  #################################

  if(isTRUE(input_fasta)) {

    # fasta files are accepted in the list format rendered by
    # the read.fasta() function from the package 'seqinr'

    # create a matrix that will contain the pairwise distances between all the sequences in the fasta file

    dist_matrix <- as.data.frame(matrix(nrow=length(seq_file),ncol=length(seq_file)))
    rownames(dist_matrix) <- names(seq_file)
    colnames(dist_matrix) <- names(seq_file)

    # Throw a warning if sequences in the fasta file are of different lengths

    if(max(lengths(seq_file)) != min(lengths(seq_file))) {

      stop("Pairwise comparisons not meaningful for sequences of different length.")

    }

    # Throw a warning if sequences in the fasta file contain non-standard characters

    if(input_seq == "aa") {

      for(i in 1:length(seq_file)) {

        for(j in 1:length(unique(seq_file[[i]]))) {

          if(!(unique(seq_file[[i]])[j] %in% c("-","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"))) {

            stop("Sequences contain non-standard characters. Accepted characters are -ARNDCQEGHILKMFPSTWYV.")

          }

        }

      }

    }

    if(input_seq == "nucl") {

      for(i in 1:length(seq_file)) {

        for(j in 1:length(unique(seq_file[[i]]))) {

          if(!(toupper(unique(seq_file[[i]])[j]) %in% c("-","A","C","G","T"))) {

            stop("Sequences contain non-standard characters. Accepted characters are -ATGC.")

          }

        }

      }

    }

    ### Create a vector seq_list_aa containing all the sequences in the fasta file translated to amino acids

    if(isTRUE(aa_dist)) {

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
              if(list(c("-")) %in% paste(seq_file[[i]][c(j*3-2,j*3-1,j*3)])) {seq_list_aa[[i]][j] <- "X"; warning("Gap detected in sequence.")}

            }

          } else {

            stop("Sequences should be aligned in open reading frame 5' -> 3'.")

          }

        }

      } else {

        stop("Sequences should be aligned in open reading frame 5' -> 3'.")

      }

    }

    ### calculate pairwise distances between all sequences in the
    ### fasta file, and save the values in the matrix
    ### (upper right matrix, values rounded to five digits)

    ##############################
    ### p-distance calculation ###
    ##############################

    if(dist_type == "P") {

      ## For input files with nucleotide sequences ##

      if(input_seq == "nucl" && (is.null(aa_dist) || isFALSE(aa_dist))) {

        # nucleotide p-distance calculation

        if(is.null(codon_pos)) {

          # Calculate the nucleotide Pdist in each pairwise comparison of the sequences
          # in seq_file

          for (i in 1:(length(names(seq_file))-1)) {

            for (j in (i+1):length(names(seq_file))) {

              dist_matrix[i,j] <- round(length(which(seq_file[[i]] != seq_file[[j]]))/length(seq_file[[i]]), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed sequence length

          if(max(codon_pos) > max(lengths(seq_file))) {

            stop("Selected codons exceed sequence length.")

          } else {

            # Calculate the nucleotide Pdist in each pairwise comparison of the sequences
            # in seq_list, using the vector codon_pos to define which codons to compare

            for (i in 1:(length(names(seq_file))-1)) {

              for (j in (i+1):length(names(seq_file))) {

                dist_matrix[i,j] <- round(length(which(seq_file[[i]][codon_pos] != seq_file[[j]][codon_pos]))/length(seq_file[[i]][codon_pos]), digits=5)

              }

            }

          }

        }

      }

      if(input_seq == "nucl" && isTRUE(aa_dist)) {

        # amino acid p-distance calculation

        if(is.null(codon_pos)) {

          # Calculate the amino acid Pdist in each pairwise comparison of the sequences
          # in seq_list_aa

          for (i in 1:(length(names(seq_file))-1)) {

            for (j in (i+1):length(names(seq_file))) {

              dist_matrix[i,j] <- round(length(which(seq_list_aa[[i]] != seq_list_aa[[j]]))/length(seq_list_aa[[i]]), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_list_aa))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the amino acid Pdist in each pairwise comparison of the sequences
            # in seq_list_aa, using the vector codon_pos to define  which codons to
            # compare

            for (i in 1:(length(names(seq_file))-1)) {

              for (j in (i+1):length(names(seq_file))) {

                dist_matrix[i,j] <- round(length(which(seq_list_aa[[i]][codon_pos] != seq_list_aa[[j]][codon_pos]))/length(seq_list_aa[[i]][codon_pos]), digits=5)

              }

            }

          }

        }

      }

      ## For input files with amino acid sequences ##

      if(input_seq == "aa") {

        # amino acid p-distance calculation

        if(is.null(codon_pos)) {

          # Calculate the amino acid Pdist in each pairwise comparison of the sequences
          # in seq_list

          for (i in 1:(length(names(seq_file))-1)) {

            for (j in (i+1):length(names(seq_file))) {

              dist_matrix[i,j] <- round(length(which(seq_file[[i]] != seq_file[[j]]))/length(seq_file[[i]]), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_file))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the amino acid Pdist in each pairwise comparison of the sequences
            # in seq_file, using the vector codon_pos to define which codons to compare

            for (i in 1:(length(names(seq_file))-1)) {

              for (j in (i+1):length(names(seq_file))) {

                dist_matrix[i,j] <- round(length(which(seq_file[[i]][codon_pos] != seq_file[[j]][codon_pos]))/length(seq_file[[i]][codon_pos]), digits=5)

              }

            }

          }

        }

      }

    }

    #####################################
    ### Grantham distance calculation ###
    #####################################

    if(dist_type == "G") {

      # Define a Grantham distance matrix

      g_dist_matrix <- rbind.data.frame(

        c(  0, 112, 111, 126, 195,  91, 107,  60,  86,  94,  96, 106,  84, 113,  27,  99,  58, 148, 112,  64),
        c(112,   0,  86,  96, 180,  43,  54, 125,  29,  97, 102,  26,  91,  97, 103, 110,  71, 101,  77,  96),
        c(111,  86,   0,  23, 139,  46,  42,  80,  68, 149, 153,  94, 142, 158,  91,  46,  65, 174, 143, 133),
        c(126,  96,  23,   0, 154,  61,  45,  94,  81, 168, 172, 101, 160, 177, 108,  65,  85, 181, 160, 152),
        c(195, 180, 139, 154,   0, 154, 170, 159, 174, 198, 198, 202, 196, 205, 169, 112, 149, 215, 194, 192),
        c( 91,  43,  46,  61, 154,   0,  29,  87,  24, 109, 113,  53, 101, 116,  76,  68,  42, 130,  99,  96),
        c(107,  54,  42,  45, 170,  29,   0,  98,  40, 134, 138,  56, 126, 140,  93,  80,  65, 152, 122, 121),
        c( 60, 125,  80,  94, 159,  87,  98,   0,  98, 135, 138, 127, 127, 153,  42,  56,  59, 184, 147, 109),
        c( 86,  29,  68,  81, 174,  24,  40,  98,   0,  94,  99,  32,  87, 100,  77,  89,  47, 115,  83,  84),
        c( 94,  97, 149, 168, 198, 109, 134, 135,  94,   0,   5, 102,  10,  21,  95, 142,  89,  61,  33,  29),
        c( 96, 102, 153, 172, 198, 113, 138, 138,  99,   5,   0, 107,  15,  22,  98, 145,  92,  61,  36,  32),
        c(106,  26,  94, 101, 202,  53,  56, 127,  32, 102, 107,   0,  95, 102, 103, 121,  78, 110,  85,  97),
        c( 84,  91, 142, 160, 196, 101, 126, 127,  87,  10,  15,  95,   0,  28,  87, 135,  81,  67,  36,  21),
        c(113,  97, 158, 177, 205, 116, 140, 153, 100,  21,  22, 102,  28,   0, 114, 155, 103,  40,  22,  50),
        c( 27, 103,  91, 108, 169,  76,  93,  42,  77,  95,  98, 103,  87, 114,   0,  74,  38, 147, 110,  68),
        c( 99, 110,  46,  65, 112,  68,  80,  56,  89, 142, 145, 121, 135, 155,  74,   0,  58, 177, 144, 124),
        c( 58,  71,  65,  85, 149,  42,  65,  59,  47,  89,  92,  78,  81, 103,  38,  58,   0, 128,  92,  69),
        c(148, 101, 174, 181, 215, 130, 152, 184, 115,  61,  61, 110,  67,  40, 147, 177, 128,   0,  37,  88),
        c(112,  77, 143, 160, 194,  99, 122, 147,  83,  33,  36,  85,  36,  22, 110, 144,  92,  37,   0,  55),
        c( 64,  96, 133, 152, 192,  96, 121, 109,  84,  29,  32,  97,  21,  50,  68, 124,  69,  88,  55,   0)

      )

      colnames(g_dist_matrix) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
      rownames(g_dist_matrix) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

      ## For input files with nucleotide sequences ##

      # Print error message if amino acid calculation is not specified

      if(input_seq == "nucl" && (is.null(aa_dist) || isFALSE(aa_dist))) {

        stop("Grantham distances can only be calculated for amino acid sequences. Please load amino acid sequences or specify aa_dist = T.")

      }

      if(input_seq == "nucl" && isTRUE(aa_dist)) {

        # Grantham distance calculation

        if(is.null(codon_pos)) {

          # Calculate the Grantham distances in each pairwise comparison of
          # the sequences in seq_list_aa

          for (i in 1:(length(names(seq_file))-1)) {

            for (j in (i+1):length(names(seq_file))) {

              g_dist_ij <- vector("numeric", length=length(seq_list_aa[[i]]))

              for(k in 1:length(seq_list_aa[[i]])) {

                g_dist_ij[k] <- g_dist_matrix[which(rownames(g_dist_matrix)==seq_list_aa[[i]][k]),which(colnames(g_dist_matrix)==seq_list_aa[[j]][k])]

              }

              dist_matrix[i,j] <- round(mean(g_dist_ij), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_list_aa))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the Grantham distance in each pairwise
            # comparison of the sequences in seq_list_aa, using
            # the vector codon_pos to define which codons to compare

            for (i in 1:(length(names(seq_file))-1)) {

              for (j in (i+1):length(names(seq_file))) {

                g_dist_ij <- vector("numeric", length=length(codon_pos))

                for(k in 1:length(codon_pos)) {

                  g_dist_ij[k] <- g_dist_matrix[which(rownames(g_dist_matrix)==seq_list_aa[[i]][codon_pos[k]]),which(colnames(g_dist_matrix)==seq_list_aa[[j]][codon_pos[k]])]

                }

                dist_matrix[i,j] <- round(mean(g_dist_ij), digits=5)

              }

            }

          }

        }

      }

      ## For input files with amino acid sequences ##

      if(input_seq == "aa") {

        # Grantham distance calculation

        if(is.null(codon_pos)) {

          # Calculate the Grantham distances in each pairwise comparison of
          # the sequences in seq_file

          for (i in 1:(length(names(seq_file))-1)) {

            for (j in (i+1):length(names(seq_file))) {

              g_dist_ij <- vector("numeric", length=length(seq_file[[i]]))

              for(k in 1:length(seq_file[[i]])) {

                g_dist_ij[k] <- g_dist_matrix[which(rownames(g_dist_matrix)==seq_file[[i]][k]),which(colnames(g_dist_matrix)==seq_file[[j]][k])]

              }

              dist_matrix[i,j] <- round(mean(g_dist_ij), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_file))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the Grantham distance in each pairwise
            # comparison of the sequences in seq_file, using
            # the vector codon_pos to define which codons to compare

            for (i in 1:(length(names(seq_file))-1)) {

              for (j in (i+1):length(names(seq_file))) {

                g_dist_ij <- vector("numeric", length=length(codon_pos))

                for(k in 1:length(codon_pos)) {

                  g_dist_ij[k] <- g_dist_matrix[which(rownames(g_dist_matrix)==seq_file[[i]][codon_pos[k]]),which(colnames(g_dist_matrix)==seq_file[[j]][codon_pos[k]])]

                }

                dist_matrix[i,j] <- round(mean(g_dist_ij), digits=5)

              }

            }

          }

        }

      }

    }

    #####################################
    ### Sandberg distance calculation ###
    #####################################

    if(dist_type == "S") {

      # Calculate Sandberg distances as the Euclidian distances between all 20 amino
      # acids based on the five physicochemical z-descriptors described in Sandberg
      # et al. (1998)

      # Define a matrix with z-descriptor values for all 20 amino acids

      z_desc_matrix <- rbind.data.frame(

        c(0.24, -2.32, 0.60, -0.14, 1.30),
        c(3.52, 2.50, -3.50, 1.99, -0.17),
        c(3.05, 1.62, 1.04, -1.15, 1.61),
        c(3.98, 0.93, 1.93, -2.46, 0.75),
        c(0.84, -1.67, 3.71, 0.18, -2.65),
        c(1.75, 0.50, -1.44, -1.34, 0.66),
        c(3.11, 0.26, -0.11, -3.04, -0.25),
        c(2.05, -4.06, 0.36, -0.82, -0.38),
        c(2.47, 1.95, 0.26, 3.90, 0.09),
        c(-3.89, -1.73, -1.71, -0.84, 0.26),
        c(-4.28, -1.30, -1.49, -0.72, 0.84),
        c(2.29, 0.89, -2.49, 1.49, 0.31),
        c(-2.85, -0.22, 0.47, 1.94, -0.98),
        c(-4.22, 1.94, 1.06, 0.54, -0.62),
        c(-1.66, 0.27, 1.84, 0.70, 2.00),
        c(2.39, -1.07, 1.15, -1.39, 0.67),
        c(0.75, -2.18, -1.12, -1.46, -0.40),
        c(-4.36, 3.94, 0.59, 3.44, -1.59),
        c(-2.54, 2.44, 0.43, 0.04, -1.47),
        c(-2.59, -2.64, -1.54, -0.85, -0.02)

      )

      colnames(z_desc_matrix) <- c("z1", "z2", "z3", "z4", "z5")
      rownames(z_desc_matrix) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

      # Calculate Euclidian distances between 20 amino acids based on the
      # z-descriptor values

      s_dist_matrix <- as.matrix(dist(z_desc_matrix, method="euclidian", diag=T, upper=T))

      ## For input files with nucleotide sequences ##

      # Print error message if amino acid calculation is not specified

      if(input_seq == "nucl" && (is.null(aa_dist) || isFALSE(aa_dist))) {

        stop("Sandberg distances can only be calculated for amino acid sequences. Please load amino acid sequences or specify aa_dist = T.")

      }

      if(input_seq == "nucl" && isTRUE(aa_dist)) {

        # Create five matrices for outputting z-descriptor values
        # for all sequences in the data set

        z1_matrix <- as.data.frame(matrix(nrow=length(names(seq_file)),ncol=max(lengths(seq_list_aa))))
        rownames(z1_matrix) <- names(seq_file)
        colnames(z1_matrix) <- as.character(c(1:max(lengths(seq_list_aa))))
        z2_matrix <- as.data.frame(matrix(nrow=length(names(seq_file)),ncol=max(lengths(seq_list_aa))))
        rownames(z2_matrix) <- names(seq_file)
        colnames(z2_matrix) <- as.character(c(1:max(lengths(seq_list_aa))))
        z3_matrix <- as.data.frame(matrix(nrow=length(names(seq_file)),ncol=max(lengths(seq_list_aa))))
        rownames(z3_matrix) <- names(seq_file)
        colnames(z3_matrix) <- as.character(c(1:max(lengths(seq_list_aa))))
        z4_matrix <- as.data.frame(matrix(nrow=length(names(seq_file)),ncol=max(lengths(seq_list_aa))))
        rownames(z4_matrix) <- names(seq_file)
        colnames(z4_matrix) <- as.character(c(1:max(lengths(seq_list_aa))))
        z5_matrix <- as.data.frame(matrix(nrow=length(names(seq_file)),ncol=max(lengths(seq_list_aa))))
        rownames(z5_matrix) <- names(seq_file)
        colnames(z5_matrix) <- as.character(c(1:max(lengths(seq_list_aa))))

        # Sandberg distance calculation

        if(is.null(codon_pos)) {

          # Calculate the Sandberg distances in each pairwise comparison of
          # the sequences in seq_list_aa

          for (i in 1:(length(names(seq_file))-1)) {

            for (j in (i+1):length(names(seq_file))) {

              s_dist_ij <- vector("numeric", length=length(seq_list_aa[[i]]))

              for(k in 1:length(seq_list_aa[[i]])) {

                s_dist_ij[k] <- s_dist_matrix[which(rownames(s_dist_matrix)==seq_list_aa[[i]][k]),which(colnames(s_dist_matrix)==seq_list_aa[[j]][k])]

              }

              dist_matrix[i,j] <- round(mean(s_dist_ij), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_list_aa))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the Sandberg distance in each pairwise
            # comparison of the sequences in seq_list_aa, using
            # the vector codon_pos to define which codons to compare

            for (i in 1:(length(names(seq_file))-1)) {

              for (j in (i+1):length(names(seq_file))) {

                s_dist_ij <- vector("numeric", length=length(codon_pos))

                for(k in 1:length(codon_pos)) {

                  s_dist_ij[k] <- s_dist_matrix[which(rownames(s_dist_matrix)==seq_list_aa[[i]][codon_pos[k]]),which(colnames(s_dist_matrix)==seq_list_aa[[j]][codon_pos[k]])]

                }

                dist_matrix[i,j] <- round(mean(s_dist_ij), digits=5)

              }

            }

          }

        }

        # add z-descriptor values to the z1-5 matrices

        for(i in 1:length(names(seq_file))) {

          for(k in 1:max(lengths(seq_list_aa))) {

            z1_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list_aa[[i]][k]),1]
            z2_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list_aa[[i]][k]),2]
            z3_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list_aa[[i]][k]),3]
            z4_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list_aa[[i]][k]),4]
            z5_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_list_aa[[i]][k]),5]

          }

        }

      }

      ## For input files with amino acid sequences ##

      if(input_seq == "aa") {

        # Create five matrices for outputting z-descriptor values
        # for all sequences in the data set

        z1_matrix <- as.data.frame(matrix(nrow=length(names(seq_file)),ncol=max(lengths(seq_file))))
        rownames(z1_matrix) <- names(seq_file)
        colnames(z1_matrix) <- as.character(c(1:max(lengths(seq_file))))
        z2_matrix <- as.data.frame(matrix(nrow=length(names(seq_file)),ncol=max(lengths(seq_file))))
        rownames(z2_matrix) <- names(seq_file)
        colnames(z2_matrix) <- as.character(c(1:max(lengths(seq_file))))
        z3_matrix <- as.data.frame(matrix(nrow=length(names(seq_file)),ncol=max(lengths(seq_file))))
        rownames(z3_matrix) <- names(seq_file)
        colnames(z3_matrix) <- as.character(c(1:max(lengths(seq_file))))
        z4_matrix <- as.data.frame(matrix(nrow=length(names(seq_file)),ncol=max(lengths(seq_file))))
        rownames(z4_matrix) <- names(seq_file)
        colnames(z4_matrix) <- as.character(c(1:max(lengths(seq_file))))
        z5_matrix <- as.data.frame(matrix(nrow=length(names(seq_file)),ncol=max(lengths(seq_file))))
        rownames(z5_matrix) <- names(seq_file)
        colnames(z5_matrix) <- as.character(c(1:max(lengths(seq_file))))

        # Sandberg distance calculation

        if(is.null(codon_pos)) {

          # Calculate the Sandberg distances in each pairwise comparison of
          # the sequences in seq_file

          for (i in 1:(length(names(seq_file))-1)) {

            for (j in (i+1):length(names(seq_file))) {

              s_dist_ij <- vector("numeric", length=length(seq_file[[i]]))

              for(k in 1:length(seq_file[[i]])) {

                s_dist_ij[k] <- s_dist_matrix[which(rownames(s_dist_matrix)==seq_file[[i]][k]),which(colnames(s_dist_matrix)==seq_file[[j]][k])]

              }

              dist_matrix[i,j] <- round(mean(s_dist_ij), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed amino acid sequence length

          if(max(codon_pos) > max(lengths(seq_file))) {

            stop("Selected codons exceed amino acid sequence length.")

          } else {

            # Calculate the Sandberg distance in each pairwise
            # comparison of the sequences in seq_file, using
            # the vector codon_pos to define which codons to compare

            for (i in 1:(length(names(seq_file))-1)) {

              for (j in (i+1):length(names(seq_file))) {

                s_dist_ij <- vector("numeric", length=length(codon_pos))

                for(k in 1:length(codon_pos)) {

                  s_dist_ij[k] <- s_dist_matrix[which(rownames(s_dist_matrix)==seq_file[[i]][codon_pos[k]]),which(colnames(s_dist_matrix)==seq_file[[j]][codon_pos[k]])]

                }

                dist_matrix[i,j] <- round(mean(s_dist_ij), digits=5)

              }

            }

          }

        }

        # add z-descriptor values to the z1-5 matrices

        for(i in 1:length(names(seq_file))) {

          for(k in 1:max(lengths(seq_file))) {

            z1_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_file[[i]][k]),1]
            z2_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_file[[i]][k]),2]
            z3_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_file[[i]][k]),3]
            z4_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_file[[i]][k]),4]
            z5_matrix[i,k] <- z_desc_matrix[which(rownames(z_desc_matrix)==seq_file[[i]][k]),5]

          }

        }

      }

      # Export the z1-5 matrices as .csv files
      write.csv(z1_matrix,file=paste(path_out,"/z1_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))
      write.csv(z2_matrix,file=paste(path_out,"/z2_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))
      write.csv(z3_matrix,file=paste(path_out,"/z3_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))
      write.csv(z4_matrix,file=paste(path_out,"/z4_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))
      write.csv(z5_matrix,file=paste(path_out,"/z5_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))

    }

  }

  # Export the distance matrix as a .csv file
  write.csv(dist_matrix,file=paste(path_out,"/dist_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv", sep = ""))

}
