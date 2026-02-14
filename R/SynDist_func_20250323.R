#' SynDist() function
#'
#' \code{\link{SynDist}} identifies and quantifies synonymous variation
#' among aligned protein-coding DNA sequences, that is, nucleotide
#' substitutions that do not translate to changes in the amino acid
#' sequences, due to degeneracy of the genetic code.
#'
#' The SynDist() function takes a fasta file or a 'dada2'-style sequence
#' occurrence table (with aligned sequences as column names and samples in
#' rows) as input and identifies synonymous variation by pairwise sequence
#' comparisons.
#'
#' SynDist() can do qualitative or quantitative analysis of synonymous
#' variation. If analysis="codon" is specified, the function identifies
#' synonymous nucleotide variation and outputs tables with the number of
#' observations of synonymous nucleotide changes per base and per codon
#' among all pairwise sequence comparisons in the data set. These tables
#' also specify, for each base or codon position, the proportions of the
#' total pairwise comparisons that harbor synonymous substitutions.
#'
#' If analysis="dist", the function produces a distance matrix specifying
#' the number and proportion (p-distance) of synonymous nucleotide
#' changes in each pairwise sequence comparison in the data set. In the
#' distance matrix, synonymous p-distance is calculated as the number of
#' synonymous nucleotide changes observed in each pairwise sequence
#' comparison divided by the sequence length (number of bases). If a
#' 'dada2'-style sequence occurrence table is provided as input, the
#' SynDist() function furthermore produces two tables with the mean number
#' of synonymous variations and mean synonymous p-distances among all
#' pairwise comparisons of the sequences in each sample in the data set.
#' (Note: The means will be NA for samples that have 0 or 1 sequence(s).)
#'
#' The SynDist() function includes an option for the user to specify which
#' codons to compare. This is useful e.g. if the sequences contain gaps in
#' some codons, which should be excluded from quantitative analysis.
#'
#' SynDist() translates the supplied DNA sequences to amino acid sequences
#' using the standard genetic code and sequences must be aligned in open
#' reading frame. The function only accepts the following characters in the
#' sequences: -,a,t,g,c,A,T,G,C
#'
#' Nucleotide triplets containing gaps (-) are translated to 'X', similar to
#' stop codons. Please note that '-' are treated as unique characters in p-
#' distance calculations. The function will give warnings if gaps or stop
#' codons are detected. If you wish to exclude stop codons or gaps from
#' distance calculations, please use the codon_pos option to specify which
#' codons to compare.
#'
#' If you publish data or results produced with MHCtools, please cite both of
#' the following references:
#' Roved, J. (2022). MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J. (2024). MHCtools 1.5: Analysis of MHC sequencing data in R. In S.
#' Boegel (Ed.), HLA Typing: Methods and Protocols (2nd ed., pp. 275–295).
#' Humana Press. https://doi.org/10.1007/978-1-0716-3874-3_18
#'
#' @param seq_file is a sequence occurrence table as output by the 'dada2'
#'   pipeline, which has samples in rows and nucleotide sequences in columns.
#'   Optionally, a fasta file can be supplied as input in the format rendered
#'   by read.fasta() from the package 'seqinr'.
#' @param path_out is a user defined path to the folder where the output files
#'   will be saved.
#' @param input_fasta optional, a logical (TRUE/FALSE) that indicates whether
#'   the input file is a fasta file (TRUE) or a 'dada2'-style sequence table
#'   (NULL/FALSE). The default is NULL/FALSE.
#' @param codon_pos is optional, a vector of comma separated integers specifying
#'   which codons to include in analyses. If omitted, analyses are made using
#'   all codons. Note: With SynDist(), codon_pos should always be specified as
#'   codons, i.e. numbered nucleotide triplets in open reading frame.
#' @param analysis is used to specify the desired kind of analysis. It takes
#'   the values 'dist' for quantification of pairwise synonymous variation
#'   between sequences, or 'codon' for quantification of synonymous substitutions
#'   per nucleotide or codon position. The argument is optional with 'dist' as
#'   default.
#' @return When analysis="dist", the function produces a .csv distance matrix
#'   with the number of synonymous substitutions in each pairwise sequence
#'   comparison in the upper right matrix and the synonymous p-distance in each
#'   pairwise sequence comparison in the lower left matrix. If a sequence
#'   occurrence table is given as input file, the function additionally produces
#'   two tables with the mean number of synonymous substitutions and the mean
#'   synonymous p-distance across all pairwise sequence comparisons for each
#'   sample in the data set. If a sequence occurrence table is given as input
#'   file, the sequences are named in the output matrix by an index number that
#'   corresponds to their column number in the input file.
#'   If analysis="codon", the function produces two .csv summary tables, one with
#'   the total number of synonymous substitutions per nucleotide position across
#'   all pairwise sequence comparisons and one with the number of synonymous
#'   codon variations per codon across all pairwise sequence comparisons. Note
#'   that in the codon summary table, the synonymous codon variation does not
#'   quantify the number of nucleotide variations between the synonymous codons,
#'   since that can be derived from the nucleotide summary table. Each summary
#'   table also contains a column that specifies the proportion of the observed
#'   number of synonymous variations (per nucleotide position or codon) out of the
#'   number of pairwise sequence comparisons. E.g., if three sequences are
#'   compared and a synonymous substitution is observed for a given codon once
#'   (i.e., between two of the three sequences), that gives a proportion of
#'   synonymous observations of one out of three pairwise sequence comparisons for
#'   that codon.
#' @seealso For more information about 'dada2', visit
#'   <https://benjjneb.github.io/dada2/>
#' @examples
#' seq_file <- sequence_table_SynDist
#' path_out <- tempdir()
#' SynDist(seq_file, path_out, input_fasta=NULL,codon_pos=c(1,2,3,4,5,6,7,8),
#' analysis="dist")
#' @export

SynDist <- function(seq_file, path_out, input_fasta=NULL, codon_pos=NULL, analysis="dist") {

  print("Note: For correct analysis, sequences must be aligned in open reading frame 5' -> 3'")

    ##########################################
    ###### Sequence table as input file ######
    ##########################################

    if(is.null(input_fasta) || isFALSE(input_fasta)) {

      # The dada2 sequence table does not use sequences names, but identifies
      # sequence variants by their nuceotide sequence. Here I create a vector for
      # naming the sequences by their column number in the seq_table

      seq_names <- vector("character", length=length(colnames(seq_file)))

      seq_names <- paste0("Sequence_", formatC(seq(1:length(colnames(seq_file))), width = nchar(length(colnames(seq_file))), format = "d", flag = "0"))
      # the formatC() expression creates index numbers of the sequences with zeroes
      # padded in front, so that all numbers have the same number of digits to prevent
      # RegEx pattern matching between e.g. "Sequence_1" and "Sequence_1X".

      # Extract the sample names to a new vector

      sample_names <- rownames(seq_file)

      # Create vectors mean_diff_bases and mean_pdist that will get the mean distances per individual

      mean_diff_bases <- vector("numeric", length=length(sample_names))
      mean_pdist <- vector("numeric", length=length(sample_names))

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

      for(i in 1:length(seq_list)) {

        for(j in 1:length(unique(seq_list[[i]]))) {

          if(!(toupper(unique(seq_list[[i]])[j]) %in% c("-","A","C","G","T"))) {

            stop("Sequences contain non-standard characters. Accepted characters are -atgcATGC.")

          }

        }

      }

      ### Create a vector seq_list_aa containing all the sequences in seq_file translated to amino acids

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

      for(i in 1:length(seq_list)) {

        # fill in NaNs in seq_list_aa
        seq_list_aa[[i]] <- NaN*seq(length(seq_list[[i]])/3)

        # check whether sequence lengths are divisible by 3
        if((length(seq_list[[i]])/3) %% 1 == 0) {

          for(j in 1:(length(seq_list[[i]])/3)) {

            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Phe) {seq_list_aa[[i]][j] <- "F"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Leu) {seq_list_aa[[i]][j] <- "L"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Ile) {seq_list_aa[[i]][j] <- "I"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Met) {seq_list_aa[[i]][j] <- "M"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Val) {seq_list_aa[[i]][j] <- "V"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Ser) {seq_list_aa[[i]][j] <- "S"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Pro) {seq_list_aa[[i]][j] <- "P"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Thr) {seq_list_aa[[i]][j] <- "T"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Ala) {seq_list_aa[[i]][j] <- "A"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Tyr) {seq_list_aa[[i]][j] <- "Y"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% His) {seq_list_aa[[i]][j] <- "H"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Gln) {seq_list_aa[[i]][j] <- "Q"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Asn) {seq_list_aa[[i]][j] <- "N"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Lys) {seq_list_aa[[i]][j] <- "K"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Asp) {seq_list_aa[[i]][j] <- "D"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Glu) {seq_list_aa[[i]][j] <- "E"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Cys) {seq_list_aa[[i]][j] <- "C"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Trp) {seq_list_aa[[i]][j] <- "W"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Arg) {seq_list_aa[[i]][j] <- "R"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Gly) {seq_list_aa[[i]][j] <- "G"}
            if(toupper(paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)],collapse="")) %in% Stop) {seq_list_aa[[i]][j] <- "X"; warning("Stop codon detected in sequence.")}
            if(list(c("-")) %in% paste(seq_list[[i]][c(j*3-2,j*3-1,j*3)])) {seq_list_aa[[i]][j] <- "X"; warning("Gap detected in sequence.")}

          }

        } else {

          stop("Sequences should be aligned in open reading frame 5' -> 3'")

        }

      }

      #####################################################
      ### Identify codons with synonymous substitutions ###
      #####################################################

      ### Option to output tables with number of pairwise comparisons
      ### with synonymous nucleotide changes per base and per codon in
      ### the data set

      if(analysis == "codon") {

        # Create table for presence/absence of synonymous nucleotide
        # substitutions in all pairwise comparisons
        n_pairs <- length(seq_list)*(length(seq_list)-1)/2 # number of pairwise sequence comparisons
        n_nucl <- max(lengths(seq_list)) # sequence length in nucleotides
        SynSubTab <- matrix(data=rep(0,n_pairs*n_nucl), nrow=n_pairs, ncol=n_nucl)

        # Create table for presence/absence of synonymous codons
        # in all pairwise comparisons
        n_cod <- max(lengths(seq_list)/3) # sequence length in codons
        SynCodTab <- matrix(data=rep(0,n_pairs*n_cod), nrow=n_pairs, ncol=n_cod)

        # Create two vectors to sum the presence/absence of synonymous nucleotide
        # substitutions or codons in all pairwise comparisons
        SynSub_sum <- vector()
        SynCod_sum <- vector()
        # Create a vector that defines nucleotide triplet positions
        TripPos <- as.character(rep(c(1,2,3),n_nucl/3))

        if(is.null(codon_pos)) {

          # create a vector for counting loop iterations
          iter <- numeric(0)

          # Identify codons with synonymous nucleotide substitutions

          for (i in 1:(length(colnames(seq_file))-1)) {

            for (j in (i+1):length(colnames(seq_file))) {

              iter <- sum(iter, 1)

              # Identify codons with no amino acid change
              InvCod <- which(seq_list_aa[[i]] == seq_list_aa[[j]])

              # Identify codons with nucleotide differences
              VarCod <- vector()

              for(k in 1:(length(seq_list[[i]])/3)) {

                if(toupper(paste(seq_list[[i]][c(k*3-2,k*3-1,k*3)],collapse="")) != toupper(paste(seq_list[[j]][c(k*3-2,k*3-1,k*3)],collapse=""))) {VarCod <- append(VarCod, k)}

              }

              # Identify codons with synonymous nucleotide substitutions
              SynCod <- VarCod[which(VarCod %in% InvCod)]

              # Get the base positions that correspond to the codons in SynCod
              SynCod_bases <- vector()

              if(length(SynCod) > 0) {

                for(k in 1:(length(SynCod))) {

                  SynCod_bases <- append(SynCod_bases,c(SynCod[k]*3-2,SynCod[k]*3-1,SynCod[k]*3))

                }

              }

              # Identify the positions harboring nucleotide changes in SynCod_bases
              VarSynCod_bases <- vector()

              if(length(SynCod) > 0) {

                for(k in 1:(length(SynCod_bases))) {

                  if(toupper(seq_list[[i]][SynCod_bases[k]]) != toupper(seq_list[[j]][SynCod_bases[k]])) {VarSynCod_bases <- append(VarSynCod_bases, SynCod_bases[k])}

                }

              }

              # Add presence/absence of synonymous nucleotide substitutions to
              # table with all pairwise comparisons
              # add 1 in the positions with synonymous nucleotide substitutions
              SynSubTab[iter,VarSynCod_bases] <- 1

              # similar for the codons
              SynCodTab[iter,SynCod] <- 1

            }

          }

          # Output a vector with column sums for each table
          # Include nucleotide triplet positions in the
          # table with synonymous nucleotide substitutions
          # Output also a vector with the proportions of
          # pairwise comparisons that harbored a synonymous
          # substitution in each base position or codon
          SynSub_sum <- cbind.data.frame(TripPos,colSums(SynSubTab),colSums(SynSubTab)/n_pairs)
          rownames(SynSub_sum) <- paste0("bp",1:n_nucl)
          colnames(SynSub_sum) <- c("Triplet.Pos","No.Syn.Obs","Prop.Syn.Obs")
          SynCod_sum <- cbind.data.frame(colSums(SynCodTab),colSums(SynCodTab)/n_pairs)
          rownames(SynCod_sum) <- paste0("cod",1:n_cod)
          colnames(SynCod_sum) <- c("No.Syn.Obs","Prop.Syn.Obs")

          # Save the vectors as .csv
          write.csv(SynSub_sum,file=paste0(path_out,"/Syn_subst_sum_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))
          write.csv(SynCod_sum,file=paste0(path_out,"/Syn_codons_sum_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))

        } else {

          # Throw a warning if selected codons exceed sequence length

          if(max(codon_pos) > max(lengths(seq_list)/3)) {

            stop("Selected codons exceed sequence length.")

          } else {

            # create a vector for counting loop iterations
            iter <- numeric(0)

            # Identify codons with synonymous nucleotide substitutions

            for (i in 1:(length(colnames(seq_file))-1)) {

              for (j in (i+1):length(colnames(seq_file))) {

                iter <- sum(iter, 1)

                # Identify codons with no amino acid change
                InvCod <- which(seq_list_aa[[i]] == seq_list_aa[[j]])

                # Identify codons with nucleotide differences
                VarCod <- vector()

                for(k in 1:(length(seq_list[[i]])/3)) {

                  if(toupper(paste(seq_list[[i]][c(k*3-2,k*3-1,k*3)],collapse="")) != toupper(paste(seq_list[[j]][c(k*3-2,k*3-1,k*3)],collapse=""))) {VarCod <- append(VarCod, k)}

                }

                # Identify codons with synonymous nucleotide substitutions
                SynCod <- VarCod[which(VarCod %in% InvCod)]

                # Reduce SynCod to only include codons that are found in codon_pos
                SynCod <- SynCod[which(SynCod %in% codon_pos)]

                # Get the base positions that correspond to the codons in SynCod
                SynCod_bases <- vector()

                if(length(SynCod) > 0) {

                  for(k in 1:(length(SynCod))) {

                    SynCod_bases <- append(SynCod_bases,c(SynCod[k]*3-2,SynCod[k]*3-1,SynCod[k]*3))

                  }

                }

                # Identify the positions harboring nucleotide changes in SynCod_bases
                VarSynCod_bases <- vector()

                if(length(SynCod) > 0) {

                  for(k in 1:(length(SynCod_bases))) {

                    if(toupper(seq_list[[i]][SynCod_bases[k]]) != toupper(seq_list[[j]][SynCod_bases[k]])) {VarSynCod_bases <- append(VarSynCod_bases, SynCod_bases[k])}

                  }

                }

                # Add presence/absence of synonymous nucleotide substitutions to
                # table with all pairwise comparisons
                # add 1 in the positions with synonymous nucleotide substitutions
                SynSubTab[iter,VarSynCod_bases] <- 1

                # similar for the codons
                SynCodTab[iter,SynCod] <- 1

              }

            }

            # Output a vector with column sums for each table
            # Include nucleotide triplet positions in the
            # table with synonymous nucleotide substitutions
            # Output also a vector with the proportions of
            # pairwise comparisons that harbored a synonymous
            # substitution in each base position or codon
            SynSub_sum <- cbind.data.frame(TripPos,colSums(SynSubTab),colSums(SynSubTab)/n_pairs)
            rownames(SynSub_sum) <- paste0("bp",1:n_nucl)
            colnames(SynSub_sum) <- c("Triplet.Pos","No.Syn.Obs","Prop.Syn.Obs")
            SynCod_sum <- cbind.data.frame(colSums(SynCodTab),colSums(SynCodTab)/n_pairs)
            rownames(SynCod_sum) <- paste0("cod",1:n_cod)
            colnames(SynCod_sum) <- c("No.Syn.Obs","Prop.Syn.Obs")

            ## Adjust output to only report specified codon_pos

            # Convert codon_pos to nucleotide positions
            codon_pos_bases <- vector()

            for(k in 1:(length(codon_pos))) {

              codon_pos_bases <- append(codon_pos_bases,c(codon_pos[k]*3-2,codon_pos[k]*3-1,codon_pos[k]*3))

            }

            # Fill NA in the bases that are not specified in codon_pos_bases
            for(i in 1:n_nucl) {

              if(!(i %in% codon_pos_bases)) {

                SynSub_sum$No.Syn.Obs[i] <- NA
                SynSub_sum$Prop.Syn.Obs[i] <- NA

              }

            }

            # Fill NA in the codons that are not specified in codon_pos
            for(i in 1:n_cod) {

              if(!(i %in% codon_pos)) {

                SynCod_sum$No.Syn.Obs[i] <- NA
                SynCod_sum$Prop.Syn.Obs[i] <- NA

              }

            }

            # Save the vectors as .csv
            write.csv(SynSub_sum,file=paste0(path_out,"/Syn_subst_sum_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))
            write.csv(SynCod_sum,file=paste0(path_out,"/Syn_codons_sum_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))

          }

        }

      }

      ##############################
      ### p-distance calculation ###
      ##############################

      ### calculate pairwise distances between all sequences in the
      ### sequence table, and save the values in the dist_matrix
      ### (p-distance values rounded to five digits)

      if(analysis == "dist") {

        if(is.null(codon_pos)) {

          # Calculate the synonymous nucleotide Pdist in each pairwise comparison of the sequences
          # in seq_list

          for (i in 1:(length(colnames(seq_file))-1)) {

            for (j in (i+1):length(colnames(seq_file))) {

              # Identify codons with no amino acid change
              InvCod <- which(seq_list_aa[[i]] == seq_list_aa[[j]])

              # Identify codons with nucleotide differences
              VarCod <- vector()

              for(k in 1:(length(seq_list[[i]])/3)) {

                if(toupper(paste(seq_list[[i]][c(k*3-2,k*3-1,k*3)],collapse="")) != toupper(paste(seq_list[[j]][c(k*3-2,k*3-1,k*3)],collapse=""))) {VarCod <- append(VarCod, k)}

              }

              # Identify codons with synonymous nucleotide substitutions
              SynCod <- VarCod[which(VarCod %in% InvCod)]

              # Get the base positions that correspond to the codons in SynCod
              SynCod_bases <- vector()

              for(k in 1:(length(SynCod))) {

                SynCod_bases <- append(SynCod_bases,c(SynCod[k]*3-2,SynCod[k]*3-1,SynCod[k]*3))

              }

              # Save the number of synonymous nucleotide changes in the upper right part of the dist_matrix
              dist_matrix[i,j] <- length(which(seq_list[[i]][SynCod_bases] != seq_list[[j]][SynCod_bases]))

              # Calculate P-distance as the number of synonymous nucleotide changes divided by total sequence length
              # Save it in the lower left part of the dist_matrix
              dist_matrix[j,i] <- round(length(which(seq_list[[i]][SynCod_bases] != seq_list[[j]][SynCod_bases]))/length(seq_list[[i]]), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed sequence length

          if(max(codon_pos) > max(lengths(seq_list)/3)) {

            stop("Selected codons exceed sequence length.")

          } else {

            # Calculate the synonymous nucleotide Pdist in each pairwise comparison of the sequences
            # in seq_list, using the vector codon_pos to define which codons to compare

            for (i in 1:(length(colnames(seq_file))-1)) {

              for (j in (i+1):length(colnames(seq_file))) {

                # Identify codons with no amino acid change
                InvCod <- which(seq_list_aa[[i]] == seq_list_aa[[j]])

                # Identify codons with nucleotide differences
                VarCod <- vector()

                for(k in 1:(length(seq_list[[i]])/3)) {

                  if(toupper(paste(seq_list[[i]][c(k*3-2,k*3-1,k*3)],collapse="")) != toupper(paste(seq_list[[j]][c(k*3-2,k*3-1,k*3)],collapse=""))) {VarCod <- append(VarCod, k)}

                }

                # Identify codons with synonymous nucleotide substitutions
                SynCod <- VarCod[which(VarCod %in% InvCod)]

                # Reduce SynCod to only include codons that are found in codon_pos
                SynCod <- SynCod[which(SynCod %in% codon_pos)]

                # Get the base positions that correspond to the codons in SynCod
                SynCod_bases <- vector()

                for(k in 1:(length(SynCod))) {

                  SynCod_bases <- append(SynCod_bases,c(SynCod[k]*3-2,SynCod[k]*3-1,SynCod[k]*3))

                }

                # Save the number of synonymous nucleotide changes among the codons in codon_pos in the upper right part of the dist_matrix
                dist_matrix[i,j] <- length(which(seq_list[[i]][SynCod_bases] != seq_list[[j]][SynCod_bases]))

                # Calculate P-distance as the number of synonymous nucleotide changes divided by total number of bases in the codons in codon_pos
                # Save it in the lower left part of the dist_matrix
                dist_matrix[j,i] <- round(length(which(seq_list[[i]][SynCod_bases] != seq_list[[j]][SynCod_bases]))/(3*length(codon_pos)), digits=5)

              }

            }

          }

        }

        # Export the distance matrix as a .csv file
        write.csv(dist_matrix,file=paste0(path_out,"/Syndist_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))

        ### Calculate mean distances for each sample in the data set

        for (i in 1:length(sample_names)) {

          # Fetch column numbers for the sequences in sample i in the sequence
          # table

          z <- which(seq_file[i,] > 0)

          if(length(z) < 2) {

            mean_diff_bases[i] <- NA
            mean_pdist[i] <- NA

          } else {

            # Create two vectors md_bases and md_pdist

            md_bases <- vector("numeric", length=length(z))
            md_pdist <- vector("numeric", length=length(z))

            # Generate a list of all pairwise combinations of the elements in z

            pwc <- combn(sort(z),2,simplify=T)

            # For each combination in pwc

            for(j in 1:length(pwc[1,]))  {

              # Extract the distance in each pairwise comparison of the sequences
              # in seq_list from the dist_matrix, using the numbers from pwc as
              # indices to extract the values from the matrix

              # number of synonymous bases is extracted from the upper right
              md_bases[j] <- dist_matrix[pwc[1,j],pwc[2,j]]

              # synonymous p-distance is extracted from the lower left
              md_pdist[j] <- dist_matrix[pwc[2,j],pwc[1,j]]

            }

            # Add the means of the md vectors to the mean_diff_bases and mean_pdist vectors

            mean_diff_bases[i] <- mean(md_bases)
            mean_pdist[i] <- mean(md_pdist)

          }

        }

        # Export a table with sample names and mean distances as a .csv file

        mean_diff_bases_tab <- as.data.frame(mean_diff_bases)
        rownames(mean_diff_bases_tab) <- sample_names
        write.csv(mean_diff_bases_tab,file=paste0(path_out,"/Mean_no_syn_bases_table_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))

        # Export a table with sample names and mean distances as a .csv file

        mean_pdist_tab <- as.data.frame(mean_pdist)
        rownames(mean_pdist_tab) <- sample_names
        write.csv(mean_pdist_tab,file=paste0(path_out,"/Mean_syn_pdist_table_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))

      }

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

      for(i in 1:length(seq_file)) {

        for(j in 1:length(unique(seq_file[[i]]))) {

          if(!(toupper(unique(seq_file[[i]])[j]) %in% c("-","A","C","G","T"))) {

            stop("Sequences contain non-standard characters. Accepted characters are -atgcATGC.")

          }

        }

      }

      ### Create a vector seq_list_aa containing all the sequences in the fasta file translated to amino acids

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

      for(i in 1:length(seq_file)) {

        # fill in NaNs in seq_list_aa
        seq_list_aa[[i]] <- NaN*seq(length(seq_file[[i]])/3)

        # check whether sequence lengths are divisible by 3
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

          stop("Sequences should be aligned in open reading frame 5' -> 3'")

        }

      }

      #####################################################
      ### Identify codons with synonymous substitutions ###
      #####################################################

      ### Option to output tables with number of pairwise comparisons
      ### with synonymous nucleotide changes per base and per codon in
      ### the data set

      if(analysis == "codon") {

        # Create table for presence/absence of synonymous nucleotide
        # substitutions in all pairwise comparisons
        n_pairs <- length(seq_file)*(length(seq_file)-1)/2 # number of pairwise sequence comparisons
        n_nucl <- max(lengths(seq_file)) # sequence length in nucleotides
        SynSubTab <- matrix(data=rep(0,n_pairs*n_nucl), nrow=n_pairs, ncol=n_nucl)

        # Create table for presence/absence of synonymous codons
        # in all pairwise comparisons
        n_cod <- max(lengths(seq_file)/3) # sequence length in codons
        SynCodTab <- matrix(data=rep(0,n_pairs*n_cod), nrow=n_pairs, ncol=n_cod)

        # Create two vectors to sum the presence/absence of synonymous nucleotide
        # substitutions or codons in all pairwise comparisons
        SynSub_sum <- vector()
        SynCod_sum <- vector()
        # Create a vector that defines nucleotide triplet positions
        TripPos <- as.character(rep(c(1,2,3),n_nucl/3))

        if(is.null(codon_pos)) {

          # create a vector for counting loop iterations
          iter <- numeric(0)

          # Identify codons with synonymous nucleotide substitutions

          for (i in 1:(length(seq_file)-1)) {

            for (j in (i+1):length(seq_file)) {

              iter <- sum(iter, 1)

              # Identify codons with no amino acid change
              InvCod <- which(seq_list_aa[[i]] == seq_list_aa[[j]])

              # Identify codons with nucleotide differences
              VarCod <- vector()

              for(k in 1:(length(seq_file[[i]])/3)) {

                if(toupper(paste(seq_file[[i]][c(k*3-2,k*3-1,k*3)],collapse="")) != toupper(paste(seq_file[[j]][c(k*3-2,k*3-1,k*3)],collapse=""))) {VarCod <- append(VarCod, k)}

              }

              # Identify codons with synonymous nucleotide substitutions
              SynCod <- VarCod[which(VarCod %in% InvCod)]

              # Get the base positions that correspond to the codons in SynCod
              SynCod_bases <- vector()

              if(length(SynCod) > 0) {

                for(k in 1:(length(SynCod))) {

                  SynCod_bases <- append(SynCod_bases,c(SynCod[k]*3-2,SynCod[k]*3-1,SynCod[k]*3))

                }

              }

              # Identify the positions harboring nucleotide changes in SynCod_bases
              VarSynCod_bases <- vector()

              if(length(SynCod) > 0) {

                for(k in 1:(length(SynCod_bases))) {

                  if(toupper(seq_file[[i]][SynCod_bases[k]]) != toupper(seq_file[[j]][SynCod_bases[k]])) {VarSynCod_bases <- append(VarSynCod_bases, SynCod_bases[k])}

                }

              }

              # Add presence/absence of synonymous nucleotide substitutions to
              # table with all pairwise comparisons
              # add 1 in the positions with synonymous nucleotide substitutions
              SynSubTab[iter,VarSynCod_bases] <- 1

              # similar for the codons
              SynCodTab[iter,SynCod] <- 1

            }

          }

          # Output a vector with column sums for each table
          # Include nucleotide triplet positions in the
          # table with synonymous nucleotide substitutions
          # Output also a vector with the proportions of
          # pairwise comparisons that harbored a synonymous
          # substitution in each base position or codon
          SynSub_sum <- cbind.data.frame(TripPos,colSums(SynSubTab),colSums(SynSubTab)/n_pairs)
          rownames(SynSub_sum) <- paste0("bp",1:n_nucl)
          colnames(SynSub_sum) <- c("Triplet.Pos","No.Syn.Obs","Prop.Syn.Obs")
          SynCod_sum <- cbind.data.frame(colSums(SynCodTab),colSums(SynCodTab)/n_pairs)
          rownames(SynCod_sum) <- paste0("cod",1:n_cod)
          colnames(SynCod_sum) <- c("No.Syn.Obs","Prop.Syn.Obs")

          # Save the vectors as .csv
          write.csv(SynSub_sum,file=paste0(path_out,"/Syn_subst_sum_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))
          write.csv(SynCod_sum,file=paste0(path_out,"/Syn_codons_sum_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))

        } else {

          # Throw a warning if selected codons exceed sequence length

          if(max(codon_pos) > max(lengths(seq_file)/3)) {

            stop("Selected codons exceed sequence length.")

          } else {

            # create a vector for counting loop iterations
            iter <- numeric(0)

            # Identify codons with synonymous nucleotide substitutions

            for (i in 1:(length(seq_file)-1)) {

              for (j in (i+1):length(seq_file)) {

                iter <- sum(iter, 1)

                # Identify codons with no amino acid change
                InvCod <- which(seq_list_aa[[i]] == seq_list_aa[[j]])

                # Identify codons with nucleotide differences
                VarCod <- vector()

                for(k in 1:(length(seq_file[[i]])/3)) {

                  if(toupper(paste(seq_file[[i]][c(k*3-2,k*3-1,k*3)],collapse="")) != toupper(paste(seq_file[[j]][c(k*3-2,k*3-1,k*3)],collapse=""))) {VarCod <- append(VarCod, k)}

                }

                # Identify codons with synonymous nucleotide substitutions
                SynCod <- VarCod[which(VarCod %in% InvCod)]

                # Reduce SynCod to only include codons that are found in codon_pos
                SynCod <- SynCod[which(SynCod %in% codon_pos)]

                # Get the base positions that correspond to the codons in SynCod
                SynCod_bases <- vector()

                if(length(SynCod) > 0) {

                  for(k in 1:(length(SynCod))) {

                    SynCod_bases <- append(SynCod_bases,c(SynCod[k]*3-2,SynCod[k]*3-1,SynCod[k]*3))

                  }

                }

                # Identify the positions harboring nucleotide changes in SynCod_bases
                VarSynCod_bases <- vector()

                if(length(SynCod) > 0) {

                  for(k in 1:(length(SynCod_bases))) {

                    if(toupper(seq_file[[i]][SynCod_bases[k]]) != toupper(seq_file[[j]][SynCod_bases[k]])) {VarSynCod_bases <- append(VarSynCod_bases, SynCod_bases[k])}

                  }

                }

                # Add presence/absence of synonymous nucleotide substitutions to
                # table with all pairwise comparisons
                # add 1 in the positions with synonymous nucleotide substitutions
                SynSubTab[iter,VarSynCod_bases] <- 1

                # similar for the codons
                SynCodTab[iter,SynCod] <- 1

              }

            }

            # Output a vector with column sums for each table
            # Include nucleotide triplet positions in the
            # table with synonymous nucleotide substitutions
            # Output also a vector with the proportions of
            # pairwise comparisons that harbored a synonymous
            # substitution in each base position or codon
            SynSub_sum <- cbind.data.frame(TripPos,colSums(SynSubTab),colSums(SynSubTab)/n_pairs)
            rownames(SynSub_sum) <- paste0("bp",1:n_nucl)
            colnames(SynSub_sum) <- c("Triplet.Pos","No.Syn.Obs","Prop.Syn.Obs")
            SynCod_sum <- cbind.data.frame(colSums(SynCodTab),colSums(SynCodTab)/n_pairs)
            rownames(SynCod_sum) <- paste0("cod",1:n_cod)
            colnames(SynCod_sum) <- c("No.Syn.Obs","Prop.Syn.Obs")

            ## Adjust output to only report specified codon_pos

            # Convert codon_pos to nucleotide positions
            codon_pos_bases <- vector()

            for(k in 1:(length(codon_pos))) {

              codon_pos_bases <- append(codon_pos_bases,c(codon_pos[k]*3-2,codon_pos[k]*3-1,codon_pos[k]*3))

            }

            # Fill NA in the bases that are not specified in codon_pos_bases
            for(i in 1:n_nucl) {

              if(!(i %in% codon_pos_bases)) {

                SynSub_sum$No.Syn.Obs[i] <- NA
                SynSub_sum$Prop.Syn.Obs[i] <- NA

              }

            }

            # Fill NA in the codons that are not specified in codon_pos
            for(i in 1:n_cod) {

              if(!(i %in% codon_pos)) {

                SynCod_sum$No.Syn.Obs[i] <- NA
                SynCod_sum$Prop.Syn.Obs[i] <- NA

              }

            }

            # Save the vectors as .csv
            write.csv(SynSub_sum,file=paste0(path_out,"/Syn_subst_sum_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))
            write.csv(SynCod_sum,file=paste0(path_out,"/Syn_codons_sum_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))

          }

        }

      }

      ##############################
      ### p-distance calculation ###
      ##############################

      ### calculate pairwise distances between all sequences in the
      ### fasta file, and save the values in the dist_matrix
      ### (p-distance values rounded to five digits)

      if(analysis == "dist") {

        if(is.null(codon_pos)) {

          # Calculate the synonymous nucleotide Pdist in each pairwise comparison of the sequences
          # in seq_list

          for (i in 1:(length(seq_file)-1)) {

            for (j in (i+1):length(seq_file)) {

              # Identify codons with no amino acid change
              InvCod <- which(seq_list_aa[[i]] == seq_list_aa[[j]])

              # Identify codons with nucleotide differences
              VarCod <- vector()

              for(k in 1:(length(seq_file[[i]])/3)) {

                if(toupper(paste(seq_file[[i]][c(k*3-2,k*3-1,k*3)],collapse="")) != toupper(paste(seq_file[[j]][c(k*3-2,k*3-1,k*3)],collapse=""))) {VarCod <- append(VarCod, k)}

              }

              # Identify codons with synonymous nucleotide substitutions
              SynCod <- VarCod[which(VarCod %in% InvCod)]

              # Get the base positions that correspond to the codons in SynCod
              SynCod_bases <- vector()

              for(k in 1:(length(SynCod))) {

                SynCod_bases <- append(SynCod_bases,c(SynCod[k]*3-2,SynCod[k]*3-1,SynCod[k]*3))

              }

              # Save the number of synonymous nucleotide changes in the upper right part of the dist_matrix
              dist_matrix[i,j] <- length(which(seq_file[[i]][SynCod_bases] != seq_file[[j]][SynCod_bases]))

              # Calculate P-distance as the number of synonymous nucleotide changes divided by total sequence length
              # Save it in the lower left part of the dist_matrix
              dist_matrix[j,i] <- round(length(which(seq_file[[i]][SynCod_bases] != seq_file[[j]][SynCod_bases]))/length(seq_file[[i]]), digits=5)

            }

          }

        } else {

          # Throw a warning if selected codons exceed sequence length

          if(max(codon_pos) > max(lengths(seq_file)/3)) {

            stop("Selected codons exceed sequence length.")

          } else {

            # Calculate the synonymous nucleotide Pdist in each pairwise comparison of the sequences
            # in seq_file, using the vector codon_pos to define which codons to compare

            for (i in 1:(length(seq_file)-1)) {

              for (j in (i+1):length(seq_file)) {

                # Identify codons with no amino acid change
                InvCod <- which(seq_list_aa[[i]] == seq_list_aa[[j]])

                # Identify codons with nucleotide differences
                VarCod <- vector()

                for(k in 1:(length(seq_file[[i]])/3)) {

                  if(toupper(paste(seq_file[[i]][c(k*3-2,k*3-1,k*3)],collapse="")) != toupper(paste(seq_file[[j]][c(k*3-2,k*3-1,k*3)],collapse=""))) {VarCod <- append(VarCod, k)}

                }

                # Identify codons with synonymous nucleotide substitutions
                SynCod <- VarCod[which(VarCod %in% InvCod)]

                # Reduce SynCod to only include codons that are found in codon_pos
                SynCod <- SynCod[which(SynCod %in% codon_pos)]

                # Get the base positions that correspond to the codons in SynCod
                SynCod_bases <- vector()

                for(k in 1:(length(SynCod))) {

                  SynCod_bases <- append(SynCod_bases,c(SynCod[k]*3-2,SynCod[k]*3-1,SynCod[k]*3))

                }

                # Save the number of synonymous nucleotide changes among the codons in codon_pos in the upper right part of the dist_matrix
                dist_matrix[i,j] <- length(which(seq_file[[i]][SynCod_bases] != seq_file[[j]][SynCod_bases]))

                # Calculate P-distance as the number of synonymous nucleotide changes divided by total number of bases in the codons in codon_pos
                # Save it in the lower left part of the dist_matrix
                dist_matrix[j,i] <- round(length(which(seq_file[[i]][SynCod_bases] != seq_file[[j]][SynCod_bases]))/(3*length(codon_pos)), digits=5)

              }

            }

          }

        }

        # Export the distance matrix as a .csv file
        write.csv(dist_matrix,file=paste0(path_out,"/Syndist_matrix_", c(format(Sys.Date(),"%Y%m%d")), ".csv"))

      }

    }

}
