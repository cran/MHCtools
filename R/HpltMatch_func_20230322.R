#' HpltMatch() function
#'
#' Putative haplotypes may be identical to each other, or they may differ only by
#' incongruent or unresolved sequences. It is therefore useful to curate putative
#' haplotypes by comparing them to identify potentially overlapping types as
#' candidates for further investigation. \code{\link{HpltMatch}} calculates the
#' proportion of matching sequences between pairs of haplotypes and produces a
#' .csv table with values in a lower left matrix. If a threshold value is
#' specified, a list of haplotype matches where the proportion of matching
#' sequences exceeds the threshold will be produced.
#'
#' Note: The NestTablesXL() function provides a useful format for further
#' investigation of potentially overlapping haplotypes.
#'
#' If you publish data or results produced with MHCtools, please cite both of
#' the following references:
#' Roved, J. (2022). MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J. (2024). MHCtools 1.5: Analysis of MHC sequencing data in R. In S.
#' Boegel (Ed.), HLA Typing: Methods and Protocols (2nd ed., pp. 275–295).
#' Humana Press. https://doi.org/10.1007/978-1-0716-3874-3_18
#'
#' @param hplt_occ_matrix A binary (logical) occurrence matrix with the data set
#'   sequences in columns and the putative haplotypes in rows, as produced by the
#'   CreateHpltOccTable() function.
#' @param path_out a user defined path to the folder where the output file(s) will
#'   be saved.
#' @param threshold a numerical value between 0 and 1 (default NULL) specifying a
#'   threshold for the proportion of matching sequences between haplotypes.
#' @return  A table specifying the proportions of matching sequences between
#'   pairs of haplotypes (in a lower left matrix). If a threshold value is
#'   specified, a list of haplotype matches where the proportion of matching
#'   sequences exceeds the threshold will be printed to the console. The list will
#'   also be saved in the output path, and can be reopened in R e.g. using the
#'   readRDS() function in the base package.
#'   Note: HpltMatch() will overwrite any existing files with the same output
#'   file names in path_out.
#' @seealso \code{\link{HpltFind}}; \code{\link{CreateHpltOccTable}};
#'   \code{\link{NestTablesXL}}
#' @examples
#' hplt_occ_matrix <- hplt_occurrence_matrix
#' path_out <- tempdir()
#' HpltMatch(hplt_occ_matrix, path_out, threshold=NULL)
#' @export

HpltMatch <- function(hplt_occ_matrix, path_out, threshold=NULL) {

  # Put the names of the haplotypes into a vector called hplt_names
  hplt_names <- rownames(hplt_occ_matrix)

  # Create a matrix that will contain the distances between haplotypes
  hplt_dist_matrix <- matrix(nrow = length(hplt_names), ncol = length(hplt_names))
  colnames(hplt_dist_matrix) <- hplt_names
  rownames(hplt_dist_matrix) <- hplt_names

  # If a threshold value is specified, create a list that will contain vectors for each haplotype with the names of the putatively matching haplotypes
  if(!is.null(threshold)) {

    Matches <- vector("list", length(hplt_names))
    names(Matches) <- hplt_names

  }

  # Loop over all the haplotypes
  for (i in 1:(length(hplt_names)-1)) {

    # Fetch column numbers for the sequences in hplt i in the hplt_occ_matrix
    z <- which(hplt_occ_matrix[i,] > 0)

    # Match haplotype i against each of the others using a loop
    for (j in (i+1):length(hplt_names)) {

      # Fetch column numbers for the sequences in hplt j in the hplt_occ_matrix
      v <- which(hplt_occ_matrix[j,] > 0)

      # Create logical vectors with the congruence between hplt i and hplt j
      x <- z %in% v
      y <- v %in% z

      # Calculate the proportion of matches in comparisons between hplt i and hplt j
      # Values are entered in a lower left matrix in hplt_dist_matrix
      hplt_dist_matrix[j,i] <- (length(which(x==TRUE))+length(which(y==TRUE)))/(length(x)+length(y))

      # If a threshold value is specified, do a logical test of whether the proportion of matches between hplt i and hplt j > threshold
      if(!is.null(threshold)) {

        # If proportion of matches > threshold, append hplt j to Matches[[i]]
        Matches[[i]] <- factor()
        if(!is.na(hplt_dist_matrix[j,i]) & (hplt_dist_matrix[j,i] > threshold)) Matches[[i]] <- append(Matches[[i]], hplt_names[j])

      }

    }

  }

  # Save hplt_dist_matrix as .csv
  write.csv(hplt_dist_matrix, file=paste0(path_out,"/Hplt_dist_matrix_",c(format(Sys.Date(),"%Y%m%d")),".csv"))

  # If a threshold value is specified, print Matches and save it as .Rds
  if(!is.null(threshold)) {

    print(Matches)
    saveRDS(Matches, file=paste0(path_out,"/Hplt_matches_threshold_",threshold,"_",c(format(Sys.Date(),"%Y%m%d")),".Rds"))

  }

}
