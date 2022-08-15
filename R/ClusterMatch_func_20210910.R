#' ClusterMatch() function
#'
#' \code{\link{ClusterMatch}} is a tool for evaluating whether k-means()
#' clustering models with similar estimated values of k identify similar
#' clusters. ClusterMatch() also summarizes model stats as means for
#' different estimated values of k. It is designed to take files produced
#' by the BootKmeans() function as input, but other data can be analysed
#' if the descriptions of the data formats given below are observed
#' carefully.
#'
#' If you publish data or results produced with MHCtools, please cite both of
#' the following references:
#' Roved, J. 2022. MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J., Hansson, B., Stervander, M., Hasselquist, D., & Westerdahl, H. 2022.
#' MHCtools - an R package for MHC high-throughput sequencing data: genotyping,
#' haplotype and supertype inference, and downstream genetic analyses in non-model
#' organisms. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13645
#'
#' @param filepath a user defined path to a folder that contains the set of
#'   K-cluster files to be matched against each other. The algorithm will attempt
#'   to load all files in the folder, so it should contain only the relevant
#'   K-cluster files. If the clusters were generated using the BootKmeans()
#'   function, such a folder (named Clusters) was created by the algorithm in the
#'   output path given by the user.
#'   Each K-cluster file should correspond to the model$cluster object in kmeans()
#'   saved as a .Rdata file. Such files are generated as part of the output from
#'   BootKmeans(). ClusterMatch() assumes that the file names contain the string
#'   "model_" followed by a model number, which must match the corresponding row
#'   numbers in k_summary_table. If the data used was generated with the
#'   BootKmeans() function, the formats and numbers will match by default.
#' @param path_out a user defined path to the folder where the output files will
#'   be saved.
#' @param k_summary_table a data frame summarizing the stats of the kmeans()
#'   models that produced the clusters in the K-cluster files. If the data used
#'   was generated with the BootKmeans() function, a compatible
#'   k_summary_table was produced in the output path with the file name
#'   "k_means_bootstrap_summary_stats_<date>.csv".
#'   If other data is analysed, please observe these formatting requirements:
#'   The k_summary_table must contain the data for each kmeans() model in rows
#'   and as minimum the following columns:
#'   - k-value (colname: k.est)
#'   - residual total within sums-of-squares (colname: Tot.withinss.resid)
#'   - residual AIC (colname: AIC.resid)
#'   - residual BIC (colname: BIC.resid)
#'   - delta BIC/max BIC (colname: prop.delta.BIC)
#'   - delta BIC/k.est (colname: delta.BIC.over.k)
#'   It is crucial that the models have the same numbers in the K-cluster file
#'   names and in the k_summary_table, and that the rows of the table are ordered
#'   by the model number.
#' @return The function returns a summary table, which for each estimated number
#'   of clusters (i.e. the k-values of the models) lists:
#'   - number of models that found i clusters
#'   - mean residual total within sums-of-squares
#'   - mean residual AIC
#'   - mean residual BIC
#'   - mean delta BIC/max BIC
#'   - mean delta BIC/k
#'   - mean number of allele assignments that fall outside of the i most abundant
#'     clusters across all pairwise comparisons between the models that found i
#'     clusters
#'   - mean proportion of allele assignments that fall outside of the i most
#'     abundant clusters across all pairwise comparisons between the models that
#'     found i clusters
#'   The summary table is also saved as a .csv file in the output path.
#' @seealso \code{\link{BootKmeans}}
#' @examples
#' filepath <- system.file("extdata/ClusterMatch", package="MHCtools")
#' path_out <- tempdir()
#' k_summary_table <- k_summary_table
#' ClusterMatch(filepath, path_out, k_summary_table)
#' @importFrom "utils" "combn"
#' @importFrom "mgcv" "uniquecombs"
#' @export

ClusterMatch <- function(filepath, path_out, k_summary_table) {

  # Add a / to the end of the filepath
  filepath <- paste(filepath,"/",sep="")
  # Get the file names in the data folder
  file_names <- dir(filepath)
  # Sort the file names by model number
  file_names <- file_names[order(as.numeric(gsub("[^0-9]", "", regmatches(file_names, regexpr("model_[0-9]+", file_names)))))]
  # this sorting is done by extracting "model_[0-9]+" from the file name using a regex, then extracting the model number using gsub
  # Note: This assumes that the string "model_" followed by the model number appears in the file name.
  # The model numbers in the file names must match the corresponding ones specified in k_summary_table,
  # i.e. the models must have the same index numbers in the files and in k_summary_table.

  # Load the K-clusters
  for(i in 1:length(file_names)) {

    assign(paste("Kclusters_",i,sep=""), readRDS(file.path(filepath, file_names[i])))

  }

  # Create a vector of the observed k-estimates
  k_estimates <- sort(unique(k_summary_table$k.est))
  # Identify singletons
  k_est_singletons <- k_estimates[as.numeric(which(table(k_summary_table$k.est)==1))]

  # Create tables with the cluster assignments of the models that found each number of clusters
  for(i in 1:length(k_estimates)) {
    # create each vector tab<mink>-tab<maxk>
    assign(paste("tab",k_estimates[i],sep=""), vector())
    # assign $grp to vector
    for(j in 1:length(which(k_summary_table$k.est==k_estimates[i]))) {

      cluster_ass <- get(paste("Kclusters_",which(k_summary_table$k.est==k_estimates[i])[j],sep=""))
      assign(paste("tab",k_estimates[i],sep=""), append(get(paste("tab",k_estimates[i],sep="")), cluster_ass))

    }

    # reorganize vector to matrix
    assign(paste("tab",k_estimates[i],sep=""),matrix(get(paste("tab",k_estimates[i],sep="")), nrow=length(get(paste("tab",k_estimates[i],sep="")))/length(which(k_summary_table$k.est==k_estimates[i])), ncol=length(which(k_summary_table$k.est==k_estimates[i]))))

  }

  # Create a vector that will summarize the number of models for each value of k
  no_models <- numeric(length(k_estimates))
  # Create a vector that will summarize the mean number of allele assignments that fall outside of the k most abundant clusters accross all pairwise comparison between the models that found k clusters
  mean_no_ass_low_rank_clust <- numeric(length(k_estimates))
  # Create a vector that will summarize the mean proportion of allele assignments that fall outside of the k most abundant clusters accross all pairwise comparison between the models that found k clusters
  mean_prop_ass_low_rank_clust <- numeric(length(k_estimates))

  for(i in 1:length(k_estimates)) {

    # skip the singletons - combn() will fail on those
    if(k_estimates[i] %in% k_est_singletons) {

      mean_no_ass_low_rank_clust[i] <- NA
      mean_prop_ass_low_rank_clust[i] <- NA
      no_models[i] <- 1

    } else {

      # get all pairwise combinations of the columns of the table with the cluster assignments of the models that found k_estimates[i] clusters
      pwc <- combn(1:length(get(paste("tab",k_estimates[i],sep=""))[1,]), 2, simplify=F)

      # create a vector that will contain the number of allele assignments that fall outside of the k_estimates[i] most abundant clusters for each pairwise comparison between the models that found k_estimates[i] clusters
      no_ass_low_rank_clust <- vector(length=length(pwc))
      # create a vector that will contain the proportion of allele assignments that fall outside of the k_estimates[i] most abundant clusters for each pairwise comparison between the models that found k_estimates[i] clusters
      prop_ass_low_rank_clust <- vector(length=length(pwc))

      # Calculate the proportion of allele assignments that fall outside of the k_estimates[i] most abundant clusters for each pairwise comparison
      for(j in 1:length(pwc)) {

        # create a vector with the number of occurrences of each unique cluster in the jth pairwise comparison
        no_occ <- vector(length=max(attributes(uniquecombs(get(paste("tab",k_estimates[i],sep=""))[,pwc[[j]]]))$index))
        # the function uniquecombs() identifies unique rows in a table, which in the model tables represent unique clusters
        # the occurrences of each unique cluster are specified in the 'index' attribute of uniquecombs()
        # thus the number of occurrences of each unique cluster can be calculated from that

        for(l in 1:max(attributes(uniquecombs(get(paste("tab",k_estimates[i],sep=""))[,pwc[[j]]]))$index)) {

          no_occ[l] <- length(which(attributes(uniquecombs(get(paste("tab",k_estimates[i],sep=""))[,pwc[[j]]]))$index == l))

        }

        # sum up the total number of allele assignments that fall outside of the k_estimates[i] most abundant clusters in the jth pairwise comparison
        ifelse(length(no_occ) > k_estimates[i], no_ass_low_rank_clust[j] <- sum(sort(no_occ,decreasing=T)[(k_estimates[i]+1):length(no_occ)])*2, no_ass_low_rank_clust[j] <- 0)
        # Calculate the proportion of allele assignments that fall outside of the k_estimates[i] most abundant clusters in the jth pairwise comparison
        ifelse(length(no_occ) > k_estimates[i], prop_ass_low_rank_clust[j] <- sum(sort(no_occ,decreasing=T)[(k_estimates[i]+1):length(no_occ)])*2/length(get(paste("tab",k_estimates[i],sep=""))[,pwc[[j]]]), prop_ass_low_rank_clust[j] <- 0)

      }

      # calculate the mean number of allele assignments that fall outside of the k_estimates[i] most abundant clusters across all pairwise comparison between the models that found k_estimates[i] clusters
      mean_no_ass_low_rank_clust[i] <- mean(no_ass_low_rank_clust)
      # calculate the mean proportion of allele assignments that fall outside of the k_estimates[i] most abundant clusters across all pairwise comparison between the models that found k_estimates[i] clusters
      mean_prop_ass_low_rank_clust[i] <- mean(prop_ass_low_rank_clust)
      # assign the number of models for each value of k to the no_models vector
      no_models[i] <- length(which(k_summary_table$k.est==k_estimates[i]))

    }

  }

  # Evaluate the average residual total within group SS, AIC, and BIC in models with different estimated values of k
  mean.Tot.withinss.resid <- tapply(k_summary_table$Tot.withinss.resid,k_summary_table$k.est,mean)
  mean.AIC.resid <- tapply(k_summary_table$AIC.resid,k_summary_table$k.est,mean)
  mean.BIC.resid <- tapply(k_summary_table$BIC.resid,k_summary_table$k.est,mean)
  # Evaluate what proportion of the total BIC that is used on average by models with different estimated values of k
  mean.prop.delta.BIC <- tapply(k_summary_table$prop.delta.BIC,k_summary_table$k.est,mean)
  # Evaluate the average reduction in BIC per cluster in models with different estimated values of k
  mean.delta.BIC.over.k <- tapply(k_summary_table$delta.BIC.over.k,k_summary_table$k.est,mean)

  # Create a table that summarizes the repeatability of the allele assignment to clusters for each estimated value of k
  CM_summary <- cbind.data.frame(k_estimates, no_models, mean.Tot.withinss.resid, mean.AIC.resid, mean.BIC.resid, mean.prop.delta.BIC, mean.delta.BIC.over.k, mean_no_ass_low_rank_clust, mean_prop_ass_low_rank_clust)
  rownames(CM_summary) <- c(paste("k=",k_estimates,sep=""))
  colnames(CM_summary) <- c("k.est","No.models","Mean.Tot.withinss.resid","Mean.AIC.resid","Mean.BIC.resid","Mean.prop.delta.BIC","Mean.delta.BIC.over.k","Mean.NoAss.low.rank.clust","Mean.PropAss.low.rank.clust")
  print(CM_summary)

  # save the table as .csv
  write.csv(CM_summary, file=paste(path_out,"/ClusterMatch_summary_stats_",c(format(Sys.Date(),"%Y%m%d")),".csv",sep=""))

}
