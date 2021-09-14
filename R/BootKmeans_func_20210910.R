#' BootKmeans() function
#'
#' \code{\link{BootKmeans}} is a wrapper for the kmeans() function of the 'stats'
#' package, which allows for bootstrapping. Bootstrapping k-estimates may be
#' desirable in data sets, where the BIC- vs. k-values do not produce clear
#' inflection points ("elbows").
#'
#' BootKmeans() performs multiple runs of kmeans() scanning k-values from 1 to a
#' maximum value defined by the user. In each scan, an optimal k-value is
#' estimated using a user-defined threshold of BIC reduction. The method is an
#' automated version of visually inspecting elbow plots of BIC- vs. k-values.
#' The number of scans to be performed is defined by the user.
#'
#' For each k-estimate scan, the algorithm produces a summary of the stats incl.
#' total within SS, AIC, and BIC, an elbow plot (BIC vs. k), and a set of cluster
#' files corresponding to the estimated optimal k-value. It also produces a table
#' summarizing the stats of the final selected kmeans() models corresponding to
#' the estimated optimal k-values of each scan.
#'
#' After running BootKmeans() on a data set, it is recommended to subsequently
#' evaluate the repeatability of the bootstrapped k-estimation scans with the
#' ClusterMatch() function also included in MHCtools.
#'
#' Input data format:
#' A set of five z-matrices containing numerical values of the z-descriptors
#' (z1-z5) for each amino acid position in a sequence alignment. Each column
#' should represent an amino acid position and each row one sequence in the
#' alignment.
#'
#' If you publish data produced with MHCtools, please cite:
#' Roved, J. 2020. MHCtools: Analysis of MHC data in non-model species. Cran.
#' Roved, J., Hansson, B., Stervander, M., Hasselquist, D., Westerdahl, H. 2020.
#' Non-random association of MHC-I alleles in favor of high diversity haplotypes
#' in wild songbirds revealed by computer-assisted MHC haplotype inference using
#' the R package MHCtools. bioRxiv.
#'
#' @param z1_matrix a matrix with numerical values of the first z-descriptor for
#'   each amino acid position in all sequences in the data set.
#' @param z2_matrix a matrix with numerical values of the second z-descriptor
#'   for each amino acid position in all sequences in the data set.
#' @param z3_matrix a matrix with numerical values of the third z-descriptor for
#'   each amino acid position in all sequences in the data set.
#' @param z4_matrix a matrix with numerical values of the fourth z-descriptor
#'   for each amino acid position in all sequences in the data set.
#' @param z5_matrix a matrix with numerical values of the fifth z-descriptor for
#'   each amino acid position in all sequences in the data set.
#' @param threshold a numerical value between 0 and 1 specifying the threshold
#'   of reduction in BIC for selecting a k estimate for each kmeans clustering
#'   model. The value specifies a proportion of the max observed reduction in BIC
#'   when increasing k by 1 (default 0.01).
#' @param no_scans an integer specifying the number of k estimation scans
#'   to run (default 1,000).
#' @param max_k an integer specifying the hypothetical maximum number of clusters
#'   to detect (default 40). In each k estimation scan, the algorithm runs a
#'   kmeans() clustering model for each value of k between 1 and max_k.
#' @param iter.max an integer specifying the maximum number of iterations allowed
#'   in each kmeans() clustering model (default 1,000,000).
#' @param nstart an integer specifying the number of rows in the set of input
#'   matrices that will be chosen as initial centers in the kmeans() clustering
#'   models (default 200).
#' @param algorithm character vector, specifying the method for the kmeans()
#'   clustering function, one of c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
#'   default is "Hartigan-Wong".
#' @param path_out a user defined path to the folder where the output files
#'   will be saved.
#' @return The function produces three folders in path_out, which contain for
#'   each scan the estimated k-clusters saved as .Rdata files, an elbow plot saved
#'   as .pdf, and a stats summary table saved as a .csv file. In path_out a summary
#'   of all scans performed in the bootstrap run is also saved as .csv. This table
#'   is also shown in the console.
#'   Should alternative elbow plots be desired, they may be produced manually with
#'   the stats presented in the summary tables for each scan.
#' @note AIC and BIC are calculated from the kmeans model objects by the following
#'   formulae:
#'   - AIC = D + 2*m*k
#'   - BIC = D + log(n)*m*k
#'   in which:
#'   - m = ncol(fit$centers)
#'   - n = length(fit$cluster)
#'   - k = nrow(fit$centers)
#'   - D = fit$tot.withinss
#' @seealso \code{\link{ClusterMatch}}; \code{\link{DistCalc}}
#' @examples
#' z1_matrix <- z1_matrix
#' z2_matrix <- z2_matrix
#' z3_matrix <- z3_matrix
#' z4_matrix <- z4_matrix
#' z5_matrix <- z5_matrix
#' path_out <- tempdir()
#' BootKmeans(z1_matrix, z2_matrix, z3_matrix, z4_matrix, z5_matrix, threshold=0.01,
#' no_scans=10, max_k=20, iter.max=10, nstart=10, algorithm="Hartigan-Wong",
#' path_out=path_out)
#' @importFrom "stats" "kmeans"
#' @importFrom "grDevices" "dev.off" "pdf"
#' @importFrom "graphics" "abline" "mtext" "par" "plot"
#' @export

BootKmeans <- function(z1_matrix, z2_matrix, z3_matrix, z4_matrix, z5_matrix, threshold=0.01, no_scans=1000, max_k=40, iter.max=1000000, nstart=200, algorithm="Hartigan-Wong", path_out=path_out) {

  # Generate a matrix that combines the five z-descriptor input matrices
  z_matrix <- matrix(NA,nrow=length(rownames(z1_matrix)),ncol=5*length(colnames(z1_matrix)))
  rownames(z_matrix) <- rownames(z1_matrix)
  colnames(z_matrix) <- vector(length=5*length(colnames(z1_matrix)))

  for(i in 1:length(colnames(z1_matrix))) {

    for(j in 1:5) {

      colnames(z_matrix)[(i-1)*5+j] <- paste("c", i, "_z", j, sep="")

    }

  }

  for(i in 1:length(colnames(z1_matrix))) {

    z_matrix[,(i-1)*5+1] <- z1_matrix[,i]
    z_matrix[,(i-1)*5+2] <- z2_matrix[,i]
    z_matrix[,(i-1)*5+3] <- z3_matrix[,i]
    z_matrix[,(i-1)*5+4] <- z4_matrix[,i]
    z_matrix[,(i-1)*5+5] <- z5_matrix[,i]

  }

  # Create folders in the output path that will hold the K-clusters, K-stats tables, and elbow plots for each scan
  dir.create(paste(path_out,"/Clusters",sep=""))
  dir.create(paste(path_out,"/Kstats_tables",sep=""))
  dir.create(paste(path_out,"/plots",sep=""))

  # Create vectors that will hold k estimates and residual sums of squares, AIC, and BIC values for each scan
  k.est <- numeric(length=no_scans)
  Totss.resid <- numeric(length=no_scans)
  Tot.withinss.resid <- numeric(length=no_scans)
  Betweenss.resid <- numeric(length=no_scans)
  AIC.resid <- numeric(length=no_scans)
  BIC.resid <- numeric(length=no_scans)
  BIC.max <- numeric(length=no_scans)
  BIC.min <- numeric(length=no_scans)

  # Run clustering models using the kmeans() function from the 'stats' package
  # estimate the number of clusters from the incremental change in BIC
  for(i in 1:no_scans) {

    # create a vector that will contain the k-means stats for each scan
    assign(paste("Kstats_",i,sep=""), vector())
    # create a temporary list that will contain the cluster assignments for each value of k
    for(k in 1:max_k) {assign(paste("temp_kclusters_",k,sep=""), list())}
    # create a list that will contain the final cluster assignments for each scan
    assign(paste("Kclusters_",i,sep=""), list())
    # assign stats and cluster assignments to the respective vector and list
    for(k in 1:max_k) {
      model <- kmeans(x = z_matrix, centers = k, iter.max, nstart, algorithm)
      assign(paste("temp_kclusters_",k,sep=""), model$cluster)
      assign(paste("Kstats_",i,sep=""), append(get(paste("Kstats_",i,sep="")), c(
        k,
        model$totss,
        model$tot.withinss,
        model$betweenss,
        model$tot.withinss + 2*ncol(model$centers)*nrow(model$centers),
        model$tot.withinss + log(length(model$cluster))*ncol(model$centers)*nrow(model$centers))))
    }
    # reformat the Kstats vector to a data frame
    assign(paste("Kstats_",i,sep=""), data.frame(matrix(get(paste("Kstats_",i,sep="")), nrow=length(get(paste("Kstats_",i,sep="")))/6, byrow=T, dimnames=list(c(paste("k=",c(1:max_k),sep="")), c("k","Total.ss","Tot.within.ss","Between.ss","AIC","BIC")))))
    # estimate the number of clusters as the value of k, where the reduction in BIC is below the desired threshold of the maximal reduction observed
    # note, using min() because the values of diff(BIC) are negative
    k.est[i] <- which(diff(get(paste("Kstats_",i,sep=""))$BIC) > threshold*min(diff(get(paste("Kstats_",i,sep=""))$BIC)))[1]
    # extract the final cluster assignments for the estimated value of k and save as a list
    assign(paste("Kclusters_",i,sep=""), get(paste("temp_kclusters_",k.est[i],sep="")))
    saveRDS(get(paste("Kclusters_",i,sep="")), file=paste(path_out,"/Clusters","/Kclusters_model_",i,"_",c(format(Sys.Date(),"%Y%m%d")),".RData",sep=""))
    # extract residual sums of squares, AIC, and BIC values for each scan
    Totss.resid[i] <- get(paste("Kstats_",i,sep=""))$Total.ss[k.est[i]]
    Tot.withinss.resid[i] <- get(paste("Kstats_",i,sep=""))$Tot.within.ss[k.est[i]]
    Betweenss.resid[i] <- get(paste("Kstats_",i,sep=""))$Between.ss[k.est[i]]
    AIC.resid[i] <- get(paste("Kstats_",i,sep=""))$AIC[k.est[i]]
    BIC.resid[i] <- get(paste("Kstats_",i,sep=""))$BIC[k.est[i]]
    # extract min and max BIC values for each scan
    BIC.max[i] <- max(get(paste("Kstats_",i,sep=""))$BIC)
    BIC.min[i] <- min(get(paste("Kstats_",i,sep=""))$BIC)
    # save the Kstats table as .csv
    write.csv(get(paste("Kstats_",i,sep="")), file=paste(path_out,"/Kstats_tables","/Kstats_model_",i,"_",c(format(Sys.Date(),"%Y%m%d")),".csv",sep=""))
    # Create an elbow plot for the selected model in each scan and save as pdf
    pdf(paste(path_out,"/plots","/k_means_BIC_plot_model_",i,"_",c(format(Sys.Date(),"%Y%m%d")),".pdf",sep=""), width=5, height=5)
    par(mfrow=c(1,1))
    plot(get(paste("Kstats_",i,sep=""))$k, get(paste("Kstats_",i,sep=""))$BIC, xlab="k", ylab="Bayesian Information Criterion", type="b", main=paste("Model ",i,sep=""), bty="l")
    abline(v=k.est[i], col="red", lty=2)
    mtext(paste("k-est = ", k.est[i], sep=""), side=3, col="red")
    dev.off()

  }

  # Create a table with k estimates and residual, min, max, and delta BIC-values, and the residual sums of squares and AIC for selected model in each scan
  delta.BIC <- BIC.max-BIC.resid
  prop.delta.BIC <- delta.BIC/(BIC.max-BIC.min)
  delta.BIC.over.k <- delta.BIC/k.est
  k_summary <- cbind.data.frame(k.est,Totss.resid,Tot.withinss.resid,Betweenss.resid,AIC.resid,BIC.max,BIC.min,BIC.resid,delta.BIC,prop.delta.BIC,delta.BIC.over.k)
  rownames(k_summary) <- c(paste("model_",1:no_scans,sep=""))
  print(k_summary)

  # save the k_summary table as .csv
  write.csv(k_summary, file=paste(path_out,"/k_means_bootstrap_summary_stats_",c(format(Sys.Date(),"%Y%m%d")),".csv",sep=""))

}