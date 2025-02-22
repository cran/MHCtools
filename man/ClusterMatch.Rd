% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClusterMatch_func_20250108.R
\name{ClusterMatch}
\alias{ClusterMatch}
\title{ClusterMatch() function}
\usage{
ClusterMatch(filepath, path_out, k_summary_table)
}
\arguments{
\item{filepath}{a user defined path to a folder that contains the set of
K-cluster files to be matched against each other. The algorithm will attempt
to load all files in the folder, so it should contain only the relevant
K-cluster files. If the clusters were generated using the BootKmeans()
function, such a folder (named Clusters) was created by the algorithm in the
output path given by the user.
Each K-cluster file should correspond to the model$cluster object in kmeans()
saved as a .RData file. Such files are generated as part of the output from
BootKmeans(). ClusterMatch() assumes that the file names contain the string
"model_" followed by a model number, which must match the corresponding row
numbers in k_summary_table. If the data used was generated with the
BootKmeans() function, the formats and numbers will match by default.}

\item{path_out}{a user defined path to the folder where the output files will
be saved.}

\item{k_summary_table}{a data frame summarizing the stats of the kmeans()
models that produced the clusters in the K-cluster files. If the data used
was generated with the BootKmeans() function, a compatible
k_summary_table was produced in the output path with the file name
"k_means_bootstrap_summary_stats_<date>.csv".
If other data is analyzed, please observe these formatting requirements:
The k_summary_table must contain the data for each kmeans() model in rows
and as minimum the following columns:
- k-value (colname: k.est)
- residual total within sums-of-squares (colname: Tot.withinss.resid)
- residual AIC (colname: AIC.resid)
- residual BIC (colname: BIC.resid)
- delta BIC/max BIC (colname: prop.delta.BIC)
- delta BIC/k.est (colname: delta.BIC.over.k)
It is crucial that the models have the same numbers in the K-cluster file
names and in the k_summary_table, and that the rows of the table are ordered
by the model number.}
}
\value{
The function returns a summary table, which for each estimated number
  of clusters (i.e. the k-values of the models) lists:
  - number of models that found i clusters
  - mean residual total within sums-of-squares
  - mean residual AIC
  - mean residual BIC
  - mean delta BIC/max BIC
  - mean delta BIC/k
  - mean number of allele assignments that fall outside of the i most abundant
    clusters across all pairwise comparisons between the models that found i
    clusters
  - mean proportion of allele assignments that fall outside of the i most
    abundant clusters across all pairwise comparisons between the models that
    found i clusters
  The summary table is also saved as a .csv file in the output path.
}
\description{
\code{\link{ClusterMatch}} is a tool for evaluating whether k-means()
clustering models with similar estimated values of k identify similar
clusters. ClusterMatch() also summarizes model stats as means for
different estimated values of k. It is designed to take files produced
by the BootKmeans() function as input, but other data can be analyzed
if the descriptions of the data formats given below are observed
carefully.
}
\details{
If you publish data or results produced with MHCtools, please cite both of
the following references:
Roved, J. (2022). MHCtools: Analysis of MHC data in non-model species. Cran.
Roved, J. (2024). MHCtools 1.5: Analysis of MHC sequencing data in R. In S.
Boegel (Ed.), HLA Typing: Methods and Protocols (2nd ed., pp. 275–295).
Humana Press. https://doi.org/10.1007/978-1-0716-3874-3_18
}
\examples{
filepath <- system.file("extdata/ClusterMatch", package="MHCtools")
path_out <- tempdir()
k_summary_table <- k_summary_table
ClusterMatch(filepath, path_out, k_summary_table)
}
\seealso{
\code{\link{BootKmeans}}
}
