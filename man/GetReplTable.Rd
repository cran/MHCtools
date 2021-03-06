% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetReplTable_func_20200914.R
\name{GetReplTable}
\alias{GetReplTable}
\title{GetReplTable function}
\usage{
GetReplTable(filepath)
}
\arguments{
\item{filepath}{is a user defined path to the folder where the output files
from the ReplMatch() function have been saved.}
}
\value{
A table with the mean proportion of incongruent sequences for each
  replicate set.
}
\description{
\code{\link{GetReplTable}} uses the output files produced by the ReplMatch()
function to produce a table with the replicate sets and their respective mean
proportion of incongruent sequences.
}
\details{
If you publish data produced with MHCtools, please cite:
Roved, J. 2020. MHCtools: Analysis of MHC data in non-model species. Cran.
Roved, J., Hansson, B., Stervander, M., Hasselquist, D., Westerdahl, H. 2020.
Non-random association of MHC-I alleles in favor of high diversity haplotypes
in wild songbirds revealed by computer-assisted MHC haplotype inference using
the R package MHCtools. bioRxiv.
}
\examples{
filepath <- system.file("extdata/ReplMatchOut/", package="MHCtools")
GetReplTable(filepath)
}
\seealso{
\code{\link{ReplMatch}}; \code{\link{GetReplStats}}
}
