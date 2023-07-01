.onAttach <- function(libname, pkgname) {
  package_citation1 <- "Roved, J. 2022. MHCtools: Analysis of MHC data in non-model species. Cran. URL: https://cran.r-project.org/web/packages/MHCtools/index.html"
  package_citation2 <- "Roved, J., Hansson, B., Stervander, M., Hasselquist, D., & Westerdahl, H. 2022. MHCtools - an R package for MHC high-throughput sequencing data: genotyping, haplotype and supertype inference, and downstream genetic analyses in non-model organisms. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13645"
  packageStartupMessage("Thank you for using MHCtools!")
  packageStartupMessage("To acknowledge my work, please cite the following:")
  packageStartupMessage(paste(package_citation1))
  packageStartupMessage(paste(package_citation2))
}
