.onAttach <- function(libname, pkgname) {
  package_citation1 <- "Roved, J. (2022). MHCtools: Analysis of MHC data in non-model species. Cran. URL: https://cran.r-project.org/web/packages/MHCtools/index.html"
  package_citation2 <- "Roved, J. (2024). MHCtools 1.5: Analysis of MHC sequencing data in R. In S. Boegel (Ed.), HLA Typing: Methods and Protocols (2nd ed., pp. 275-295). Humana Press. https://doi.org/10.1007/978-1-0716-3874-3_18"
  packageStartupMessage("Thank you for using MHCtools!")
  packageStartupMessage("To acknowledge my work, please cite the following:")
  packageStartupMessage(paste(package_citation1))
  packageStartupMessage(paste(package_citation2))
}
