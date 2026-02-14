## NEWS

### MHCtools Version 1.6.0  

*February 14th, 2026*  

This update introduces the new SynDist() function to MHCtools. SynDist() performs qualitative and quantitative analyses of synonymous variation among aligned protein-coding DNA sequences, that is, nucleotide substitutions that do not translate to changes in the amino acid sequences due to degeneracy of the genetic code.

The BootKmeans() function was updated with a progress reporter, that prints loop iteration to the console while the function is running.  

### MHCtools Version 1.5.5  

*January 17th, 2025*  

This update fixed an error in the citation info in MHCtools.  

### MHCtools Version 1.5.4  

*January 15th, 2025*  

This update fixed a bug in the BootKmeans() function that caused an error when users tried to load the saved .RData cluster files using the load() function.  

In addition, a couple of notes were added to the description of BootKmeans() in the reference manual. 

The citation details were updated with the following:  
Roved, J. (2024). MHCtools 1.5: Analysis of MHC sequencing data in R. In S. Boegel (Ed.), HLA Typing: Methods and Protocols (2nd ed., pp. 275–295). Humana Press. https://doi.org/10.1007/978-1-0716-3874-3_18

### MHCtools Version 1.5.3  

*July 8th, 2023*  

This update fixes a bug in the DistCalc() function.

### MHCtools Version 1.5.2  

*June 30th, 2023*  

This update fixes a bug in the DistCalc() function and introduces a few minor edits in the reference manual.  

### MHCtools Version 1.5.1  

*March 22nd, 2023*  

This update fixes a bug in the ReplMatch() and PapaDiv() functions, which caused an error when matching the occurrence of sequences between samples. Sequences are now named using index numbers of the sequences with zeroes padded in front, e.g. "Sequence_001" to "Sequence_999". This prevents RegEx pattern matching between e.g. "Sequence_1" and "Sequence_1X". This naming concept has been implemented throughout MHCtools, for the sake of consistency. The updated functions are: CreateFas(), CreateSamplesFas(), DistCalc(), HpltFind(), ReplMatch(), and PapaDiv().

In addition, a bug was fixed in the HpltMatch() function, which caused an error when running the function with the setting threshold=NULL. 

### MHCtools Version 1.5.0  

*October 19th, 2022*  

**New functions**  
This update introduces three new functions to MHCtools, which are all related to haplotype inference:  
  
The CreateHpltOccTable() function creates a haplotype-sequence occurrence matrix from the output of HpltFind(), for easy overview of which sequences are present in which haplotypes.  
  
The HpltMatch() function compares haplotypes to help identify overlapping and potentially identical types.  
  
The NestTablesXL() function translates the output from HpltFind() to an Excel workbook, that provides a convenient overview for evaluation and curation of the inferred putative haplotypes.  

More details on these new functions are provided in the reference manual.  

**Adjustment to DistCalc()**  
The DistCalc() function has been updated so that it can handle data sets where some samples have 0 or 1 sequence(s) in a sequence occurrence table. DistCalc() will assign NA to the mean distance value for such samples in the output table.  

**Adjustment to ReplMatch()**  
The ReplMatch() function has been updated to throw a warning if all samples in a replicate set have 0 sequences. ReplMatch() will report which replicate set is problematic and suggest to remove it from the replicates table.     

**Miscellaneous**  
Finally, a few minor adjustments have been made to improve the efficiency of the code in most functions, with no change in functionality.  

### MHCtools Version 1.4.3  

*August 15th, 2022*  

This update fixes a bug in the HpltFind function, which caused an error when no alleles are shared between a chick and one of its parents (this may e.g. occur in nests with extra-pair fertilization).  

Furthermore, the updated version of HpltFind introduces the alpha-parameter, which may be used to fine-tune the haplotype analysis (details provided in the reference manual).  

### MHCtools Version 1.4.2  

*May 23rd, 2022*  

This version provides the following citation details following presentation of MHCtools in Molecular Ecology Resources:  
Roved, J., Hansson, B., Stervander, M., Hasselquist, D., & Westerdahl, H. (2022). MHCtools - an R package for MHC high-throughput sequencing data: genotyping, haplotype and supertype inference, and downstream genetic analyses in non-model organisms. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13645

### MHCtools Version 1.4.1  

*October 11th, 2021*  

This update fixes a bug in the DistCalc() function, which caused an error when attempting to calculate p-distances from a fasta file of amino acid sequences. 

### MHCtools Version 1.4.0  

*September 13th, 2021*  

**The BootKmeans() and ClusterMatch() functions**  
In this update, two new functions BootKmeans() and ClusterMatch() have been added. In MHC data analysis, it is often desirable to group alleles by their physico-chemical properties, as MHC receptors with similar properties share the repertoire of peptides they can bind. From a functional immunological perspective, alleles with similar properties may therefore be regarded as belonging to the same supertypes, which in many cases can simplify statistical analyses by reducing the number of independent variables and increase statistical power as more samples will share supertypes compared to alleles. Inference of MHC supertypes has traditionally been carried out by k-means clustering analysis on a set of z-descriptors of the physico-chemical properties of the amino acid sequences. However, the inference of relevant clusters is not always straightforward, since MHC data sets do not always produce clear inflection points (e.g. the elbow in an elbow plot). The BootKmeans() function is a wrapper for the kmeans() function of the stats package, which allows for bootstrapping of a k-means clustering analysis. BootKmeans() performs multiple runs of kmeans() and estimates optimal k-values based on a user-defined threshold of BIC reduction. The method may be seen as an automated and bootstrapped version of visually inspecting elbow plots of BIC- vs. k-values.
To evaluate which of the bootstrapped k-means models is most accurate and/or informative, the ClusterMatch() function offers a tool for evaluating whether different k-means clustering models identify similar clusters, and summarize bootstrap model stats as means for different estimated values of k. ClusterMatch() is designed to take files produced by the BootKmeans() function as input, but other data can be analyzed if the descriptions of the required data formats are observed carefully. 

### MHCtools Version 1.3.0  

*September 15th, 2020*  

**The DistCalc() function**  
In this update, the new function DistCalc() replaces CalcPdist(). The DistCalc() function calculates Grantham, Sandberg, or p-distances and produces a matrix with distances from pairwise comparisons of all sequences in an input file. Input files may be either a fasta file (in the format rendered by the read.fasta() function from the package seqinr) or a dada2-style sequence table. If a dada2-style sequence table is used as input, the function produces a table with mean distances for each sample in the data set in addition to the distance matrix. Amino acid distances may be calculated from input files with nucleotide sequences by automatic translation using the standard genetic code. The function furthermore includes an option for the user to specify which codons to compare, which is useful e.g. if conducting the analysis only on codon positions involved in specific functions such as the peptide-binding of an MHC molecule.

**Z-descriptors**  
If calculation of Sandberg distances is specified when using DistCalc(), the function additionally outputs five tables with physico-chemical z-descriptor values for each amino acid position in all sequences in the input file. These tables may be useful for further downstream analyses, such as estimation of MHC supertypes. The z-descriptor values are derived from Sandberg et al. 1998, JMed Chem. 41(14):2481-2491.

**Miscellaneous**  
This update furthermore removes the dependency on the package rlist.

### MHCtools Version 1.2.1  

*August 8th, 2019*  

This update provides a minor modification of the CalcPdist() function, so that when a dada2 sequence table is used as input, sequences are now named by an index number (corresponding to their column number in the input table) in the p-distance matrix.

### MHCtools Version 1.2.0  

*August 8th, 2019*  

This update replaces the MeanPdist() function with the new function CalcPdist(), which has more universal applications. The CalcPdist function produces a matrix with p-distances from pairwise comparisons of all sequences in an input file. If a dada2 sequence table is used as input, the function furthermore produces a table with mean p-distances for each sample in the data set. Optionally, amino acid p-distances may be calculated, based on the standard genetic code. The function furthermore includes an option for the user to specify which codons to compare, which is useful e.g. if conducting the analysis only on codon positions involved in specific functions such as the peptide-binding of an MHC molecule. Input files may be either a fasta file (fasta files are accepted in the format rendered by e.g. the read.fasta() function from the package seqinr) or a dada2 sequence table.

### MHCtools Version 1.1.1  

*February 4th, 2019* 

This update fixes bugs in the HpltFind and PapaDiv functions. Furthermore, a minor modification of the functions GetHpltTable and GetReplTable ensures that results are presented according to nest or replicate number, respectively (previous versions presented them in alphanumeric order).

### MHCtools Version 1.1.0  

*October 23rd, 2017*  

This update adds the new function MeanPdist() to the R package MHCtools. The MeanPdist() function calculates the mean p-distance from pairwise comparisons of the sequences in each sample in a data set. The function includes an option for the user to specify which codons to compare, which is useful e.g. if conducting the analysis only on codon positions involved in specific functions such as the peptide-binding of an MHC molecule.

### Welcome to MHCtools Version 1.0.0

*September 29th, 2017*  

This new R package contains nine useful functions for analysis of major histocompatibility complex (MHC) data in non-model species. The functions are tailored for amplicon data sets that have been filtered using the dada2 pipeline (for more information visit <https://benjjneb.github.io/dada2/>), but even other data sets can be analyzed, if the data tables are formatted according to the description in each function.  
The ReplMatch() function matches replicates in data sets in order to evaluate genotyping success.  
The GetReplTable() and GetReplStats() functions perform such an evaluation.  
The HpltFind() function infers putative haplotypes from families in the data set.  
The GetHpltTable() and GetHpltStats() functions evaluate the accuracy of the haplotype inference.  
The PapaDiv() function compares parent pairs in the data set and calculate their joint MHC diversity, taking into account sequence variants that occur in both parents.  
The CreateFas() function creates a fasta file with all the sequences in the data set.  
The CreateSamplesFas() function creates a fasta file for each sample in the data set.  
