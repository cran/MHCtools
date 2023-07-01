## README  

### Welcome to the R package MHCtools  

MHCtools contains fifteen tools for bioinformatics processing and analysis of major histocompatibility complex (MHC) data. The functions are tailored for amplicon sequencing data sets that have been filtered using the dada2 method (Callahan et al. 2016; for more information visit <https://benjjneb.github.io/dada2/>), but even other data sets can be analyzed. Each of the functions are described below. For usage examples, please inspect the help pages for each function.  
  
### Evolutionary and functional differences between sequences  

The DistCalc() function is a useful tool for calculating distances from pairwise sequence comparisons. It offers calculation of Grantham distances (Grantham 1974), Sandberg distances (Sandberg et al. 1998), or simple p-distances (i.e., proportion of variable codons) in pairwise comparisons of aligned sequences. When calculating Sandberg distances, the function additionally outputs five tables with physico-chemical z-descriptor values (Sandberg et al. 1998) for each amino acid position in all sequences in the data set. These tables may be useful for further downstream analyses, such as estimation of MHC supertypes.  

The DistCalc() function takes a fasta file or a dada2-style sequence occurrence table (with aligned sequences as column names and samples in rows) as input and produces a matrix with pairwise distances for all sequences in the data set. If a dada2-style sequence occurrence table is provided as input, the DistCalc() function furthermore produces a table with the mean distances from all pairwise comparisons of the sequences in each sample in the data set.  

The DistCalc() function includes an option for the user to specify which codons to compare, which is useful e.g. if conducting the analysis only on codons that are involved in specific functions, such as the peptide-binding of an MHC molecule. It also includes an option to calculate amino acid distances directly from protein-coding DNA sequences using the standard genetic code.  
  
### MHC supertype inference 

In MHC data analysis, it is often desirable to group alleles by their physico-chemical properties, as MHC receptors with similar properties share the repertoire of peptides they can bind. From a functional immunological perspective, alleles with similar properties may therefore be regarded as belonging to the same supertypes, which in many cases can be exploited to simplify statistical analyses by reducing the number of independent variables and increase statistical power (i.e., as relatively more samples will share supertypes compared to alleles). Inference of MHC supertypes has traditionally been carried out by k-means clustering analysis on a set of z-descriptors of the physico-chemical properties of the amino acid sequences (Sandberg et al. 1998), and the DistCalc() function described above offers to produce such descriptors from amino acid sequences. However, the inference of relevant clusters of MHC alleles is not always straightforward, since MHC data sets might not always produce clear inflection points (e.g. the elbow in an elbow plot). 

As a solution to this problem, the BootKmeans() function offers bootstrapping of k-means clustering analysis for greatly improved confidence in the estimated clusters - i.e., the MHC supertypes. BootKmeans() is a wrapper for the kmeans() function of the stats package and performs multiple runs of kmeans() while estimating optimal k-values based on a set threshold for the step-wise reduction in BIC. The method may be seen as an automated and bootstrapped version of visually inspecting elbow plots of BIC- vs. k-values.

To evaluate which in a set of bootstrapped k-means models is most accurate and/or informative, the ClusterMatch() function offers a tool for evaluating whether different k-means clustering models identify similar clusters, and summarize bootstrap model stats as means for different estimated values of k. ClusterMatch() is designed to take files produced by the BootKmeans() function as input, but other data can be analysed if the descriptions of the required data formats are observed carefully. 
  
### MHC haplotype inference  

The HpltFind() function is designed to automatically infer MHC haplotypes from the genotypes of parents and offspring in families in non-model species, where MHC sequence variants cannot be identified as belonging to individual loci. Knowing the haplotypes in such species can be a valuable source of information for several reasons, e.g.:  
  
* The MHC alleles on each haplotype are often tightly linked, and therefore the combined effects of the alleles on a haplotype ultimately affect the fitness of individuals.  
* Knowing the haplotypes allows for statistical analyses of effects of individual MHC alleles, as the presence of other linked alleles can be controlled for.  
* Knowing the haplotypes can be valuable in studies that attempt to sequence the entire genomic region of the MHC.  
  
The HpltFind() function outputs a set of R lists containing for each family the putative haplotypes, the names of sequences that could not be resolved with certainty in each parent, the names of the sequences that were incongruent in the genotypes of the family, and the mean proportion of incongruent sequences (which is a measure of the haplotype inference success and largely influenced by the exactness of the genotyping experiment).  
  
To evaluate the output, the GetHpltTable() function will use the output files to produce a table with the mean proportion of incongruent sequences for each family. If the mean proportion of incongruent sequences is generally low, but certain families have many incongruent sequences, biological reasons may be causing the mismatches, e.g. extra-pair fertilizations or recombination events. The GetHpltStats() function will use the output files to calculate the mean of the mean proportion of incongruent sequences across all families in the data set.  
  
The CreateHpltOccTable() function reads the R lists output by the HpltFind() function and creates a binary (logical) haplotype-sequence occurrence matrix, which provides an easy overview of which sequences that are present in which haplotypes.  
  
As the HpltFind() function infers haplotypes based on the segregation patterns within families, some inferred putative haplotypes may be identical to each other (i.e., reoccur in different families), or they may be in a gray zone if they differ only by incongruent or unresolved sequences. It is therefore useful to compare putative haplotypes to help identify overlapping (i.e., potentially identical) ones as candidates for further investigation. However, the number of pairwise comparisons easily gets very large. The HpltMatch() function calculates the proportion of matching sequences between pairs of haplotypes and outputs a table with the resulting values. In addition, if a threshold is specified, a list of haplotype matches where the proportion of matching sequences exceeds the threshold will be produced. This is highly convenient for reducing the number of haplotype pairs for further investigation to a manageable size.  
  
The NestTablesXL() function reads the R lists output by the HpltFind() function and translates them to an Excel workbook, providing a convenient overview for evaluation and curating of the inferred putative haplotypes. The workbook contains separate tabs for each family (nest) in the data set and provides an overview of the genotypes of the samples in each family (in form of a sequence occurrence matrix) and the inferred haplotypes, including lists of unresolved and incongruent sequences. The tables are highly convenient for investigating of segregation patterns of individual sequences, e.g. when comparing potentially identical haplotypes or tracking sequence segregation patterns between families or across generations.  
  
### Parent pair diversity  

In studies of mate choice, fitness, and heritability, it can be of interest to analyze similarities and differences between the MHC genotypes of parents in families. The PapaDiv() function calculates the joint diversity in parent pairs, taking into account alleles that are shared between the parents.  

The PapaDiv() function outputs a set of R lists containing for the joint diversity of each parent pair, the proportion of sequences that are shared between the parents, the diversity of each of the parents, the observed sequence variants in each parent, the matched sequence variants, and the incongruent sequence variants in each parent.  

For downstream data analysis, the PapaDiv() function produces a summary table with the names of the parents in a pair, their respective MHC diversities, and the joint parent pair diversity.  
  
### Replicate matching  

In amplicon filtering it is often necessary to compare technical replicates in order to estimate the accuracy of a genotyping experiment. This may be done both to optimize filtering settings and to estimate repeatability. The function ReplMatch() is designed to automatically compare technical replicates in an amplicon filtering data set.  

The ReplMatch() function outputs a set of R lists containing for each replicate set the observed sequence variants, the names of the sequences that were incongruent in the replicates, and the mean proportion of incongruent sequences (if 100% matches are expected between the replicates, this is equivalent to an error rate in the genotyping process).  

To evaluate the output, the GetReplTable() function will use the output files to produce a table with the replicate sets and their respective mean proportion of incongruent sequences. The GetReplStats() function will use the output files to calculate the number of replicate sets with zero incongruent sequences, the proportion of replicate sets with zero incongruent sequences, the mean of the mean proportion of incongruent sequences across all replicate sets, and the repeatability of the sequencing experiment.  
  
### Exporting FASTA files  

CreateFas() and CreateSamplesFas() are two simple, but useful tools. The CreateFas() function creates a fasta file with all the sequences from a dada2-style sequence occurrence table. The CreateSamplesFas() function similarly creates a unique fasta file for each sample in the table.  
  
### Citation
  
If you publish data or results produced with MHCtools, please cite both of the following references: 
Roved, J. 2022. MHCtools: Analysis of MHC data in non-model species. Cran. https://cran.r-project.org/web/packages/MHCtools/index.html
Roved, J., Hansson, B., Stervander, M., Hasselquist, D., & Westerdahl, H. 2022. MHCtools - an R package for MHC high-throughput sequencing data: genotyping, haplotype and supertype inference, and downstream genetic analyses in non-model organisms. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13645  

### References  

Callahan, B.J., Mcmurdie, P.J., Rosen, M.J., Han, A.W., Johnson, A.J.A., Holmes, S.P. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nat. Methods 13.  
Grantham R. 1974. Amino acid difference formula to help explain protein evolution. Science 185:862-864.  
Sandberg M, Eriksson L, Jonsson J, Sjostrom M, Wold S. 1998. New chemical descriptors relevant for the design of biologically active peptides. A multivariate characterization of 87 amino acids. JMed Chem. 41(14):2481-2491.  
  
*Copyright Jacob Roved*  
