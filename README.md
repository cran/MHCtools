## README

**Welcome to the R package MHCtools**

This package contains a set of ten tools for analysis of major histocompatibility complex (MHC) data in non-model species. The functions are tailored for amplicon data sets that have been filtered using the 'dada2' method (for more information visit <https://benjjneb.github.io/dada2>). Each of the functions are described below. For usage examples, please inspect the help pages for each function.
***
**Exporting FASTA files**

The CreateFas function() creates a fasta file with all the sequences in the data set.
The CreateSamplesFas() function creates a fasta file for each sample in the data set.
The input to both functions is the sequence table output by dada2.
***
**Replicate matching in amplicon filtering**

In amplicon filtering it is sometimes valuable to compare technical replicates in order to estimate the accuracy of a genotyping experiment. This may be done both to optimize filtering settings and to estimate repeatability before publication. The function ReplMatch() is designed to automatically compare technical replicates in an amplicon filtering data set.

The ReplMatch() function outputs a set of R lists containing for each replicate set the observed sequence variants, the names of the sequences that were incongruent in the replicates, and the mean proportion of incongruent sequences (if 100% matches are expected between the replicates, this is equivalent of an error rate in the sequencing process).

To evaluate the output, the GetReplTable() function will use the output files to produce a table with the replicate sets and their respective mean proportion of incongruent sequences. The GetReplStats() function will use the output files to calculate the number of replicate sets with zero incongruent sequences, the proportion of replicate sets with zero incongruent sequences, the mean of the mean proportion of incongruent sequences across all replicate sets, and the repeatability of the sequencing experiment.
***
**MHC Haplotype Inference**

The HpltFind() function is designed to automatically infer MHC haplotypes from the genotypes of parents and offspring in families in non-model species, where MHC sequence variants cannot be identified as belonging to individual loci. Knowing the haplotypes in such species can be a valuable source of information for several reasons, e.g.:

* The MHC alleles on each haplotype are often tightly linked, and therefore the combined effects of the alleles on a haplotype ultimately affect the fitness of individuals.
* Knowing the haplotypes allows for statistical analyses of effects of individual MHC alleles, as the presence of other linked alleles can be controlled for. (This underlines one of the major issues with the popular but problematic practice of grouping alleles into functional supertypes to gain statistical power.)
* Knowing the haplotypes can be very valuable in studies that attempt to sequence the entire genomic region of the MHC.

The HpltFind() function outputs a set of R lists containing for each family (or nest) the putative haplotypes, the names of sequences that could not be resolved with certainty in each parent, the names of the sequences that were incongruent in the genotypes of the family, and the mean proportion of incongruent sequences (which is a measure of the haplotype inference success and largely influenced by the exactness of the genotyping experiment).

To evaluate the output, the GetHpltTable() function will use the output files to produce a table with the mean proportion of incongruent sequences for each family. If the mean proportion of incongruent sequences is generally low, but certain nests have many incongruent sequences, biological reasons may be causing the mismatches, e.g. extra-pair fertilizations or recombination events. The GetHpltStats() function will use the output files to calculate the mean of the mean proportion of incongruent sequences across all nests in the data set.
***
**Parent pair diversity**

The PapaDiv() function is useful for heritability analyses in non-model species, where one wants to estimate the heritability of MHC diversity. The function calculates the joint diversity in parent pairs, taking into account alleles that are shared between the parents.

The PapaDiv() function outputs a set of R lists containing for the joint diversity of each parent pair, the proportion of sequences that are shared between the parents, the diversity of each of the parents, the observed sequence variants in each parent, the matched sequence variants, and the incongruent sequence variants in each parent.

For downstream data analysis, the PapaDiv() function produces a summary table with the names of the parents in a pair, their respective MHC diversities, and the joint parent pair diversity.
***
**P-distance calculations**

The CalcPdist() function is useful for calculating p-distances from pairwise sequence comparisons.

The CalcPdist() function takes a dada2 sequence table or a fasta file as input and produces a matrix with all pairwise p-distances for all sequences in the data set. If a dada2 sequence table is provided as input, the function furthermore produces a table with the mean p-distance from all pairwise comparisons of the sequences in each sample in the data set. 

It includes an option for the user to specify which codons to compare, which is useful e.g. if conducting the analysis only on codon positions involved in specific functions such as the peptide-binding of an MHC molecule. It also includes an option to calculate amino acid p-distances from protein-coding sequences using a standard genetic code.
***
*Copyright Jacob Roved*
