## NEWS

September 29 2017

Welcome to the R package MHCtools

This new R package contains nine useful functions for analysis of major histocompatibility complex (MHC) data in non-model species. The functions are tailored for amplicon data sets that have been filtered using the 'dada2' pipeline (for more information visit <https://benjjneb.github.io/dada2>), but even other data sets can be analyzed, if the data tables are formatted according to the description in each function.
The ReplMatch() function matches replicates in data sets in order to evaluate genotyping success.
The GetReplTable() and GetReplStats() functions perform such an evaluation.
The HpltFind() function infers putative haplotypes from families in the data set. 
The GetHpltTable() and GetHpltStats() functions evaluate the accuracy of the haplotype inference.
The PapaDiv() function compares parent pairs in the data set and calculate their joint MHC diversity, taking into account sequence variants that occur in both parents.
The CreateFas() function creates a fasta file with all the sequences in the data set.
The CreateSamplesFas() function creates a fasta file for each sample in the data set.

Version 1.0.0 is the initial release of the R package MHCtools. 


October 23rd 2017

MHCtools Version 1.1.0 

This update adds the new function MeanPdist() to the R package 'MHCtools'. The MeanPdist() function calculates the mean p-distance from pairwise comparisons of the sequences in each sample in a data set. The function includes an option for the user to specify which codons to compare, which is useful e.g. if conducting the analysis only on codon positions involved in specific functions such as the peptide-binding of an MHC molecule.


February 4th 2019

MHCtools Version 1.1.1 

This update fixes bugs in the HpltFind and PapaDiv functions. Furthermore, a minor modification of the functions GetHpltTable and GetReplTable ensures that results are presented according to nest or replicate number, respectively (previous versions presented them in alphanumeric order).


August 8th 2019

MHCtools Version 1.2.0

This update replaces the MeanPdist() function with the new function CalcPdist(), which has more universal applications. The CalcPdist function produces a matrix with p-distances from pairwise comparisons of all sequences in an input file. If a dada2 sequence table is used as input, the function furthermore produces a table with mean p-distances for each sample in the data set. Optionally, amino acid p-distances may be calculated, based on the standard genetic code. The function furthermore includes an option for the user to specify which codons to compare, which is useful e.g. if conducting the analysis only on codon positions involved in specific functions such as the peptide-binding of an MHC molecule. Input files may be either a fasta file (fasta files are accepted in the format rendered by e.g. the read.fasta() function from the package 'seqinr') or a dada2 sequence table.

August 8th 2019

MHCtools Version 1.2.1

This update provides a minor modification of the CalcPdist() function, so that when a dada2 sequence table is used as input, sequences are now named by an index number (corresponding to their column number in the input table) in the p-distance matrix.
