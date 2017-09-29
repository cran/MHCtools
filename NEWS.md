## NEWS

September 29 2017

Welcome to the R package MHCtools

This new R package contains nine useful functions for analysis of major histocompatibility complex (MHC) data in non-model species. The functions are tailored for amplicon data sets that have been filtered using the 'dada2' pipeline (for more information visit <https://benjjneb.github.io/dada2>), but even other data sets can be analysed, if the data tables are formatted according to the description in each function.
The ReplMatch() function matches replicates in data sets in order to evaluate genotyping success.
The GetReplTable() and GetReplStats() functions perform such an evaluation.
The HpltFind() function infers putative haplotypes from families in the data set. 
The GetHpltTable() and GetHpltStats() functions evaluate the accuracy of the haplotype inference.
The PapaDiv() function compares parent pairs in the data set and calculate their joint MHC diversity, taking into account sequence variants that occur in both parents.
The CreateFas() function creates a fasta file with all the sequences in the data set.
The CreateSamplesFas() function creates a fasta file for each sample in the data set.

Version 1.0.0 is the initial release of the R package MHCtools. 
