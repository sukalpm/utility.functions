# utility.functions

R package containing general-purpose utility functions
Function list - 

1. readCellRangerFolder - Import CellRanger-processed 10X single-cell RNA sequencing data into R
2. prepRankedGeneList - Given a differential expression table, prepare a ranked list of genes by fold change
3. genesetEnrichment - Given a differential expression table and a list of genesets, perform geneset enrichment analysis using a Mann-Whitney test
4. fisher_enrichment - Given a list of differentially expressed genes, and a list of genesets, perform hypergeometric enrichment analysis on the genesets
5. rank_matrix - Given a matrix, efficiently rank all the entries
6. plotVolcano - Given a differential expression table, generate a volcano plot depicting the results
7. excludeSexChrGenes - Given a list of genes or a differential expression table, exclude all genes associated with the human X or Y chromosome


Authored by Sukalp Muzumdar
