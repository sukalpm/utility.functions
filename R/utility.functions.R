#Package utility.functions authored by Sukalp Muzumdar, CSHL

#Utility function to read in 10X single-cell data (output folder from CellRanger)
#'readCellRangerFolder
#'Prepare a ranked list of genes from a given DESeq2 output data.frame
#'
#' @param data.dir Folder containing matrix.tsv, genes/features.tsv, and barcodes.tsv (output from CellRanger pipeline)
#' @param proj.name Meaningful sample name for Seurat object
#' @param min.cells Subset feature matrix to only features detected in this many cells
#' @param min.features Only include cells with these many detected features
#' @return A Seurat object which can be used for further downstream analyses
#' @export
#'
#' @examples
#' readCellRangerFolder(path_to_cell_ranger_output, "Single-cell dataset", min.cells = 2, min.features = 100)

readCellRangerFolder <- function(data.dir, proj.name, min.cells = 3, min.features = 200){
  if (dir.exists(data.dir) == FALSE){
    return -1
  }
  else{
    count.data <- Read10X(data.dir = data.dir, gene.column = 1, strip.suffix = TRUE)

    data <- CreateSeuratObject(counts = count.data, project = proj.name,
                               min.cells = min.cells, min.features = min.features)
    return (data)
  }
}

#Prep gene list for M-W-style DE
#' prepRankedGeneList
#'Prepare a ranked list of genes from a given DESeq2 output data.frame
#'
#' @param filt_de_tbl Output from DESeq2 with all genes with padj = NA excluded
#'
#' @return A list of ranks with gene names as the labels.
#' @export
#'
#' @examples
#' prepRankedGeneList(DE_Seq2_output)
prepRankedGeneList <- function(de_tbl = filt_de_tbl){
  l2fc <- filt_de_tbl %>% dplyr::select(log2FoldChange) %>% "[["(1)
  names(l2fc) <- filt_de_tbl %>% dplyr::select(gene) %>% "[["(1)
  rank_l2fc <- rank(l2fc, ties.method = "average")
  return(rank_l2fc)
}


#Mann-Whitney-style enrichment on whole DE table
#' genesetEnrichment
#'
#'Perform enrichment analysis, given a ranked list of genes, a list of genesets, and appropriate thresholds
#' @param ranked_list A ranked list of genes based on a differential expression table (ranked by log2 FC or adj. p-values)
#' @param gene_set_list A list of gene sets to be evaluated for enrichment in the ranked gene list
#' @param max_size_gene_set The maximum size of geneset to be evaluated (default = 100)
#' @param min_size_gene_set The minimum size of geneset to be evaluated (default = 10)
#' @param fdr_correct Whether to perform Benjamini-Hochberg multiple hypothesis correction (FDR correction) (boolean - default = TRUE)
#' @param fdr_threshold The FDR threshold to use (ignored if fdr_correct = FALSE, default = 0.05)
#' @param hypothesis The hypothesis to be tested - one of "greater", "less", or "two.sided" (default = "two.sided")
#'
#' @return
#' A dataframe with one row per geneset detailing the ranks of shared and not-shared genes and the
#' p-value of the enrichment analysis based on the ranked list provided.
#' @export
#'
#' @examples
#' genesetEnrichment(prepRankedGeneList(DE_Seq2_output), GO_genesets, hypothesis = "greater")
genesetEnrichment <-
  function(ranked_list, gene_set_list, max_size_gene_set = 100, min_size_gene_set = 10, fdr_correct = TRUE,
           fdr_threshold = 0.05, hypothesis = "two.sided"){


    gene_list <- names(ranked_list)

    subset_gene_set_list <- list()

    x <- lapply(names(gene_set_list), function(i){
      subset_gene_set_list[[i]] <<- gene_set_list[[i]][gene_set_list[[i]] %in% gene_list]
    })

    keep <- lapply(seq_along(subset_gene_set_list),
                   function(i){
                     ifelse(
                       ((length(subset_gene_set_list[[i]]) > min_size_gene_set) &
                          (length(subset_gene_set_list[[i]]) < max_size_gene_set)), TRUE, FALSE)

                   })

    filt_gene_set <- subset_gene_set_list[unlist(keep)]

    raw <- data.frame()

    x <- lapply(seq_along(filt_gene_set), function(i){

      raw <<- rbind(raw, data.frame(geneset = names(filt_gene_set)[i],
                                    n_shared = length(ranked_list[unlist(filt_gene_set[i])]),
                                    n_not_shared = length(ranked_list[-c(which(names(ranked_list) %in% unlist(filt_gene_set[i])))]),
                                    avg_rank_shared = mean(ranked_list[unlist(filt_gene_set[i])]),
                                    avg_rank_not_shared = mean(ranked_list[-c(which(names(ranked_list) %in% unlist(filt_gene_set[i])))]),
                                    pval = wilcox.test(ranked_list[unlist(filt_gene_set[i])],
                                                       ranked_list[-c(which(names(ranked_list) %in% unlist(filt_gene_set[i])))],
                                                       alternative = hypothesis)$p.value))
    })

    raw <- raw %>% mutate(padj = p.adjust(pval, "BH"))
    return(raw)
  }

#Fisher-style enrichment for a subset of genes (signature etc.)
#' fisher_enrichment
#'
#'Perform enrichment analysis using the hypergeometric (Fisher) test, given a list of selected genes, the other genes in the "universe",
#'a list of genesets, and appropriate thresholds
#' @param genesets The genesets to evaluate for enrichment
#' @param de_genes The selected genes
#' @param not_de_genes The remaining (unselected) genes in the universe
#' @param min_size_gene_set The minimum size of the geneset for enrichment analysis (default = 10)
#' @param max_size_gene_set The maximum size of the geneset for enrichment analysis (default = 101)
#'
#' @return
#' A dataframe with one row per geneset detailing the p-value of the enrichment analysis based on the genesets and gene lists provided.
#' @export
#'
#' @examples
#' fisher_enrichment(GO_genesets, DE_genes_list, not_DE_genes_list)
fisher_enrichment <- function(genesets, de_genes, not_de_genes, min_size_gene_set = 10, max_size_gene_set = 101){

  total <- union(de_genes, not_de_genes)

  subset_genesets <-
    genesets[unlist(lapply(seq_along(genesets),
                           function(i){
                             ifelse((length(genesets[[i]]) > min_size_gene_set) &
                                      (length(genesets[[i]]) < max_size_gene_set), TRUE,FALSE)
                           }))]

  x <- lapply(names(subset_genesets), function(i){
    subset_genesets[[i]] <<- genesets[[i]][genesets[[i]] %in% total]
  })


  res_enrich <- data.frame()
  x <- lapply(X = seq_along(subset_genesets), FUN = function(i){
    res_enrich <<- rbind(res_enrich,
                         data.frame(geneset = names(subset_genesets)[i],
                                    len_gene_set = length(subset_genesets[[i]]),
                                    n_shared = sum(de_genes %in% subset_genesets[[i]]),
                                    pval = 1 - phyper(sum(de_genes %in% subset_genesets[[i]]) - 1,
                                                      length(de_genes),
                                                      length(not_de_genes),
                                                      length(subset_genesets[[i]]))))
  })
  return(res_enrich)
}

#Quick rank matrix
#' rank_matrix
#'
#'Function to rank a matrix quickly while replacing any NA values
#'
#' @param M The matrix to be ranked
#' @param na_value Value to be substituted if any NAs are encountered (default = 0.5)
#'
#' @return The matrix \code{M} with values replaced by normalized ranks ranging from 0 to 1.
#' @export
#'
#' @examples
#' rank_matrix(matrix = M, na_value = 0.5)
rank_matrix = function(M, na_value = 0.5) {
  is_na = is.na(M)
  result = matrixStats::colRanks(M, dim = c(length(M),1), ties.method="average")
  result = result / (length(result) - sum(is_na))
  dim(result) = dim(M)
  result[is_na] = na_value
  dimnames(result) = dimnames(M)
  return(result)
}

#' plotVolcano
#'
#'Plot a volcano plot from a given differential expression table (x-axis = Log2(fold change), and y-axis = -Log2(adjusted p-value))
#' @param res_tbl Output from differential expression analysis (e.g. from DESeq2)
#' @param title Title of plot
#' @param padj_thresh Adjusted p-value threshold for marking genes as significantly DE (default = 0.05)
#' @param l2fc_thresh Log2 FC threshold for marking genes as significantly DE (default = 1.2)
#' @param max_overlaps Passed to ggrepel to determine density of labeling genes (default = 20)
#'
#' @return
#' A volcano plot of the given DE results table.
#' @export
#'
#' @examples
#' plotVolcano(DE_Seq2_output, "Volcano plot - disease vs. control")
plotVolcano <- function(res_tbl, title, padj_thresh = 0.05, l2fc_thresh = 1.2, max_overlaps = 20){
  library(ggExtra)
  library(ggrepel)
  p <- res_tbl %>% filter(!is.na(padj)) %>% ggplot() +
    aes(x = log2FoldChange, y = -log2(padj)) +
    geom_point(aes(color = abs(log2FoldChange) > l2fc_thresh & padj < padj_thresh), size = 2.5) +
    geom_vline(xintercept = -1 * l2fc_thresh, color = "blue", linetype = "dotted") +
    geom_vline(xintercept = l2fc_thresh, color = "blue", linetype = "dotted") +
    geom_hline(yintercept = -log2(padj_thresh), color = "red", linetype = "dotted") +
    ggtitle(title) +
    geom_text_repel(max.overlaps = max_overlaps, aes(x = log2FoldChange, y = -log2(padj), label = gene),
                    data = res_tbl %>%
                      filter(abs(log2FoldChange) > l2fc_thresh & -log2(padj) > -log2(l2fc_thresh))) +
    scale_color_manual(values = setNames(c("orangered4", "steelblue"), c(TRUE, FALSE))) +
    xlab("Log2 (FC)") +
    ylab("-Log (adj. P)") +
    theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = 0.1),
          legend.position = "none",
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 15))

  p <- ggMarginal(p, type = "densigram")
  return(p)
}

#' excludeSexChrGenes
#'
#'Remove sex-chromosome genes from a table with gene names as a column (e.g. from a DE table)
#' @param de_tbl Output from differential expression analysis (e.g. from DESeq2)
#' @param gene_chr_mapping Mapping of HGNC symbols to chromosome names/numbers (e.g. as obtained from biomaRt)
#'
#' @return
#' Table with genes on the X- and Y-chromosome removed.
#' @export
#'
#' @examples
#' excludeSexChrGenes(DE_Seq2_output, biomaRt_mapping)
excludeSexChrGenes <- function(de_tbl, gene_chr_mapping){
  m <- merge(de_tbl, gene_chr_mapping %>% dplyr::select(hgnc_symbol, chromosome_name) %>% dplyr::rename(gene = hgnc_symbol), by = "gene")
  return(m %>% filter(chromosome_name != "X" & chromosome_name != "Y"))
}
