% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.functions.R
\name{genesetEnrichment}
\alias{genesetEnrichment}
\title{genesetEnrichment}
\usage{
genesetEnrichment(
  ranked_list,
  gene_set_list,
  max_size_gene_set = 100,
  min_size_gene_set = 10,
  fdr_correct = TRUE,
  fdr_threshold = 0.05,
  hypothesis = "two.sided"
)
}
\arguments{
\item{ranked_list}{A ranked list of genes based on a differential expression table (ranked by log2 FC or adj. p-values)}

\item{gene_set_list}{A list of gene sets to be evaluated for enrichment in the ranked gene list}

\item{max_size_gene_set}{The maximum size of geneset to be evaluated (default = 100)}

\item{min_size_gene_set}{The minimum size of geneset to be evaluated (default = 10)}

\item{fdr_correct}{Whether to perform Benjamini-Hochberg multiple hypothesis correction (FDR correction) (boolean - default = TRUE)}

\item{fdr_threshold}{The FDR threshold to use (ignored if fdr_correct = FALSE, default = 0.05)}

\item{hypothesis}{The hypothesis to be tested - one of "greater", "less", or "two.sided" (default = "two.sided")}
}
\value{
A dataframe with one row per geneset detailing the ranks of shared and not-shared genes and the
p-value of the enrichment analysis based on the ranked list provided.
}
\description{
Perform enrichment analysis, given a ranked list of genes, a list of genesets, and appropriate thresholds
}
\examples{
genesetEnrichment(prepRankedGeneList(DE_Seq2_output), GO_genesets, hypothesis = "greater")
}
