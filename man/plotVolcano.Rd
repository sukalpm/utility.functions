% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.functions.R
\name{plotVolcano}
\alias{plotVolcano}
\title{plotVolcano}
\usage{
plotVolcano(
  res_tbl,
  title,
  padj_thresh = 0.05,
  l2fc_thresh = 1.2,
  max_overlaps = 20
)
}
\arguments{
\item{res_tbl}{Output from differential expression analysis (e.g. from DESeq2)}

\item{title}{Title of plot}

\item{padj_thresh}{Adjusted p-value threshold for marking genes as significantly DE (default = 0.05)}

\item{l2fc_thresh}{Log2 FC threshold for marking genes as significantly DE (default = 1.2)}

\item{max_overlaps}{Passed to ggrepel to determine density of labeling genes (default = 20)}
}
\value{
A volcano plot of the given DE results table.
}
\description{
Plot a volcano plot from a given differential expression table (x-axis = Log2(fold change), and y-axis = -Log2(adjusted p-value))
}
\examples{
plotVolcano(DE_Seq2_output, "Volcano plot - disease vs. control")
}
