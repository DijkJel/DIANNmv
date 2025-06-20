% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsea_plots.R
\name{plot_gsea_volcano}
\alias{plot_gsea_volcano}
\title{Plots GSEA data as volcano plot}
\usage{
plot_gsea_volcano(
  gsea,
  padj_cutoff = 0.05,
  label = "sig",
  top_n = NULL,
  up_color = "red3",
  down_color = "dodgerblue",
  ns_color = "grey70"
)
}
\arguments{
\item{gsea}{gsea Output from \link{perform_GSEA}.}

\item{padj_cutoff}{padj_cutoff The maximum p.adjust value allowed for inclusion of the pathway.}

\item{label}{Specifies which points to label. The default is 'sig', labeling
all significant points. Entering a value for top_n limits the labeling to
the top_n up- and top_n down-regulated proteins based on the p.adj.
When providing a vector with protein names, only those points are labeled.}

\item{top_n}{Specifies how many significant points to label.}

\item{up_color}{Character string specifying the color for significant points with positive NES.}

\item{down_color}{Character string specifying the color for significant points with positive NES.}

\item{ns_color}{Character string specifying the color for non-significant points.}
}
\value{
A ggplot2 object
}
\description{
Plots GSEA data as volcano plot
}
\examples{
\dontrun{
gsea <- perform_GSEA(res, genesets)
volcano <- plot_gsea_volcano(gsea, top_n = 10)
}
}
