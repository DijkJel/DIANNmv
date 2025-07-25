% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsea_plots.R
\name{plot_gsea_bubbleplot}
\alias{plot_gsea_bubbleplot}
\title{Plots GSEA data as bubble plot. Mainly useful to compare multiple gsea comparisons.}
\usage{
plot_gsea_bubbleplot(
  gsea,
  ...,
  sample_names = NULL,
  padj_cutoff = 0.05,
  top_n = Inf,
  remove_prefix = F,
  max_name_length = Inf,
  break_names = T
)
}
\arguments{
\item{gsea}{Output from \link{perform_GSEA}.}

\item{...}{Additional gsea outputs to be included in visualization.}

\item{sample_names}{Character vector specifying the comparisons that were made in the GSEA. Optional.}

\item{padj_cutoff}{padj_cutoff The maximum p.adjust value allowed for inclusion of the pathway.}

\item{top_n}{he maximum number of pathways included. Takes the top_n pathways
with the lowest p.adj values.}

\item{remove_prefix}{Boolean specifying to remove the prefix from pathway names.}

\item{max_name_length}{Numeric value specifying the max length of pathway names.}

\item{break_names}{Boolean value that indicates if pathway names should be shown
on two lines when the length exceedds the max_name_length parameter}
}
\value{
A ggplot object
}
\description{
Plots GSEA data as bubble plot. Mainly useful to compare multiple gsea comparisons.
}
\examples{
\dontrun{

gsea1 <- perform_GSEA(res1, genesets)
gsea2 <- perform_GSEA(res2, genesets)

#Show data of a single gsea as one column
bubbleplot <- plot_gsea_bubbleplot(gsea1, top_n = 10, sample_names = 'Comparison_1')

# Show data of multiple GSEAs as multiple columns.
bubbleplot <- plot_gsea_bubbleplot(gsea1, gsea2, top_n = 10, sample_names = c('test1', 'test2'))
}
}
