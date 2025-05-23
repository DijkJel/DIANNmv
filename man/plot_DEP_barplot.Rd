% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_DEP_barplot}
\alias{plot_DEP_barplot}
\title{Make bar plot showing number of differential proteins between conditions.}
\usage{
plot_DEP_barplot(res, comparisons = "all", names = NULL)
}
\arguments{
\item{res}{A data frame with results from get_DEPresults()}

\item{comparisons}{A character vector specifying which comparisons
to include. The default 'all' considers all comparisons present in res.
The bars are plotted in the order of specified comparisons.}

\item{names}{(optional). Specifies the labels for the different bars.
Defaults to the comparison names in colnames(res)}
}
\value{
A ggplot bar plot
}
\description{
Make bar plot showing number of differential proteins between conditions.
}
\examples{
se <- prepare_se(report.pg_matrix, expDesign)
res <- get_DEPresults(se, type  = 'all')

plot_DEP_barplot(res)
plot_DEP_barplot(res, comparisons = c('neg_ctrl_vs_motif1',
                                       'neg_ctrl_vs_motif2'))


}
