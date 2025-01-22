#' Experimental design for example data
#'
#' Necessary for creating SummarizedExperiment object with MS data
#'
#' @format ## `expDesign`
#' A data frame with 6 rows and 3 columns:
#' \describe{
#'   \item{label}{column name in report.pg_matrix specifying MS sample}
#'   \item{condition}{The different experimental conditions}
#'   \item{replicate}{The different replicates per condition}
#' }
#' @source 241220_V118_LK_DNA_PD
'expDesign'

#' Output file of DIANN with raw proteinGroup intensities on protein level.
#'
#' Contains 9 samples with raw intensities after analysis with DIANN.
#'
#' @format ## `report.pg_matrix`
#' A data frame with 5948 rows and 13 columns.
#'
#' @source 241220_V118_LK_DNA_PD
'report.pg_matrix'

#' Output file of DIANN with raw intensities on peptide level.
#'
#' Contains 9 samples with raw peptide intensities after analysis with DIANN.
#'
#' @format ## `report.pr_matrix`
#' A data frame with 97329 rows and 19 columns.
#'
#' @source 241220_V118_LK_DNA_PD
'report.pr_matrix'
