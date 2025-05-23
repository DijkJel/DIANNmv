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


#' Data frame with potential contaminants from MaxQuant.
#'
#'
#' @format ## `contaminants_maxquant`
#' A data frame with 246 rows and 1 columns.
#'
#' @source MaxQuant
'contaminants_maxquant'

#' List with calculated theoretical tryptic peptides for mouse and human for iBAQ
#' Peptides are calculated using 'cleaver' package with 'trypsin' set as enzyme,
#' zero miscleavages, and peptide length 7-30.
#'
#' @format ## `ibaq_peptides`
#' A list of length 2
#'
#' @source Uniprot
'ibaq_peptides'
