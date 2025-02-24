#' Adds a Potential.contaminant column to pg_matrix based on MaxQuant
#' contaminants.txt
#'
#' @param pg_matrix The report.pg_matrix output file from DIANN
#'
#' @return A data frame with a added Potential.contaminant column.
#' @export
#'
#' @examples
#' pg_matrix <- add_contaminants(report.pg_matrix)
add_contaminants = function(pg_matrix){

  pg_matrix$Potential.contaminant = sapply(pg_matrix$Protein.Group, function(x){

    genes = strsplit(x, ';')[[1]]
    is_contaminant = ifelse(any(genes %in% DIANNmv::contaminants_maxquant$uniprot_ID), '+', '')
    return(is_contaminant)
  })

  return(pg_matrix)
}

#' Prepare summarizedExperiment object from diann report.pg_matrix file
#'
#' @param pg_matrix the report.pg_matrix file from DIANN
#' @param expDesign A data frame with the experimental design. Should contain at least 'label', 'condition', and 'replicate' columns.
#' @param pr_matrix Optional argument. If the report.pr_matrix file from DIANN is provided, peptide information will be added to output.
#' @param missing_thr Integer specifying which proteinGroups are filtered out based on missing values.
#' @param min_peptides An integer specifing the cutoff for razor/unique peptides. The default is 0.
#' @param impute Specifies which imputatation method to use (default: knn). No imputation is done when entering 'none'.
#' Imputation is done using the 'MinProb' value from 'MsCoreUtils' package. See ?DEP::impute for all options.
#' @param remove_contaminants A logical value specifying if potential
#' contaminants should be removed from the pg_matrix.
#' @import DEP grid SummarizedExperiment
#'
#' @return A summarized Experiment object
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix,
#'                 expDesign, missing_thr = 1,
#'                 impute = 'knn') # creates se with missing values imputed
#'
#' se <- prepare_se(report.pg_matrix,
#'                 expDesign,
#'                  missing_thr = 1,
#'                  impute = 'none') # creates se without imputing missing values
prepare_se = function(pg_matrix, expDesign, pr_matrix = NULL, missing_thr = 0,
                      min_peptides = 0, impute = 'knn', remove_contaminants = TRUE){


  if (!is.null(pr_matrix)){
    pep = get_nPep_prMatrix(pr_matrix)
    pg_matrix = add_peptide_numbers(pg_matrix, pep)
  }


  pg_matrix = add_contaminants(pg_matrix)
  if(remove_contaminants){pg_matrix = pg_matrix[pg_matrix$Potential.contaminant != '+',]}
  pg_uniq = DEP::make_unique(pg_matrix, 'Genes', 'Protein.Group', delim = ';')

  pat = paste(expDesign$label, collapse = '|')
  lfq = grep(pat, colnames(pg_uniq))
  se = DEP::make_se(pg_uniq, lfq, expDesign)

  if('n_total' %in% colnames(pg_matrix)){
    pep = pg_matrix[,grep('npep', colnames(pg_matrix))]
    dimnames(pep) = list(rownames(se), colnames(se))
    SummarizedExperiment::rowData(se)$npep_total = pg_matrix$n_total
    SummarizedExperiment::assay(se, 'peptide_info') = pep
    se = se[SummarizedExperiment::rowData(se)$npep_total > min_peptides,]
  }


  se_filt = DEP::filter_missval(se, missing_thr)
  se = DEP::normalize_vsn(se_filt)
  p = DEP::plot_normalization(se_filt, se)
  grid::grid.newpage()
  grid::grid.draw(p)

  if(!impute == 'none'){se = DEP::impute(se, impute)}

  return(se)
}
