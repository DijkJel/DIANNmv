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

#' Adds iBAQ values based on LFQ values to SummarizedExperiment
#'
#' @param se SummarizedExperiment object in prepare_se
#'
#' @import SummarizedExperiment
#'
#' @return A Summarized Experiment object with extra assay.
#' @export
#'
#' @examples
#' \dontrun{
#' se <- add_maxLFQ_iBAQ(se)
#' }

add_maxLFQ_iBAQ = function(se){

  ibaq_peps = SummarizedExperiment::rowData(se)[,'ibaq_peptides']
  intensities = 2^as.matrix(assay(se))
  max_lfq_ibaq = intensities / ibaq_peps

  SummarizedExperiment::assay(se, 'maxLFQ_iBAQ') = max_lfq_ibaq
  return(se)
}


#' Tidies sample names and parses an experimental design.
#'
#' @param pg_matrix The report.pg_matrix file
#' @param pr_matrix The report.pr_matrix file
#'
#' @details
#' This functions tidies sample names and prepares and experimental design,
#' assuming that the structure of the sample names is <prefix_condition_replicate>.
#' An example of the prepared sample names and expDesign is printed in the console.
#'
#'
#' @return A list with a tidied pg_matrix, pr_matrix, and parsed expDesign.
#' @export
#'
#' @examples
#' tidy_data <- prepare_diann_data(report.pg_matrix, report.pr_matrix)
#' pg_matrix <- tidy_data$pg_matrix
#' pr_matrix <- tidy_data$pr_matrix
prepare_diann_data = function(pg_matrix, pr_matrix){

  cn = colnames(pg_matrix)[5:ncol(pg_matrix)]
  cn = sapply(cn, function(x){strsplit(x, '_|\\.')[[1]]})

  max_len = max(lengths(cn))
  vals = sapply(1:max_len, function(x){vals = sapply(cn, function(y){y[x]})})
  vals = vals[,!apply(vals, 2, function(x){length(unique(x)) == 1}), drop = F]

  cn = apply(vals, 1, function(x){paste(na.omit(x), collapse = '_')})

  colnames(pg_matrix)[5:ncol(pg_matrix)] = cn
  colnames(pr_matrix)[11:ncol(pr_matrix)] = cn

  ed = data.frame(label = cn,
                  condition = gsub('_\\d+$', '', cn),
                  replicate = gsub('.*_(\\d+)$', '\\1', cn))

  rownames(ed) = ed$label
  print(ed)
  return(list(pg_matrix = pg_matrix, pr_matrix = pr_matrix, expDesign = ed))
}

#' Prepare summarizedExperiment object from diann report.pg_matrix file
#'
#' @param pg_matrix the report.pg_matrix file from DIANN
#' @param expDesign A data frame with the experimental design. Should contain at least 'label', 'condition', and 'replicate' columns.
#' @param pr_matrix Optional argument. If the report.pr_matrix file from DIANN is provided, peptide information will be added to output.
#' @param missing_thr Integer specifying which proteinGroups are filtered out based on missing values.
#' @param min_peptides An integer specifing the cutoff for razor/unique peptides. The default is 0.
#' @param impute Specifies which imputatation method to use (default: knn). No imputation is done when entering 'none'.
#' See details for options.
#' @param mixed_cutoff Either 'empirally' or a value between 0-1. For details, see \link{mixed_imputation}
#' @param remove_contaminants A logical value specifying if potential
#' contaminants should be removed from the pg_matrix.
#'
#' @details
#' For standard imputation options, see ?DEP::impute. For mixed imputation, see \link{mixed_imputation}
#'
#'
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
#' # creates se without imputing missing values.
#' se <- prepare_se(report.pg_matrix,
#'                 expDesign,
#'                  missing_thr = 1,
#'                  impute = 'none')
prepare_se = function(pg_matrix, expDesign, pr_matrix = NULL, missing_thr = 0,
                      min_peptides = 0, impute = 'knn',
                      mixed_cutoff = 'empirically', remove_contaminants = TRUE){

  if (!is.null(pr_matrix) & !('n_total' %in% colnames(pg_matrix))){
    pep = get_nPep_prMatrix(pr_matrix)
    pg_matrix = add_peptide_numbers(pg_matrix, pep)
  }


  pg_matrix = add_contaminants(pg_matrix)

  if(remove_contaminants){pg_matrix = pg_matrix[pg_matrix$Potential.contaminant != '+',]}
  pg_uniq = DEP::make_unique(pg_matrix, 'Genes', 'Protein.Group', delim = ';')

  pat = paste(expDesign$label, collapse = '|')
  lfq = grep(pat, colnames(pg_uniq))
  se = DEP::make_se(pg_uniq, lfq, expDesign)


  if ('ibaq_peptides' %in% colnames(pg_matrix)){
    ibaq = pg_matrix[,grep('_iBAQ$', colnames(pg_matrix))]
    dimnames(ibaq) = list(rownames(se), colnames(se))
    SummarizedExperiment::rowData(se)$ibaq_peptides = pg_matrix$ibaq_peptides
    SummarizedExperiment::assay(se, 'iBAQ') = ibaq
    se = add_maxLFQ_iBAQ(se)
  }

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

  if(impute == 'mixed'){se = mixed_imputation(se, mixed_cutoff)}
  else if(!impute == 'none'){se = DEP::impute(se, impute)}

  colnames(se) = SummarizedExperiment::colData(se)$ID

  if(!is.null(pr_matrix)){
    colnames(pr_matrix)[11:ncol(pr_matrix)] = SummarizedExperiment::colData(se)$ID
    S4Vectors::metadata(se)$pr_matrix = pr_matrix
  }

  return(se)
}

#' Prepares summarizedExperiment object to work with functions from DEP package.
#'
#' @param se The summarizedExperiment object created with \link{prepare_se}
#'
#' @return A summarizedExperiment object.
#'
#' @import SummarizedExperiment
#'
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign, report.pr_matrix)
#' se_dep <- use_dep(se)
use_dep = function(se){
  colnames(se) = SummarizedExperiment::colData(se)$ID
  return(se)
}
