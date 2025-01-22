#' Prepare summarizedExperiment object from diann report.pg_matrix file
#'
#' @param pg_matrix the report.pg_matrix file from DIANN
#' @param expDesign A data frame with the experimental design. Should contain at least 'label', 'condition', and 'replicate' columns
#' @param missing_thr Integer specifying which proteinGroups are filtered out based on missing values.
#' @param impute A logical value specifying if missing values should be imputed. Imputation is done using the 'MinProb' value from 'MsCoreUtils' package.
#'
#' @import DEP
#' @import grid
#'
#' @return A summarized Experiment object
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix,
#'                 expDesign, missing_thr = 1,
#'                 impute = TRUE) # creates se with missing values imputed
#'
#' se <- prepare_se(report.pg_matrix,
#'                 expDesign,
#'                  missing_thr = 1,
#'                  impute = FALSE) # creates se without imputing missing values
prepare_se = function(pg_matrix, expDesign, missing_thr = 0, impute = T){

  pg_uniq = DEP::make_unique(pg_matrix, 'Genes', 'Protein.Group', delim = ';')
  pat = paste(expDesign$label, collapse = '|')
  lfq = grep(pat, colnames(pg_uniq))
  se = DEP::make_se(pg_uniq, lfq, expDesign)
  se_filt = DEP::filter_missval(se, missing_thr)
  se = DEP::normalize_vsn(se_filt)
  p = DEP::plot_normalization(se_filt, se)
  grid::grid.newpage()
  grid::grid.draw(p)

  if(impute){se = DEP::impute(se, 'MinProb')}

  return(se)
}
