#' Returns boolean masks specifying MAR/MNAR imputation
#'
#' @param matrix A matrix with log-transformed intensities
#' @param cutoff The cutoff that specifies MAR vs MNAR. See details.
#'
#' @details
#' A value is considered as MNAR when the mean value of the replicates of a
#' condition is below a threshold. This can be a fixed value between 0-1. E.g.
#' '0.1' specifies that this threshold is the 10th percentile of all values over
#' the replicates of a condition.
#'
#' When set to 'empirically', this threshold is determined based on the data.
#' In this case, all intensities that are the only non-missing value within
#' replicates are collected, and the median value of this set is used as cutoff.
#'
#' @import stats
#'
#' @return A list with two matrices: one for MAR, one for MNAR.
#' @export
#'
#' @examples
#' \dontrun{
#' masks <- create_imputation_mask(data_matrix) # Default option.
#' masks <- create_imputation_mask(data_matrix, cutoff = 0.1) # Sets the
#' #MAR/MNAR cutoff at the 10th percentile of observed values
#' }
create_imputation_mask = function(matrix, cutoff = 'empirically'){

  if (cutoff == 'empirically'){cutoff = get_detection_limit(matrix)}
  else {cutoff = stats::quantile(as.numeric(matrix), cutoff, na.rm = T)}

  mnar_mask = apply(matrix, 1, function(x){

    if (all(is.na(x))){rep(T, length(x))}
    else if (mean(x, na.rm = T) < cutoff){is.na(x)}
    else (rep(F, length(x)))
  })

  mar_mask = apply(matrix, 1, function(x){

    if (all(is.na(x))){rep(F, length(x))}
    else if (mean(x, na.rm = T) >= cutoff){is.na(x)}
    else (rep(F, length(x)))
  })

  return(list(mnar = t(mnar_mask), mar = t(mar_mask)))
}


#' Perform mixed imputation over a matrix containing data of a single condition.
#'
#' @param matrix A matrix with intensities. Should contain all replicates of
#' a single condition.
#' @param matrix_masks Boolean masks returning from \link{create_imputation_mask}
#'
#' @import MsCoreUtils
#'
#' @return A matrix with all complete observations after imputation.
#' @export
#'
#' @examples
#' \dontrun{
#' masks <- create_imputation_mask(data_matrix)
#' data <- perform_mixed_imputation(data_matrix, masks)
#' }
perform_mixed_imputation = function(matrix, matrix_masks){

  mnar_matrix = MsCoreUtils::impute_MinProb(matrix)
  mar_matrix = MsCoreUtils::impute_knn(matrix)

  matrix[matrix_masks$mnar] = mnar_matrix[matrix_masks$mnar]
  matrix[matrix_masks$mar] = mar_matrix[matrix_masks$mar]

  return(matrix)
}

#' Find the cutoff for MAR/MNAR classification
#'
#' @param data A data matrix with replicates of a single condition
#'
#' @return A numeric value
#' @export
#'
#' @examples
#' \dontrun{
#' cutoff <- get_detection_limit(data)
#' }
get_detection_limit = function(data){

  single_values = apply(data, 1, function(x){
    if (sum(!is.na(x)) == 1){x[!is.na(x)]}
    else {NA}
  })

  single_values = median(single_values, na.rm = T)
  return(single_values)
}


#' Perform mixed imputation on a data matrix
#'
#' @param data A matrix with intensity values
#' @param conditions Character vector specifying the experimental conditions
#' @param cutoff The cutoff used for MAR/MNAR classification. See \link{create_imputation_mask}
#'
#' @return A matrix without missing values
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign, report.pr_matrix, impute = 'none')
#' data_missing <- as.matrix(SummarizedExperiment::assay(se)) # Intensity matrix with missing values
#' conditions <- unique(se$conditions) # The experimental conditions
#' data_full <- mixed_imputation_matrix(data_missing, conditions, cutoff = 'empirically')
mixed_imputation_matrix = function(data, conditions, cutoff = 'empirically'){

  data_split = lapply(conditions, function(x){data[,grep(x, colnames(data))]})
  masks = lapply(data_split, function(x){
    masks = create_imputation_mask(x, cutoff)
    #imputed_data = perform_mixed_imputation(x, masks)
  })

  mnar_mask = lapply(masks, '[[', 1)
  mnar_mask = do.call(cbind, mnar_mask)
  mar_mask = lapply(masks, '[[', 2)
  mar_mask = do.call(cbind, mar_mask)

  imputed_data = perform_mixed_imputation(data, list(mnar = mnar_mask, mar = mar_mask))

  #imputed_data = do.call(cbind, imputed_data)
  return(imputed_data)
}


#' Perform mixed imputation on summarizedExperiment object
#'
#' @param se The summarizedExperiment object
#' @param cutoff The cutoff used for MAR/MNAR classification. See \link{create_imputation_mask}
#'
#' @import SummarizedExperiment
#'
#' @return A summarizedExperiment object with complete cases after mixed imputation.
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign, report.pr_matrix, impute = 'none')
#' se <- mixed_imputation(se)
mixed_imputation = function(se, cutoff = 'empirically'){

  data = as.matrix(SummarizedExperiment::assay(se))
  conditions = unique(se$condition)

  data = mixed_imputation_matrix(data, conditions, cutoff)
  SummarizedExperiment::assay(se) = data

  return(se)
}

