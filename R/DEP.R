#' Relabel the significance column in the DEP results data frame
#'
#' @param res Results data frame after DEP analysis.
#' @param pval_cutoff The p.adj significance cutoff
#' @param fc_cutoff The absolute log2 fold change significance cutoff
#'
#'
#' @return The results data frame with adjusted significance column
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign)
#' res <- get_DEPresults(se, type = 'all')
#' res <- recode_sig_col(res, pval_cutoff = 0.01, fc_cutoff = 1)
recode_sig_col = function(res, pval_cutoff = 0.05, fc_cutoff = 1){

  data = res
  ratio_col = grep('ratio', colnames(data), value = T)
  padj_col = grep('p.adj', colnames(data), value = T)
  sig_col = grep('sig', colnames(data), value = T)


  r_cols = as.list(data[,ratio_col])
  p_cols = as.list(data[,padj_col])

  sig_df = as.data.frame(mapply(function(x, y){abs(x) > fc_cutoff & y < pval_cutoff}, r_cols, p_cols))
  if (ncol(sig_df) == 1){colnames(sig_df) = gsub('ratio', 'significant', ratio_col)}
  sig_df$significant = apply(sig_df, 1, any)
  colnames(sig_df) = gsub('ratio', 'significant', colnames(sig_df))

  data[,sig_col] = sig_df
  return(data)
}


#' Perform differential protein expression analysis
#'
#' @param se SummarizedExperiment object with protein expression data.
#' @param condition1 One of the conditions to be tested present in expDesign.
#' @param condition2 One of the conditions to be tested present in expDesign.
#' @param ref_condition One of the conditions in expDesign that is used as
#' reference. Must be used in combination with type = 'control'.
#' @param alpha the p.adjust significance cutoff.
#' @param lfc The log2 fold change significance cutoff.
#' @param type Type of comparison to be made. Options are: 'manual', 'control',
#'  and 'all'.
#' @param tests Character vector specifying which contrasts to test when using
#' "type = 'manual'"
#' @param fdr_type Type of fdr correction. Options are 'fdrtool' and
#' 'BH' (Benjamini-Hochberg). Default is 'BH'
#' @param missing_thr Optional numeric value. Checks if one of the two conditions
#' in 1-vs-1 comparisons has enough values. Adds '_overimpute' column to output.
#'
#' @import DEP
#' @import SummarizedExperiment
#' @import stats
#'
#' @return A data frame with the results of differential protein expression
#' analysis
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign)
#' res <- get_DEPresults(se, 'motif1', 'neg_ctrl', type = 'manual')


get_DEPresults = function(se, condition1 = NULL, condition2 = NULL, ref_condition = NULL, tests = NULL,
                          alpha = 0.05, lfc = 1, type = 'manual', fdr_type = 'BH', missing_thr = NA){


  get_imputation_column = function(se, conditions = NULL){
    imputation_mask = SummarizedExperiment::assay(se, 'imputation_mask')
    conditions = paste0(conditions, '_')
    pat = paste(conditions, collapse = '|')
    if (!is.null(conditions)){imputation_mask = imputation_mask[,grep(pat, colnames(imputation_mask))]}
    is_imputed = apply(imputation_mask, 1, any)

    return(is_imputed)
  }

  get_missing_column = function(se, conditions = NULL, missing_thr = NA){

    imputation_mask = SummarizedExperiment::assay(se, 'imputation_mask')
    missing_mat = sapply(conditions, function(x){
      conditions = paste0(x, '_')
      imp = imputation_mask[,grep(conditions, colnames(imputation_mask))]
      rowSums(imp) > missing_thr
    })

    passes_missing = rowSums(missing_mat) == 2
    return(passes_missing)
  }

  if (type == 'manual'){
    if (!is.null(tests)){

      pats = lapply(tests, function(x){strsplit(x, '_vs_')[[1]]})
      imputation_mask = sapply(pats, function(x){get_imputation_column(se, x)})
      colnames(imputation_mask) = paste0(tests, '_imputed')

      if (!is.na(missing_thr)){
        missing_mask = sapply(pats, function(x){get_missing_column(se, x, missing_thr)})
        colnames(missing_mask) = paste0(tests, '_overimputed')
      }


      dep = DEP::test_diff(se, type = 'manual', test = tests)
      dep <- DEP::add_rejections(dep, alpha = alpha, lfc = lfc)
      res = DEP::get_results(dep)
      res = res[order(res$name),]

      pval_cols = grep('p.val', colnames(res))
      padj_cols = grep('p.adj', colnames(res))
      padj_names = grep('p.adj', colnames(res), value = T)

      if (fdr_type == 'BH'){
        padj_mat = apply(res[,pval_cols, drop = F], 2, function(x){stats::p.adjust(x, method = 'BH')})
        colnames(padj_mat) = colnames(res)[padj_cols]
        if (ncol(padj_mat) == 1){padj_mat = as.numeric(padj_mat)}
        res[,padj_cols] = padj_mat
        colnames(res)[padj_cols] = padj_names
        res = recode_sig_col(res, alpha, lfc)
      }
    }

    else {
      pat = paste(condition1, condition2, sep = '|')
      se = se[,grep(pat, colnames(se))]
      test = paste0(condition1, '_vs_', condition2)

      imputation_mask = as.matrix(get_imputation_column(se, c(condition1, condition2)), ncol = 1)
      colnames(imputation_mask) = paste0(test, '_imputed')

      if (!is.na(missing_thr)){
        missing_mask = as.matrix(get_missing_column(se, c(condition1, condition2), missing_thr), ncol = 1)
        colnames(missing_mask) = paste0(test, '_overimputed')
      }


      dep = DEP::test_diff(se, type = 'manual', test = test)
      dep <- DEP::add_rejections(dep, alpha = alpha, lfc = lfc)
      res = DEP::get_results(dep)
      res = res[order(res$name),]

      res$bh = stats::p.adjust(res[,grep('p.val', colnames(res))], method = 'BH')


      fc_col = grep('ratio', colnames(res))
      rd = as.data.frame(rowData(dep))
      diff_col = grep('diff', colnames(rd))
      res[,fc_col] = rd[,diff_col]

      if (fdr_type == 'BH'){
        res[,grep('p.adj', colnames(res))] = res$bh
        res = recode_sig_col(res, alpha, lfc)
      }

    }
  }

  else if (type == 'control'){

    conditions = unique(SummarizedExperiment::colData(se)$condition)
    conditions = conditions[conditions != ref_condition]
    conditions = lapply(conditions, c, ref_condition)
    imputation_mask = sapply(conditions, function(x){get_imputation_column(se, x)})
    cnames = sapply(conditions, paste, collapse = '_vs_')
    colnames(imputation_mask) = paste0(cnames, '_imputed')

    if (!is.na(missing_thr)){
      missing_mask = sapply(conditions, function(x){get_missing_column(se, x, missing_thr)})
      colnames(missing_mask) = paste0(cnames, '_overimputed')
    }


    dep = DEP::test_diff(se, type = "control", control = ref_condition)
    res = DEP::add_rejections(dep, alpha, lfc)
    res = DEP::get_results(res)
    res = res[order(res$name),]

    pval_cols = grep('p.val', colnames(res))
    padj_cols = grep('p.adj', colnames(res))

    if (fdr_type == 'BH'){
      padj_mat = apply(res[,pval_cols], 2, function(x){stats::p.adjust(x, method = 'BH')})
      res[,padj_cols] = padj_mat
      res = recode_sig_col(res, alpha, lfc)
    }

  }

  else if (type == 'all'){
    dep = DEP::test_diff(se, type = 'all')
    res = DEP::add_rejections(dep, alpha, lfc)
    res = DEP::get_results(res)

    pats = grep('ratio', colnames(res), value = T)
    pats = gsub('_ratio', '', pats)
    pats = lapply(pats, function(x){strsplit(x, '_vs_')[[1]]})
    imputation_mask = sapply(pats, function(x){get_imputation_column(se, x)})
    cnames = sapply(pats, paste, collapse = '_vs_')
    colnames(imputation_mask) = paste0(cnames, '_imputed')

    if (!is.na(missing_thr)){
      missing_mask = sapply(pats, function(x){get_missing_column(se, x, missing_thr)})
      colnames(missing_mask) = paste0(cnames, '_overimputed')
    }


    pval_cols = grep('p.val', colnames(res))
    padj_cols = grep('p.adj', colnames(res))

    if (fdr_type == 'BH'){
      padj_mat = apply(res[,pval_cols, drop = F], 2, function(x){stats::p.adjust(x, method = 'BH')})
      res[,padj_cols] = padj_mat
      res = recode_sig_col(res, alpha, lfc)
    }
  }

  if ('baseMean_mpi' %in% names(SummarizedExperiment::rowData(se))){
    mpi = assay(se, 'median_peptide_intensities')
    mpi = rowMeans(mpi, na.rm = T)
    res$baseMean_mpi = mpi
  }

  imputation_mask = imputation_mask[res$name,,drop = F]
  res = cbind(res, imputation_mask)

  if (!is.na(missing_thr)){
    res = cbind(res, missing_mask)
  }

  return(res)
}

