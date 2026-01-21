#' Adds standard changes to ggplot theme
#'
#' @param plot a ggplot object
#'
#' @import ggplot2
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(var1 = 1:10, var2 = sample(1:100, 10))
#' p <- ggplot(df, aes(x = var1, y = var2)) + geom_point()
#' p <- add_standardTheme(p)
add_standardTheme  = function(plot){

  plot + ggplot2::theme(axis.text = element_text(color = 'black', size = 18),
               axis.title = element_text(color = 'black', size = 18))
}

#' Prepare data for volcano plots
#'
#' @param res Data frame with results from get_DEPresults.
#' @param pval_cutoff The p.adj significance cutoff.
#' @param fc_cutoff The log2 fold change significance cutoff.
#' @param label Specifies which points to label. The default is 'sig', labeling
#' all significant points. Entering a value for top_n limits the labeling to
#' the top_n up- and down-regulated proteins based on the p.adj.
#' When providing a vector with protein names, onlythose points are labeled.
#'
#' @param top_n Specifies how many significant points to label.
#' @param up_color Specifies the color of the significantly
#' up-regulated proteins.
#' @param down_color Specifies the color of the significantly down-regulated
#' proteins.
#' @param ns_color Specifies the color of non-significant proteins.
#' @param remove_overimputed Boolean value. If set to TRUE, proteins with too many
#' imputed values in 1-vs-1 comparisons are removed. Requires 'missing_thr' set
#' in \link{get_DEPresults}.
#'
#' @return A list with a data frame for each comparison present in res.
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign)
#' res <- get_DEPresults(se, type = 'all')
#' data <- prepare_volcano_data(res)
#'


prepare_volcano_data = function(res, pval_cutoff = 0.05, fc_cutoff = 1, remove_overimputed = F, label = 'sig', top_n = NULL,
                                up_color = 'red3', down_color = 'dodgerblue', ns_color = 'grey70'){

  if (remove_overimputed){
    if (!any(grepl('_overimputed', colnames(res)))){stop('Overimputed column not present in data. Run get_DEPresults() with missing_thr set.')}
  }

  res = recode_sig_col(res, pval_cutoff, fc_cutoff)
  samples = grep('ratio', colnames(res), value = TRUE)
  samples = gsub('_ratio', '', samples)

  data_list = lapply(samples, function(x){

    df = res[,c('name', grep(x, colnames(res), value = TRUE))]
    sig_col = grep('sig', colnames(df))
    ratio_col = grep('ratio', colnames(df))
    padj_col = grep('p.adj', colnames(df))
    impute_col = grep('_imputed', colnames(df))
    if (remove_overimputed){
      overimpute_col = grep('_overimputed', colnames(df))
      df = df[!df[,overimpute_col],]
    }

    df$color = up_color
    df$color = ifelse(df[,ratio_col] < 0, down_color, df$color)
    df$color = ifelse(df[,sig_col], df$color, ns_color)
    df$shape = ifelse(df[,impute_col], 1, 19)

    if (length(label) == 1 & 'sig' %in% label){
      df$label = ifelse(df[[sig_col]], df$name, '')

      if (!is.null(top_n)){
        data_sig = df[df[[sig_col]],]
        data_sig = data_sig[order(data_sig[,padj_col],decreasing = F),]
        sig_pos = data_sig[data_sig[[ratio_col]] > 0,]
        sig_neg = data_sig[data_sig[[ratio_col]] < 0,]

        labels_to_keep = c(sig_pos[1:top_n, 'name'], sig_neg[1:top_n,'name'])
        #labels_to_keep = data_sig[c(1:top_n, (nrow(data_sig) - (top_n -1)):nrow(data_sig)),'name']
        df$label = ifelse(df$name %in% labels_to_keep, df$label, '')
      }
    }
    else{
      df$label = ifelse(df[,'name'] %in% label, df[,'name'], '')
    }

    return(df)

  })

  names(data_list) = samples
  return(data_list)
}
#' create volcano plots for all comparisons present in the results file
#'
#' @param res A data frame with results from get_DEPresults()
#' @param pval_cutoff The p.adj significance cutoff
#' @param fc_cutoff The log2 fold change significance cutoff
#' @param label Specifies which points to label. The default is 'sig', labeling
#' all significant points. Entering a value for top_n limits the labeling to
#' the top_n up- and top_n down-regulated proteins based on the p.adj.
#' When providing a vector with protein names, only those points are labeled.
#' @param top_n Specifies how many significant points to label.
#' @param up_color Specifies the color of the significantly
#' up-regulated proteins.
#' @param down_color Specifies the color of the significantly down-regulated
#' proteins.
#' @param ns_color Specifies the color of non-significant proteins.
#' @param specify_imputed Boolean specifying whether proteins with imputed values
#' need to be indicated with open circles, versus closed circles for complete cases.
#' @param remove_overimputed Boolean value. If set to TRUE, proteins with too many
#' imputed values in 1-vs-1 comparisons are removed. Requires 'missing_thr' set
#' in \link{get_DEPresults}.
#' @import ggplot2 ggrepel
#'
#' @return A single ggplot object (1 comparison) or a list with ggplot objects.
#' @export
#'
#' @examples
#'
#' library(ggplot2)
#' library(ggrepel)
#'
#' se <- prepare_se(report.pg_matrix, expDesign)
#' res <- get_DEPresults(se, type = 'all')
#' vol <- plotVolcano(res, top_n = 10) # Labels the top 10 upregulated and
#' # top10 downregulated proteins based on fdr.
#' vol <- plotVolcano(res, label = c('SMAD3', 'SMAD4')) # Only labels SMAD3/4
#' vol <- plotVolcano(res, up_color = 'green', down_color = 'yellow') #
#' # Gives a very ugly volcano plot
plotVolcano = function(res, pval_cutoff = 0.05, fc_cutoff = 1, label = 'sig', top_n = NULL,
                       up_color = 'red3', down_color = 'dodgerblue', ns_color = 'grey70',
                       specify_imputed = T, remove_overimputed = F){

  data = prepare_volcano_data(res, pval_cutoff, fc_cutoff, remove_overimputed, label, top_n, up_color, down_color, ns_color)
  titles = names(data)

  plot_list = lapply(seq(data), function(x){

    title = titles[x]
    data = data[[x]]

    data$color = factor(data$color, levels = c(ns_color, down_color, up_color))
    data = data[order(data$color, decreasing = T),]

    ratio_col = grep('ratio', colnames(data), value = T)
    padj_col = grep('p.adj', colnames(data), value = T)
    impute_col = grep('_imputed', colnames(data), value = T)

    data$category = ifelse(data[,impute_col], 'imputed', 'not-imputed')

    if (specify_imputed){

      p = ggplot2::ggplot(data, aes(x = .data[[ratio_col]], y = -log10(.data[[padj_col]]), label = label)) +
        ggplot2::geom_point(ggplot2::aes(shape = .data[['category']]), color = data$color) +
        ggplot2::geom_hline(yintercept = -log10(pval_cutoff), color = 'red', linetype = 'dashed', linewidth = 1) +
        ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), color = 'red', linetype = 'dashed', linewidth = 1) +
        ggrepel::geom_text_repel(max.overlaps = Inf, min.segment.length = 0.01) +
        ggplot2::ggtitle(title) + ggplot2::theme_classic() +
        ggplot2::scale_shape_manual(values = c('imputed' = 21, 'not-imputed' = 19))
    }

    else {
      p = ggplot2::ggplot(data, aes(x = .data[[ratio_col]], y = -log10(.data[[padj_col]]), label = label)) +
        ggplot2::geom_point(color = data$color) +
        ggplot2::geom_hline(yintercept = -log10(pval_cutoff), color = 'red', linetype = 'dashed', linewidth = 1) +
        ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), color = 'red', linetype = 'dashed', linewidth = 1) +
        ggrepel::geom_text_repel(max.overlaps = Inf, min.segment.length = 0.01) +
        ggplot2::ggtitle(title) + ggplot2::theme_classic()
    }

    p = add_standardTheme(p)
    return(p)
  })

  names(plot_list) = titles
  if (length(plot_list) == 1){plot_list = plot_list[[1]]}
  return(plot_list)
}


#' Prepare the data for making MA plots
#'
#' @param res A data frame with results from get_DEPresults()
#' @param label Specifies which points to label. The default is 'sig', labeling
#' all significant points. Entering a value for top_n limits the labeling to
#' the top_n up- and top_n down-regulated proteins based on the p.adj.
#' When providing a vector with protein names, only those points are labeled.
#' @param ns_color Specifies the color of non-significant proteins.
#' @param up_color Specifies the color of the significantly
#' up-regulated proteins.
#' @param down_color Specifies the color of the significantly down-regulated
#' proteins.
#'
#'
#' @return A list with data frames.
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign, report.pr_matrix)
#' se <- add_median_peptide_intensity(se)
#' res <- get_DEPresults(se, type = 'all')
#' data <- prepare_MA_data(res)
prepare_MA_data = function(res, label = NULL, ns_color = 'grey70', up_color = 'red3', down_color = 'dodgerblue'){


  samples = grep('ratio', colnames(res), value = TRUE)
  samples = gsub('_ratio', '', samples)

  data_list = lapply(samples, function(x){

    df = res[,c('name', grep(x, colnames(res), value = TRUE), 'baseMean_mpi')]

    ratio_col = grep('ratio', colnames(df), value = T)
    sig_col = grep('sig', colnames(df))

    if(!is.null(label)){df$label = ifelse(df$name %in% label, df$name, '')}
    else{df$label = ''}

    df = df[order(res$significant),]
    df$color = up_color
    df$color = ifelse(df[,ratio_col] < 0, down_color, df$color)
    df$color = ifelse(df[,sig_col], df$color, ns_color)

    return(df)

  })

  names(data_list) = samples
  return(data_list)
}




#' Makes MA plots for each comparison in the results data frame
#'
#' @param res A data frame with results from get_DEPresults()
#' @param label Specifies which points to label. The default is 'sig', labeling
#' all significant points. Entering a value for top_n limits the labeling to
#' the top_n up- and top_n down-regulated proteins based on the p.adj.
#' When providing a vector with protein names, only those points are labeled.
#' @param ns_color Specifies the color of non-significant proteins.
#' @param up_color Specifies the color of the significantly
#' up-regulated proteins.
#' @param down_color Specifies the color of the significantly down-regulated
#' proteins.
#'
#' @import ggplot2 ggrepel
#'
#' @return A single ggplot object (1 comparison) or a list with ggplot objects.
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign, report.pr_matrix)
#' se <- add_median_peptide_intensity(se)
#' res <- get_DEPresults(se, type = 'all')
#' MAplots <- plot_MA(res)
plot_MA = function(res, label = NULL, ns_color = 'grey70', up_color = 'red3', down_color = 'dodgerblue'){

  data = prepare_MA_data(res, label, ns_color, up_color, down_color)

  plot_list = lapply(seq(data), function(x){

    title = names(data)[[x]]
    data = data[[x]]

    ratio = grep('ratio', colnames(data), value = T)

    p = ggplot2::ggplot(data, ggplot2::aes(x = log10(.data[['baseMean_mpi']]), y = .data[[ratio]], label = .data[['label']])) +
      ggplot2::geom_point(color = data[['color']]) +
      ggplot2::theme_classic() +
      ggrepel::geom_text_repel(max.overlaps = Inf, min.segment.length = 0.01,
                               box.padding = 0.5, fontface = 'bold', color = 'black')

    p = add_standardTheme(p)
    return(p)
  })

  names(plot_list) = names(data)
  if(length(plot_list) == 1){plot_list = plot_list[[1]]}

  return(plot_list)
}

#' Plots a Venn diagram showing overlapping significant proteins
#'
#' @param res A data frame with results from get_DEPresults()
#' @param ... Additional results data frames (optional)
#' @param comparisons A character vector specifying which comparisons
#' to include. The default 'all' considers all comparisons present in res.
#' @param colors The colors used for the Venn diagram.
#' @param names (optional). Specifies the labels for the different subsets.
#' Defaults to the comparison names in colnames(res)
#'
#' @import VennDiagram
#'
#' @return A Venn diagram plot
#' @export
#'
#' @examples
#'
#' se <- prepare_se(report.pg_matrix, expDesign)
#' res <- get_DEPresults(se, type = 'all')
#' plot_venn_diagram(res) # Default option
#' plot_venn_diagram(res, colors = c('green', 'blue', 'red')) # Uses non-default
#' # colors
#' plot_venn_diagram(res, comparisons = c('neg_ctrl_vs_motif1',
#'                                        'neg_ctrl_vs_motif2')) # Only
#'                                        # includes two comparisons
#'
#'
plot_venn_diagram = function(res, ..., comparisons = 'all',
                             colors = c('#b3e2cd', '#fdcdac', '#cbd5e8',
                                        '#f4cae4', '#e6f5c9'),
                             names = NULL){


  if (!('all' %in% comparisons)){
    pat = paste(comparisons, collapse = '|')
    res = res[,c(1, 2, grep(pat, colnames(res)))]
  }

  samples = grep('ratio', colnames(res), value = TRUE)
  samples = gsub('_ratio', '', samples)
  sig_cols = paste0(samples, '_significant')

  df = res[,sig_cols]
  rownames(df) = res$name

  venn_list = lapply(df, function(x){rownames(df)[x]})

  if (is.null(names)){names(venn_list) = samples}
  else {names(venn_list) = names}


  colors = colors[1:length(venn_list)]

  vd = VennDiagram::venn.diagram(venn_list, filename = NULL,
                                 disable.logging = T, fill = colors,
                                 fontfamily = 'sans', fontface = 'bold',
                                 cat.fontfamily = 'sans', lty = 0)

  grid::grid.newpage()
  grid::grid.draw(vd)
}

#' Plots upset plot showsing overlapping significant proteins between comparisons
#'
#' @param res A data frame with results from get_DEPresults()
#' @param ... Additional results data frames (optional)
#' @param comparisons A character vector specifying which comparisons
#' to include. The default 'all' considers all comparisons present in res.
#' @param names (optional). Specifies the labels for the different subsets.
#' Defaults to the comparison names in colnames(res)
#'
#' @details
#' The layout of the resulting upset plot can be easily edited to your liking.
#' Check the vignette of ggupset to see how individual components can be altered.
#'
#'
#' @import ggplot2 ggupset
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign)
#' res <- get_DEPresults(se, type = 'all')
#' vd <- plot_upset(res)
#'
plot_upset = function(res, ..., comparisons = 'all', names = NULL){


  if (!('all' %in% comparisons)){
    pat = paste(comparisons, collapse = '|')
    res = res[,c(1, 2, grep(pat, colnames(res)))]
  }

  samples = grep('ratio', colnames(res), value = TRUE)
  samples = gsub('_ratio', '', samples)
  sig_cols = paste0(samples, '_significant')

  df = res[,sig_cols]
  rownames(df) = res$name

  venn_list = lapply(df, function(x){rownames(df)[x]})

  if (is.null(names)){names(venn_list) = samples}
  else {names(venn_list) = names}

  total_genes = unique(unlist(venn_list))
  set_member = sapply(venn_list, function(x){total_genes %in% x})
  set_member = apply(set_member, 1, function(x){colnames(set_member)[x]})

  df = data.frame(gene = total_genes)
  df$member_list = set_member

  upset_plot = ggplot2::ggplot(df, ggplot2::aes(x = .data[['member_list']])) +
    ggplot2::geom_bar() +
    ggupset::scale_x_upset() +
    ggplot2::theme_classic()

  upset_plot = add_standardTheme(upset_plot)

  return(upset_plot)
}


#' Prepare data for plot_DEP_barplot function.
#'
#' @param res A data frame with results from get_DEPresults()
#' @param comparisons A character vector specifying which comparisons
#' to include. The default 'all' considers all comparisons present in res.
#' The bars are plotted in the order of specified comparisons.
#' @param names (optional). Specifies the labels for the different bars.
#' Defaults to the comparison names in colnames(res)
#'
#' @import reshape2
#'
#' @return A data frame with the number of significant proteins per condition.
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign)
#' res <- get_DEPresults(se, type  = 'all')
#' data <- prepare_barplot_data(res) # All conditions
#' data <- prepare_barplot_data(res, comparisons = c('neg_ctrl_vs_motif1',
#'                                                   'neg_ctrl_vs_motif2')) #
#' # Only two comparisons included.
#'
prepare_barplot_data = function(res, comparisons = 'all', names = NULL){

  if (!('all' %in% comparisons)){
    pat = paste(comparisons, collapse = '|')
    res = res[,c(1, 2, grep(pat, colnames(res)))]
  }

  samples = grep('ratio', colnames(res), value = TRUE)
  samples = gsub('_ratio', '', samples)
  sig_cols = paste0(samples, '_significant')
  ratio_cols = paste0(samples, '_ratio')

  df_sig = as.list(res[,sig_cols, drop = F])
  df_ratio = as.list(res[,ratio_cols, drop = F])

  up_list = mapply(function(x, y){rownames(res)[x & (y > 0)]}, df_sig, df_ratio)
  down_list = mapply(function(x, y){rownames(res)[x & (y < 0)]}, df_sig, df_ratio)

  if (is.list(up_list)){up = lengths(up_list)} else {up = length(up_list)}
  if (is.list(down_list)){down = lengths(down_list)} else {down = length(down_list)}

  df = data.frame(Group = samples,
                  up = up,
                  down = down)

  if (!('all' %in% comparisons)){
    df$Group = factor(df$Group, levels = comparisons)
    df = df[order(df$Group),]
  }

  if(!is.null(names)){df$Group = factor(names, levels = names)}
  df_long = reshape2::melt(df, value.name = 'DEP', id.vars = 'Group', variable.name = 'direction')

  return(df_long)
}

#' Make bar plot showing number of differential proteins between conditions.
#'
#' @param res A data frame with results from get_DEPresults()
#' @param comparisons A character vector specifying which comparisons
#' to include. The default 'all' considers all comparisons present in res.
#' The bars are plotted in the order of specified comparisons.
#' @param names (optional). Specifies the labels for the different bars.
#' Defaults to the comparison names in colnames(res)
#'
#' @import ggplot2
#'
#' @return A ggplot bar plot
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign)
#' res <- get_DEPresults(se, type  = 'all')
#'
#' plot_DEP_barplot(res)
#' plot_DEP_barplot(res, comparisons = c('neg_ctrl_vs_motif1',
#'                                        'neg_ctrl_vs_motif2'))
#'
#'
plot_DEP_barplot = function(res, comparisons = 'all', names = NULL){

  data = prepare_barplot_data(res, comparisons, names)

  p = ggplot2::ggplot(data, ggplot2::aes(x = .data[['Group']], y = .data[['DEP']], fill = .data[['direction']])) +
    ggplot2::geom_bar(stat = 'identity') +
    ggplot2::geom_text(ggplot2::aes(label = .data[['DEP']]), position = position_stack(vjust = 0.5), color = 'white', size = 8) +
    ggplot2::scale_fill_manual(values = c('up' = 'red3', 'down' = 'dodgerblue')) +
    ggplot2::labs(x = '') +
    ggplot2::theme_classic()

  p = add_standardTheme(p)
  return(p)
}


