#' Prepares data for the different GSEA visualization options
#'
#' @param gsea Output from \link{perform_GSEA}.
#' @param padj_cutoff The maximum p.adjust value allowed for inclusion of the pathway.
#' @param top_n The maximum number of pathways included. Takes the top_n pathways
#' with the lowest p.adj values.
#' @param remove_prefix Boolean specifying to remove the prefix from pathway names.
#' @param max_name_length Numeric value specifying the max length of pathway names.
#'
#' @import stringr
#'
#'
#' @return A data frame with filtered GSEA data.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' gsea <- perform_GSEA(res, genesets)
#' data <- prepare_gsea_data(gsea) # all significant pathways
#' data <- prepare_gsea_data(gsea, top_n = 10) # 10 pathways with lowest p.adjust values
#' (minimally padj < 0.05)
#' }
prepare_gsea_data = function(gsea, padj_cutoff = 0.05, top_n = Inf, remove_prefix = F, max_name_length = Inf){

  remove_every_other_vowel_except_start <- function(word) {
    chars <- stringr::str_split(word, "")[[1]]
    n <- length(chars)
    protected <- rep(FALSE, n)
    is_vowel <- stringr::str_detect(chars, "[aeiouAEIOU]")
    vowel_pos <- which(is_vowel)

    # Protect vowel at start of word (first character)
    if (n > 0 && is_vowel[1]) {
      protected[1] <- TRUE
    }

    # Protect vowels after underscore
    for (i in 2:n) {
      if (chars[i-1] == "_" && is_vowel[i]) {
        protected[i] <- TRUE
      }
    }

    # Get non-protected vowel positions
    non_protected_vowel_pos <- vowel_pos[!protected[vowel_pos]]

    # Keep every other non-protected vowel (remove every second)
    # e.g. keep 1st, 3rd, 5th; remove 2nd, 4th, etc.
    keep <- seq(1, length(non_protected_vowel_pos), by = 2)
    to_remove <- non_protected_vowel_pos[-keep]

    # Replace those vowels with empty strings
    chars[to_remove] <- ""
    paste(chars, collapse = "")
  }

  gsea = gsea[gsea$padj < padj_cutoff,]
  gsea = gsea[order(gsea$padj),]
  gsea = gsea[1:ifelse(is.infinite(top_n), nrow(gsea), min(top_n, nrow(gsea))),]

  if (remove_prefix){gsea$pathway = gsub("^[^_]*_", "", gsea$pathway)}

  pw_names = sapply(gsea$pathway, function(x){

    pw_len = nchar(x)
    pw_name = x
    if (nchar(pw_name) > max_name_length){
      pw_name = remove_every_other_vowel_except_start(pw_name)
    }
    if (nchar(pw_name) > max_name_length){
      pw_name = stringr::str_sub(pw_name, 1, max_name_length)
      pw_name = paste0(pw_name, '-')
    }

    return(pw_name)
  })

  gsea$pathway = pw_names

  return(gsea)
}

#' Plots GSEA results as bar plot.
#'
#' @param gsea Output from \link{perform_GSEA}.
#' @param pos_color Character vector specifying which color to use for activated pathways.
#' @param neg_color Character vector specifying which color to use for repressed pathways.
#' @param padj_cutoff The maximum p.adjust value allowed for inclusion of the pathway.
#' @param top_n The maximum number of pathways included. Takes the top_n pathways
#' with the lowest p.adj values.
#' @param remove_prefix Boolean specifying to remove the prefix from pathway names.
#' @param max_name_length Numeric value specifying the max length of pathway names.
#'
#' @import ggplot2
#'
#' @return A ggplot bar plot object
#' @export
#'
#' @examples
#' \dontrun{
#' gsea <- perform_GSEA(res, genesets)
#'
#' # Bar plot of 10 pathways with lowest padj values.
#' barplot <- plot_gsea_barplot(gsea, top_n = 10)
#'
#' #Bar plot of 10 pathways with lowest values and different colors for bars.
#' barplot <- plot_gsea_barplot(gsea, pos_color = 'red3', neg_color = 'dodgerblue', top_n = 10)
#' }
plot_gsea_barplot = function(gsea, pos_color = 'gold1', neg_color = 'darkblue', padj_cutoff = 0.05,
                             top_n = Inf, remove_prefix = F, max_name_length = Inf){

  data = prepare_gsea_data(gsea, padj_cutoff, top_n, remove_prefix, max_name_length)
  data$group = ifelse(data$NES < 0, 'Negative', 'Positive')

  barplot = ggplot2::ggplot(data, ggplot2::aes(x = reorder(.data[['pathway']], .data[['NES']]), y = .data[['NES']], fill = .data[['group']])) +
    ggplot2::geom_bar(stat = 'identity', width = 0.6) +
    ggplot2::geom_text(data = data[data$NES > 0,], ggplot2::aes(x = .data[['pathway']], y = 0, label = paste0(.data[['pathway']], '   '), hjust = 1), size = 3.5) +
    ggplot2::geom_text(data = data[data$NES < 0,], ggplot2::aes(x = .data[['pathway']], y = 0, label = paste0(.data[['pathway']], '   '), hjust = 0), size = 3.5) +
    ggplot2::scale_fill_manual(values = c(Negative = neg_color, Positive = pos_color)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 18, color = 'black'),
          axis.title.x = element_text(size = 18, color = 'black'),
          legend.position = 'none') +
    ggplot2::labs(x = 'Normalized Enrichment Score') +
    ggplot2::theme(panel.background = element_rect(colour = "black", linewidth=2)) +
    ggplot2::ylim(c(max(abs(data$NES)) * -1.1,max(abs(data$NES)) * 1.1)) +
    ggplot2::coord_flip()

  return(barplot)
}

#' Plots GSEA data as bubble plot. Mainly useful to compare multiple gsea comparisons.
#'
#' @param gsea Output from \link{perform_GSEA}.
#' @param ... Additional gsea outputs to be included in visualization.
#' @param sample_names Character vector specifying the comparisons that were made in the GSEA. Optional.
#' @param padj_cutoff padj_cutoff The maximum p.adjust value allowed for inclusion of the pathway.
#' @param top_n he maximum number of pathways included. Takes the top_n pathways
#' with the lowest p.adj values.
#' @param remove_prefix Boolean specifying to remove the prefix from pathway names.
#' @param max_name_length Numeric value specifying the max length of pathway names.
#'
#' @import ggplot2
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#'
#' gsea1 <- perform_GSEA(res1, genesets)
#' gsea2 <- perform_GSEA(res2, genesets)
#'
#' #Show data of a single gsea as one column
#' bubbleplot <- plot_gsea_bubbleplot(gsea1, top_n = 10, sample_names = 'Comparison_1')
#'
#' # Show data of multiple GSEAs as multiple columns.
#' bubbleplot <- plot_gsea_bubbleplot(gsea1, gsea2, top_n = 10, sample_names = c('test1', 'test2'))
#' }

plot_gsea_bubbleplot = function(gsea, ..., sample_names = NULL, padj_cutoff = 0.05,
                                top_n = Inf, remove_prefix = F, max_name_length = Inf){

  l = list(gsea, ...)
  if (!is.null(sample_names)){
    names(l) = sample_names
  } else {
    names(l) = paste0('sample_', 1:length(l))
  }

  l = lapply(l, function(x){prepare_gsea_data(gsea, padj_cutoff, top_n, remove_prefix, max_name_length)})

  poi = unique(unlist(lapply(l, function(x){x$pathway})))
  df = data.frame(pathway = poi)


  data_l = lapply(seq(l), function(x){
    data = l[[x]]
    df = merge(df, as.data.frame(data[,1:7]), all.x = T, by = 'pathway')
    df$sample = names(l)[x]
    return(df)
  })

  data = do.call(rbind, data_l)
  data$sample = factor(data$sample, levels = names(l))

  bubbleplot = ggplot2::ggplot(data, ggplot2::aes(x = .data[['sample']], y = reorder(.data[['pathway']], .data[['NES']]), color = .data[['NES']], size = -log10(.data[['padj']]))) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient2(low = 'dodgerblue4', mid = 'white', high = 'red3') +
    ggplot2::theme_minimal() +
    ggplot2::theme(aspect.ratio = 5/1) +
    ggplot2::xlab('') +
    ggplot2::ylab('')

  return(bubbleplot)
}


#' Plots GSEA data as a faceted dot plot showing NES, padj, and set size.
#'
#' @param gsea gsea Output from \link{perform_GSEA}.
#' @param padj_cutoff padj_cutoff The maximum p.adjust value allowed for inclusion of the pathway.
#' @param top_n The maximum number of pathways included. Takes the top_n pathways
#' with the lowest p.adj values.
#' @param remove_prefix Boolean specifying to remove the prefix from pathway names.
#' @param max_name_length Numeric value specifying the max length of pathway names.
#'
#' @import ggplot2
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' gsea <- perform_GSEA(res, genesets)
#' dotplot <- plot_gsea_dotplot(gsea, top_n = 10)
#' }
plot_gsea_dotplot = function(gsea, padj_cutoff = 0.05,
                             top_n = Inf, remove_prefix = F, max_name_length = Inf){

  data = data = prepare_gsea_data(gsea, padj_cutoff, top_n, remove_prefix, max_name_length)
  data$group = ifelse(data$NES < 0, 'Repressed', 'Activated')
  data$group = factor(data$group, levels = c('Repressed', 'Activated'))


  complex_dotplot = ggplot2::ggplot(data, ggplot2::aes(x = reorder(.data[['pathway']], .data[['NES']]), y = .data[['NES']], color = -log10(.data[['padj']]), size = .data[['size']])) +
    ggplot2::geom_point() +
    ggplot2::coord_flip () +
    ggplot2::facet_grid(~group, scales = 'free') +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.spacing = unit(2, 'lines'),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          axis.text = element_text(color = 'black'),
          axis.title = element_text(color = 'black'),
          plot.title = element_text(hjust = 0.5)) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(y = 'Normalized Enrichment Score', x = '')

  return(complex_dotplot)
}

#' Plots GSEA data as volcano plot
#'
#' @param gsea gsea Output from \link{perform_GSEA}.
#' @param padj_cutoff padj_cutoff The maximum p.adjust value allowed for inclusion of the pathway.
#' @param label Specifies which points to label. The default is 'sig', labeling
#' all significant points. Entering a value for top_n limits the labeling to
#' the top_n up- and top_n down-regulated proteins based on the p.adj.
#' When providing a vector with protein names, only those points are labeled.
#' @param top_n Specifies how many significant points to label.
#' @param remove_prefix Boolean specifying to remove the prefix from pathway names.
#' @param max_name_length Numeric value specifying the max length of pathway names.
#' @param up_color Character string specifying the color for significant points with positive NES.
#' @param down_color Character string specifying the color for significant points with positive NES.
#' @param ns_color Character string specifying the color for non-significant points.
#'
#' @import ggplot2 ggrepel
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' gsea <- perform_GSEA(res, genesets)
#' volcano <- plot_gsea_volcano(gsea, top_n = 10)
#' }
plot_gsea_volcano = function(gsea, padj_cutoff = 0.05, label = 'sig', top_n = NULL,
                             up_color = 'red3', down_color = 'dodgerblue', ns_color = 'grey70'){


  data = as.data.frame(gsea[,1:7])
  data$color = apply(data, 1, function(x){

    x = as.numeric(x[c(3,6)])
    if(x[1] > padj_cutoff){ns_color}
    else if (x[2] < 0){down_color}
    else if (x[2] > 0){up_color}
  })

  data$color = factor(data$color, levels = c(ns_color, down_color, up_color))
  data = data[order(data$color, decreasing = T),]

  data$significant = ifelse(data$color == 'grey70', F, T)
  if (label == 'sig'){data$label = ifelse(data$significant, data$pathway, '')}
  else if (label == ''){data$label = ''}
  else {data$label = data$pathway[data$pathway %in% label]}

  data = data[order(data$padj),]


  if(!is.null(top_n)){

  up_data = data[data$NES > 0,]
  down_data = data[data$NES < 0,]

  to_keep_up = up_data[1:min(top_n, nrow(up_data)),'pathway']
  to_keep_down = down_data[1:min(top_n, nrow(down_data)),'pathway']
  to_keep = c(to_keep_up, to_keep_down)

  data$label = ifelse(data$label %in% to_keep, data$label, '')
  }

  p = ggplot2::ggplot(data, ggplot2::aes(x = .data[['NES']], y = -log10(.data[['padj']]), label = label)) +
    ggplot2::geom_point(color = data$color) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), color = 'red', linetype = 'dashed', linewidth = 1) +
    ggrepel::geom_text_repel(max.overlaps = Inf, min.segment.length = 0.01) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = 'Normalized Enrichment Score')

  p = add_standardTheme(p)

}
