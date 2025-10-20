#' Find the amino acid locations of peptides within a protein sequence
#'
#' @param peptide_sequences A character vector with the peptide sequences.
#' @param full_sequence A character vector with the sequence of the parent protein.
#'
#' @import stringr
#'
#' @return A data frame with the peptide sequences, their start position and end
#' position in the parent protein.
#' @export
#'
#' @examples
#'
#' pep_sequence <- 'AARTGGGPLR'
#' full_sequence <- 'ISIFQKRSRWCAARTGGGPLRDPEIYLWWW'
#' positions <- get_peptide_position_in_protein(pep_sequence, full_sequence)
#'
get_peptide_position_in_protein = function(peptide_sequences, full_sequence){

  start_stop = sapply(peptide_sequences, function(x){stringr::str_locate(full_sequence, x)})
  start_stop = as.data.frame(t(start_stop))
  colnames(start_stop) = c('start', 'end')
  start_stop$seq_length = start_stop[,2] - start_stop[,1]
  out = cbind(sequence = peptide_sequences, start_stop)
}


#' Prepares the data for protein coverage plots.
#'
#' @param pr_matrix The report.pr_matrix file.
#' @param genes The gene names of the proteins that you would like to plot.
#' @param positions A numeric vector indiciting the amino acid positions that
#' you want to highlight.
#' @param zoom A numeric with two values (start and stop) indicating the specific
#' part of the protein that you want to plot.
#' @param fasta Data frame with protein information when not using a default option.
#' See \link{get_ibaq_peptides}
#' @param organism Specifies which default protein database to use.
#' 'hs' for human, or 'mm' for mouse.
#' @param combine_overlap Boolean specifying whether to combine miscleaved peptides.
#' @param dodge_labels Boolean specifying wheter position labels should be plotted
#' on different heights so that they do not overlap in the plot.
#'
#' @import cleaver reshape2 ivs stats stringr
#'
#' @return A list containing a data frame with information about the tryptic
#' peptides, and a data frame with the AA positions that are highlighted.
#' @export
#'
#' @examples
#' coverage_data <- prepare_peptide_data(report.pr_matrix,  'SMAD4', positions = c(100, 120))
#'
prepare_peptide_data = function(pr_matrix, genes, positions = NULL, zoom = NULL, fasta = NULL, organism = 'hs', combine_overlap = T, dodge_labels = T){

  fasta_list = DIANNmv::ibaq_peptides
  if (is.null(fasta)){
    if(organism == 'hs'){fasta = fasta_list$hs}
    else {fasta = fasta_list$mm}
  }


  intensities = summarize_peptide_intensities(pr_matrix)
  intensities = as.data.frame(intensities)
  scaled_intensities = intensities
  scaled_intensities$theoretical = NA
  scaled_intensities = t(scale(t(log10(scaled_intensities)), scale = F))

  intensities$gene = pr_matrix[match(rownames(intensities), pr_matrix$Stripped.Sequence), 'Genes']
  intensities$uniprotID = pr_matrix[match(rownames(intensities), pr_matrix$Stripped.Sequence), 'Protein.Group']
  intensities = intensities[intensities$gene %in% genes,]

  intensities$sequence = rownames(intensities)

  df_list = split(intensities, intensities$gene)

  locations = lapply(df_list, function(x){
    upid =  x[1,'uniprotID']
    gene = pr_matrix[match(upid, pr_matrix$Protein.Group),'Genes']
    protein_sequence = fasta[fasta$uniprotID == x[1,'uniprotID'],'sequence']
    ibaq_sequences = cleaver::cleave(protein_sequence, enzym = 'trypsin')[[1]]

    ibaq_sequences = ibaq_sequences[nchar(ibaq_sequences) >= 7 & nchar(ibaq_sequences) <= 30]

    df = data.frame(sequence = ibaq_sequences,
                    intensity = NA,
                    centered_intensity = NA,
                    sample = 'theoretical',
                    gene = gene,
                    uniprotID = upid)

    df_locations = get_peptide_position_in_protein(df$sequence, protein_sequence)
    df = merge(df, df_locations, by = 'sequence')

    locations = get_peptide_position_in_protein(x$sequence, protein_sequence)
    out = merge(x, locations, by.x = 'sequence', by.y = 'sequence', all = T)

    out = reshape2::melt(out, id.vars = c('sequence','uniprotID', 'gene', 'start', 'end', 'seq_length'), value.name = 'intensity', variable.name = 'sample')

    if (combine_overlap){

      out$interval = ivs::iv(start = out$start, end = out$end)
      out$group = as.factor(format(ivs::iv_identify_group(out$interval)))
      out_agg = stats::aggregate(intensity ~ group + sample, data = out, FUN = sum, na.action = NULL, na.rm = T)
      pep_ranges = as.character(out_agg$group)
      pep_ranges = gsub('\\[|\\,|\\)', '', pep_ranges)
      out_agg$start = as.numeric(sapply(pep_ranges, function(x){strsplit(x, ' ')[[1]][1]}))
      out_agg$end = as.numeric(sapply(pep_ranges, function(x){strsplit(x, ' ')[[1]][2]}))
      out_agg$seq_length = out_agg$end - out_agg$start
      out_agg$full_length = nchar(protein_sequence)
      out_agg$sequence = apply(out_agg[,c('start', 'end')], 1, function(x){stringr::str_sub(protein_sequence, x[1], x[2])})
      out_agg$gene = gene
      out_agg$uniprotID = upid
      out_agg = out_agg[,-1]
      out = out_agg[,c('sequence','uniprotID', 'gene', 'start', 'end', 'seq_length', 'sample', 'intensity')]
    }


    df = df[,c('sequence','uniprotID', 'gene', 'start', 'end', 'seq_length', 'sample', 'intensity', 'centered_intensity')]

    out[,'intensity'][out[,'intensity'] == 0] = NA
    out$intensity = log10(out$intensity)
    int_wide = reshape2::acast(out, sequence ~ sample, value.var = 'intensity')
    int_center = t(scale(t(int_wide), scale = F))
    int_center_long = reshape2::melt(int_center, value.name = 'centered_intensity')
    out = merge(out, int_center_long, by.x = c('sequence', 'sample'), by.y = c('Var1', 'Var2'), all.x = T, all.y = F)

    out = rbind(out, df)
    out$full_length_start = 1
    out$full_length_end = nchar(protein_sequence)
    return(out)
  })

  locations = do.call(rbind, locations)
  locations[locations$sample == 'theoretical', c('intensity', 'centered_intensity')] = 0
  locations$condition = gsub('_\\d+$', '', locations$sample)
  locations$x_start = as.numeric(as.factor(locations$sample))
  locations$x_end = as.numeric(as.factor(locations$sample)) + 0.3

  if (!is.null(positions)){

    if(length(genes) > 1){
      warning('AA positions only work when a single gene name is provided.')
    }
    else{
      protein_sequence = fasta[fasta$uniprotID == locations[1, 'uniprotID'], 'sequence']
      aa = sapply(positions, function(x){stringr::str_sub(protein_sequence, x, x)})
      aa_labels = data.frame(position = positions, aa = aa)
      aa_labels$label = paste0(aa_labels$aa, aa_labels$position)

      if (dodge_labels){
        aa_labels$x_position = rep(c(-0.25, -0.5, -0.75), length.out = nrow(aa_labels))
      }
      else{
        aa_labels$x_position = -0.5
      }
    }
  }
  else{
    aa_labels = NULL
  }

  if (!is.null(zoom)){

    min_border = min(zoom)
    max_border = max(zoom)

    bool_filter = apply(locations[,c('start', 'end')], 1, function(x){
      start = x[1]
      end = x[2]
      if (end <= max_border & start >= min_border){T}
      else if (min_border < end & min_border > start){T}
      else if (max_border < end & max_border > start){T}
      else {F}
    })

    locations = locations[bool_filter,]
    locations$full_length_start = min(locations$start)
    locations$full_length_end = max(locations$end)

  }

  out_list = list(data = locations, aa_labels = aa_labels)

  return(out_list)
}

#' Plots the coverage of proteins.
#' @param se The SummarizedExperiment object from \link{prepare_se}.
#' @param genes The gene names of the proteins that you would like to plot.
#' @param positions A numeric vector indiciting the amino acid positions that
#' you want to highlight.
#' @param zoom A numeric with two values (start and stop) indicating the specific
#' part of the protein that you want to plot.
#' @param fasta Data frame with protein information when not using a default option.
#' See \link{get_ibaq_peptides}
#' @param organism Specifies which default protein database to use.
#' 'hs' for human, or 'mm' for mouse.
#' @param combine_overlap Boolean specifying whether to combine miscleaved peptides.
#' @param dodge_labels Boolean specifying wheter position labels should be plotted
#' on different heights so that they do not overlap in the plot.
#' @param scaling Boolean value indicating whether the peptide intensities should
#' be centered over the different replicates.
#' @param condition_order Character vector that specifies the order of the conditions
#' on the y-axis.
#'
#' @import ggplot2 S4Vectors
#'
#' @return A ggplot2 object showing the found peptides for a protein.
#' @export
#'
#' @examples
#' se <- prepare_se(report.pg_matrix, expDesign, report.pr_matrix)
#' smad4 <- plot_protein_coverage(se, 'SMAD4', positions = c(100, 150))
#'
#' # Zoom in on first 200 AA of protein
#' smad4 <- plot_protein_coverage(se, 'SMAD4', positions = c(100, 150), zoom = c(1, 200))
plot_protein_coverage = function(se, genes, positions = NULL, zoom = NULL, fasta = NULL, organism = 'hs', combine_overlap = T, dodge_labels = T, scaling = 'centered', condition_order = NULL){

  pr_matrix = S4Vectors::metadata(se)$pr_matrix
  rd = as.data.frame(rowData(se))
  pr_names = rd[match(pr_matrix$Protein.Group, rd$Protein.Group), 'name']
  pr_matrix$Genes = pr_names

  data_list = prepare_peptide_data(pr_matrix, genes, positions, zoom, fasta, organism, combine_overlap, dodge_labels)
  data = data_list$data
  #data = na.omit(data)
  aa_labels = data_list$aa_labels

  full_length = data[!duplicated(paste0(data$sample, data$gene)), c('sample', 'gene', 'x_start', 'x_end', 'full_length_start', 'full_length_end')]
  full_length$x_end = full_length$x_end - 0.2

  if (!is.null(condition_order)){
    ed = as.data.frame(colData(se))
    ed_order = sapply(condition_order, function(x){grep(x, ed$condition)})
    sample_order = ed[ed_order, 'ID']
    data$sample = factor(data$sample, levels = c('theoretical', sample_order))
  }

  if (scaling == 'centered'){

    p = ggplot2::ggplot(data[!is.na(data$centered_intensity),], ggplot2::aes(x = sample)) +
      ggplot2::geom_rect(data = full_length, ggplot2::aes(xmin = .data[['x_start']]-((.data[['x_end']] - .data[['x_start']])/2), xmax = .data[['x_end']]-((.data[['x_end']] - .data[['x_start']])/2), ymin = 1, ymin = .data[['full_length_start']], ymax = .data[['full_length_end']]), fill = 'grey')+
      ggplot2::geom_rect(ggplot2::aes(xmin = .data[['x_start']]-((.data[['x_end']] - .data[['x_start']])/2), xmax = .data[['x_end']]-((.data[['x_end']] - .data[['x_start']])/2), ymin = start, ymax = end, fill = .data[['centered_intensity']]), alpha = 1) +
      ggplot2::geom_rect(data = data[data$sample == 'theoretical',], ggplot2::aes(xmin = .data[['x_start']]-((.data[['x_end']] - .data[['x_start']])/2), xmax = .data[['x_end']]-((.data[['x_end']] - .data[['x_start']])/2), ymin = start, ymax = end), fill = 'black') +
      ggplot2::theme_classic() +
      ggplot2::facet_grid(~gene, scales = 'free') +
      ggplot2::scale_fill_gradient2(low = 'dodgerblue', mid = 'white', high = 'red3') +
      ggplot2::scale_x_discrete(limits = levels(data$sample)) +
      ggplot2::coord_flip()
  }
  else{
    min_value = min(data[data$sample != 'theoretical', 'intensity'], na.rm = T)

    p = ggplot2::ggplot(data[!is.na(data$intensity),], ggplot2::aes(x = sample)) +
      ggplot2::geom_rect(data = full_length, ggplot2::aes(xmin = .data[['x_start']]-((.data[['x_end']] - .data[['x_start']])/2), xmax = .data[['x_end']]-((.data[['x_end']] - .data[['x_start']])/2), ymin = .data[['full_length_start']], ymax = .data[['full_length_end']]), fill = 'grey')+
      ggplot2::geom_rect(ggplot2::aes(xmin = .data[['x_start']]-((.data[['x_end']] - .data[['x_start']])/2), xmax = .data[['x_end']]-((.data[['x_end']] - .data[['x_start']])/2), ymin = start, ymax = end, fill = .data[['intensity']]), alpha = 1) +
      ggplot2::geom_rect(data = data[data$sample == 'theoretical',], ggplot2::aes(xmin = .data[['x_start']]-((.data[['x_end']] - .data[['x_start']])/2), xmax = .data[['x_end']]-((.data[['x_end']] - .data[['x_start']])/2), ymin = start, ymax = end), fill = 'black') +
      ggplot2::theme_classic() +
      ggplot2::facet_grid(~gene, scales = 'free') +
      ggplot2::scale_fill_viridis_c(limits = c(min_value, NA)) +
      ggplot2::scale_x_discrete(limits = levels(data$sample)) +
      ggplot2::coord_flip()
  }

  if (!is.null(aa_labels)){
    p = p +
      ggplot2::geom_hline(yintercept = aa_labels$position, color = 'black', linewidth = 0.5, linetype = 'dashed') +
      ggplot2::geom_text(data = aa_labels, ggplot2::aes(y = .data[['position']], x = .data[['x_position']], label = .data[['label']]), ) +
      ggplot2::coord_flip(xlim = c(1, length(unique(data$sample))), clip = 'off') +
      #coord_cartesian(xlim = c(1, length(unique(data$sample))), clip = 'off') +
      ggplot2::theme(plot.margin = ggplot2::margin(b = 6, unit = "lines"),
                     axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 35)))
  }

  p = add_standardTheme(p)
  p = p + ggplot2::labs(x = 'Sample', y = 'Position')

  return(p)
}
