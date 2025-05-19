#' Perform in silico tryptic digest on uniprot fasta file
#'
#' @param fasta_location Path to fasta file
#' @param miscleaves Allowed number of miscleavages. Defaults to 0 (standard for iBAQ)
#' @param enzyme The enzyme used for digestion. Defaults to 'trypsin'. See 'cleaver'
#' package vignette for more options.
#'
#' @import cleaver stringr seqinr
#'
#' @return A data frame with the number of theoretically observable peptides for
#' each entry in the fasta file.
#' @export
#'
#' @examples
#' \dontrun{
#' no_ibaq_peptides <- get_ibaq_peptides('path/to/fasta.fasta')
#' }
get_ibaq_peptides = function(fasta_location, miscleaves = 0, enzyme = 'trypsin'){


  fasta = seqinr::read.fasta(fasta_location, seqtype = 'AA', as.string = TRUE, whole.header = F)
  mat = sapply(fasta, function(x){

    seq = x[1]
    name = attr(x, 'name')
    id = strsplit(name, '\\|')[[1]][2]
    npep = cleaver::cleave(seq, 'trypsin', missedCleavages = 0)[[1]]
    pep_lengths = sapply(npep, str_length)
    npep = npep[pep_lengths < 31 & pep_lengths > 6]
    npep = length(npep)

    c(seq, name, id, npep)
  })

  df = as.data.frame(t(mat))
  colnames(df) = c('sequence', 'name', 'uniprotID', 'npep')
  df$npep = as.numeric(df$npep)
  return(df)
}

#' Calculate iBAQ values from raw intensities
#'
#' @param pr_matrix The report.pr_matrix file
#' @param ibaq_stats Dataframe with proteinIDs and associated number of iBAQ peptides obtained with \link{get_ibaq_peptides}.
#' Required when not using pre-calculated iBAQ peptides
#' @param organism Specifies which organism to use if using pre-calculated iBAQ peptides. 'hs' for human, 'mm' for mouse.
#'
#' @return A dataframe with iBAQ values per proteinGroup per sample and the number
#' of iBAQ peptides.
#' @export
#'
#' @examples
#' ibaq_values <- calculate_iBAQ(report.pr_matrix, organism = 'hs') # For human
#' ibaq_values <- calculate_iBAQ(report.pr_matrix, organism = 'mm') # For mouse
calculate_iBAQ = function(pr_matrix, ibaq_stats = NULL, organism = 'hs'){

  if (is.null(ibaq_stats)){
    ibaq_peptides <- DIANNmv::ibaq_peptides
    ibaq_stats <- ibaq_peptides[[organism]]
  }

  protein_intensities = get_intensities_prMatrix(pr_matrix)
  df = merge(protein_intensities, ibaq_stats[,3:4], by.x = 'protein', by.y = 'uniprotID')
  ibaq_vals = df[,-c(1,ncol(df))] / df[,ncol(df)]
  colnames(ibaq_vals) = paste0(colnames(ibaq_vals), '_iBAQ')
  ibaq_vals$ibaq_peptides = df$npep
  rownames(ibaq_vals) = df$protein
  return(ibaq_vals)
}


#' Add iBAQ intensities to report.pg_matrix file.
#'
#' @param pg_matrix The report.pg_matrix file.
#' @param pr_matrix The report.pr_matrix file.
#' @param ibaq_stats Dataframe with proteinIDs and associated number of iBAQ peptides obtained with \link{get_ibaq_peptides}.
#' Required when not using pre-calculated iBAQ peptides
#' @param organism Specifies which organism to use if using pre-calculated iBAQ peptides. 'hs' for human, 'mm' for mouse.
#'
#' @return A report.pg_matrix data frame with iBAQ columns added.
#' @export
#'
#' @examples
#' pg <- add_iBAQ(report.pg_matrix, report.pr_matrix, organism = 'hs')
add_iBAQ = function(pg_matrix, pr_matrix, ibaq_stats = NULL, organism = 'hs'){

  ibaq_vals = calculate_iBAQ(pr_matrix, ibaq_stats, organism)
  pg = merge(pg_matrix, ibaq_vals, by.x = 'Protein.Group', by.y = 0)
  return(pg)
}
