#' Combine different peptide variants into a single intensity
#'
#' @description
#' Peptide intensities can be split out over different entries based on
#' different charge states or modifications. This function combines these
#' variants into one intensity per peptide sequence.
#'
#' @param pr_matrix The report.pr_matrix file
#' @param peptide_level The column in pr_matrix that is used to aggregate.
#'
#' @import stats
#'
#' @return a matrix with a single peptide intensity per sample for each entry in the peptide_level column
#' @export
#'
#' @examples
#'
#' mat <- summarize_peptide_intensities(report.pr_matrix)
summarize_peptide_intensities = function(pr_matrix, peptide_level = 'Stripped.Sequence'){

  ints = as.matrix(pr_matrix[,11:ncol(pr_matrix)])
  ints[is.na(ints)] = 0
  ints = stats::aggregate(ints, list(peptide = pr_matrix[,peptide_level]), 'sum')
  rownames(ints) = ints[,1]
  ints = as.matrix(ints[,-1])
  ints[ints == 0] = NA
  ints = as.data.frame(ints)

  return(ints)
}


#' Get the numer of razor/unique peptides per proteinGroup per sample
#'
#' @param pr_matrix The report.pr_matrix file
#' @param id_column The IDs used in the output file
#' @param peptide_level The column specifying at which level intensities are
#' aggregated. See \link{summarize_peptide_intensities}
#'
#' @import stats
#'
#' @return a matrix with the count of razor/unique peptides
#' @export
#'
#' @examples
#' mat <- get_nPep_prMatrix(report.pr_matrix)
#'
get_nPep_prMatrix = function(pr_matrix, id_column = 'Protein.Group', peptide_level = 'Stripped.Sequence'){

  ints = summarize_peptide_intensities(pr_matrix, peptide_level)
  ints$protein = pr_matrix[match(rownames(ints), pr_matrix[,peptide_level]), id_column]

  npep = stats::aggregate(ints[,-ncol(ints)], list(protein = ints$protein), FUN = function(x){length(na.omit(x))})
  n_total = tapply(rownames(ints), ints$protein, length)
  npep$n_total = n_total[match(npep$protein, names(n_total))]


  #pr_matrixs = merge(pr_matrix, npep, by.x = id_column, by.y = 'protein', suffixes = c('_intensity', '_npep'))

  return(npep)
}

#' Calculate protein intensities from peptide information
#'
#' @param pr_matrix The report.pr_matrix file
#' @param id_column The IDs used in the output file
#' @param peptide_level he column specifying at which level intensities are
#' aggregated. See \link{summarize_peptide_intensities}
#'
#' @import stats
#'
#' @return A matrix with the summed peptides intenties of razor/unique peptides
#' @export
#'
#' @examples
#' mat <- get_intensities_prMatrix(report.pr_matrix)
get_intensities_prMatrix = function(pr_matrix, id_column = 'Protein.Group', peptide_level = 'Stripped.Sequence'){

  ints = summarize_peptide_intensities(pr_matrix, peptide_level)
  ints$protein = pr_matrix[match(rownames(ints), pr_matrix[,peptide_level]), id_column]
  npep = stats::aggregate(ints[,-ncol(ints)], list(protein = ints$protein), FUN = function(x){sum(na.omit(x))})

  return(npep)
}

#' Add the peptide number information to pg_matrix
#'
#' @param pg_matrix The report.pg_matrix file
#' @param peptide_numbers The output from \link{get_nPep_prMatrix}
#' @param id_column The column in pr_matrix that is used to match the
#' peptide_numbers and pg_matrix. Should be identical to the column used in
#' get_nPEP_prMatrix
#'
#' @return a pg_matrix data frame with added peptide number information
#' @export
#'
#' @examples
#' peptide_numbers <- get_nPep_prMatrix(report.pr_matrix)
#' pg_matrix <- add_peptide_numbers(report.pg_matrix, peptide_numbers)
add_peptide_numbers = function(pg_matrix, peptide_numbers, id_column = 'Protein.Group'){

  pg_matrix = merge(pg_matrix, peptide_numbers, by.x = id_column, by.y = 'protein', suffixes = c('', '_npep'))
  return(pg_matrix)
}


#' Calculates the median peptide intensity per proteinGroup per sample.
#'
#' @description
#' The median peptide intensity (MPI) can be used instead of iBAQ as proxy
#' for protein abundancy.
#'
#'
#' @param pr_matrix The report.pr_matrix file
#' @param id_column The column in pr_matrix that is used as identifier in the
#' output file
#' @param peptide_level The column in pr_matrix that is used to aggregate
#' peptide data
#'
#' @import stats
#'
#' @return a matrix with median peptide intensities per sample per proteinGroup
#' @export
#'
#' @examples
#' mpi <- get_median_intensities_prMatrix(report.pr_matrix)
get_median_intensities_prMatrix = function(pr_matrix, id_column = 'Protein.Group', peptide_level = 'Stripped.Sequence'){

  ints = summarize_peptide_intensities(pr_matrix, peptide_level)
  ints$protein = pr_matrix[match(rownames(ints), pr_matrix[,peptide_level]), id_column]

  pep_median = stats::aggregate(ints[,-ncol(ints)], list(protein = ints$protein), FUN = function(x){stats::median(x, na.rm = T)})
  return(pep_median)
}

#' Adds median peptide intensities to summerizedExperiment object
#'
#' @param se SummerizedExperiment object returned from DEP::make_se()
#' @param pr_matrix The report.pr_matrix file
#' @import SummarizedExperiment
#'
#' @return a summerizedExperiment object with peptide intensities as extra assay
#' @export
#'
#' @examples
#'
#' se <- prepare_se(report.pg_matrix, expDesign)
#' se <- add_median_peptide_intensity(se, report.pr_matrix)
add_median_peptide_intensity = function(se, pr_matrix){

  pi = get_median_intensities_prMatrix(pr_matrix)
  #pi = peptide_intensities
  pi = pi[pi[,1] %in% SummarizedExperiment::rowData(se)$Protein.Group,]
  pi = pi[match(SummarizedExperiment::rowData(se)$Protein.Group, pi[,1]),]
  rownames(pi) = rownames(se)
  pi = pi[,-1]

  if (!all(colnames(pi) == SummarizedExperiment::colData(se)$label)){
    colnames(pi) = colnames(se)
    warning('Colnames of peptide file do not match colnames SE. Same sample order for both is assumed.')
  }

  SummarizedExperiment::assay(se, 'median_peptide_intensities') = pi
  SummarizedExperiment::rowData(se)$baseMean_mpi = rowMeans(SummarizedExperiment::assay(se, 'median_peptide_intensities'), na.rm = T)
  return(se)
}



