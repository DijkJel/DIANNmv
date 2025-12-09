# Prepares the data for protein coverage plots.

Prepares the data for protein coverage plots.

## Usage

``` r
prepare_peptide_data(
  pr_matrix,
  genes,
  positions = NULL,
  zoom = NULL,
  fasta = NULL,
  organism = "hs",
  combine_overlap = T,
  dodge_labels = T
)
```

## Arguments

- pr_matrix:

  The report.pr_matrix file.

- genes:

  The gene names of the proteins that you would like to plot.

- positions:

  A numeric vector indiciting the amino acid positions that you want to
  highlight.

- zoom:

  A numeric with two values (start and stop) indicating the specific
  part of the protein that you want to plot.

- fasta:

  Data frame with protein information when not using a default option.
  See
  [get_ibaq_peptides](https://dijkjel.github.io/DIANNmv/reference/get_ibaq_peptides.md)

- organism:

  Specifies which default protein database to use. 'hs' for human, or
  'mm' for mouse.

- combine_overlap:

  Boolean specifying whether to combine miscleaved peptides.

- dodge_labels:

  Boolean specifying wheter position labels should be plotted on
  different heights so that they do not overlap in the plot.

## Value

A list containing a data frame with information about the tryptic
peptides, and a data frame with the AA positions that are highlighted.

## Examples

``` r
coverage_data <- prepare_peptide_data(report.pr_matrix,  'SMAD4', positions = c(100, 120))
```
