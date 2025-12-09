# Calculate iBAQ values from raw intensities

Calculate iBAQ values from raw intensities

## Usage

``` r
calculate_iBAQ(pr_matrix, ibaq_stats = NULL, organism = "hs")
```

## Arguments

- pr_matrix:

  The report.pr_matrix file

- ibaq_stats:

  Dataframe with proteinIDs and associated number of iBAQ peptides
  obtained with
  [get_ibaq_peptides](https://dijkjel.github.io/DIANNmv/reference/get_ibaq_peptides.md).
  Required when not using pre-calculated iBAQ peptides

- organism:

  Specifies which organism to use if using pre-calculated iBAQ peptides.
  'hs' for human, 'mm' for mouse.

## Value

A dataframe with iBAQ values per proteinGroup per sample and the number
of iBAQ peptides.

## Examples

``` r
ibaq_values <- calculate_iBAQ(report.pr_matrix, organism = 'hs') # For human
ibaq_values <- calculate_iBAQ(report.pr_matrix, organism = 'mm') # For mouse
```
