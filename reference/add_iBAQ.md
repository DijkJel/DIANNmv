# Add iBAQ intensities to report.pg_matrix file.

Add iBAQ intensities to report.pg_matrix file.

## Usage

``` r
add_iBAQ(pg_matrix, pr_matrix, ibaq_stats = NULL, organism = "hs")
```

## Arguments

- pg_matrix:

  The report.pg_matrix file.

- pr_matrix:

  The report.pr_matrix file.

- ibaq_stats:

  Dataframe with proteinIDs and associated number of iBAQ peptides
  obtained with
  [get_ibaq_peptides](https://dijkjel.github.io/DIANNmv/reference/get_ibaq_peptides.md).
  Required when not using pre-calculated iBAQ peptides

- organism:

  Specifies which organism to use if using pre-calculated iBAQ peptides.
  'hs' for human, 'mm' for mouse.

## Value

A report.pg_matrix data frame with iBAQ columns added.

## Examples

``` r
pg <- add_iBAQ(report.pg_matrix, report.pr_matrix, organism = 'hs')
```
