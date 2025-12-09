# Tidies sample names and parses an experimental design.

Tidies sample names and parses an experimental design.

## Usage

``` r
prepare_diann_data(pg_matrix, pr_matrix, no_samples = NULL)
```

## Arguments

- pg_matrix:

  The report.pg_matrix file

- pr_matrix:

  The report.pr_matrix file

- no_samples:

  Optional numeric value specifying the number of samples. If not
  provided, a educated guess will be made.

## Value

A list with a tidied pg_matrix, pr_matrix, and parsed expDesign.

## Details

This functions tidies sample names and prepares and experimental design,
assuming that the structure of the sample names is
\<prefix_condition_replicate\>. An example of the prepared sample names
and expDesign is printed in the console.

## Examples

``` r
tidy_data <- prepare_diann_data(report.pg_matrix, report.pr_matrix)
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: number of items to replace is not a multiple of replacement length
#>                 label condition replicate
#> neg_ctrl_3 neg_ctrl_3  neg_ctrl         3
#> motif1_1     motif1_1    motif1         1
#> motif1_2     motif1_2    motif1         2
#> motif1_3     motif1_3    motif1         3
#> motif2_1     motif2_1    motif2         1
#> motif2_2     motif2_2    motif2         2
#> motif2_3     motif2_3    motif2         3
pg_matrix <- tidy_data$pg_matrix
pr_matrix <- tidy_data$pr_matrix
```
