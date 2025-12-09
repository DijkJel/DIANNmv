# Add the peptide number information to pg_matrix

Add the peptide number information to pg_matrix

## Usage

``` r
add_peptide_numbers(pg_matrix, peptide_numbers, id_column = "Protein.Group")
```

## Arguments

- pg_matrix:

  The report.pg_matrix file

- peptide_numbers:

  The output from
  [get_nPep_prMatrix](https://dijkjel.github.io/DIANNmv/reference/get_nPep_prMatrix.md)

- id_column:

  The column in pr_matrix that is used to match the peptide_numbers and
  pg_matrix. Should be identical to the column used in get_nPEP_prMatrix

## Value

a pg_matrix data frame with added peptide number information

## Examples

``` r
peptide_numbers <- get_nPep_prMatrix(report.pr_matrix)
pg_matrix <- add_peptide_numbers(report.pg_matrix, peptide_numbers)
```
