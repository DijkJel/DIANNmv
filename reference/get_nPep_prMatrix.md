# Get the numer of razor/unique peptides per proteinGroup per sample

Get the numer of razor/unique peptides per proteinGroup per sample

## Usage

``` r
get_nPep_prMatrix(
  pr_matrix,
  id_column = "Protein.Group",
  peptide_level = "Stripped.Sequence"
)
```

## Arguments

- pr_matrix:

  The report.pr_matrix file

- id_column:

  The IDs used in the output file

- peptide_level:

  The column specifying at which level intensities are aggregated. See
  [summarize_peptide_intensities](https://dijkjel.github.io/DIANNmv/reference/summarize_peptide_intensities.md)

## Value

a matrix with the count of razor/unique peptides

## Examples

``` r
mat <- get_nPep_prMatrix(report.pr_matrix)
```
