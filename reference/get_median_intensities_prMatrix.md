# Calculates the median peptide intensity per proteinGroup per sample.

The median peptide intensity (MPI) can be used instead of iBAQ as proxy
for protein abundancy.

## Usage

``` r
get_median_intensities_prMatrix(
  pr_matrix,
  id_column = "Protein.Group",
  peptide_level = "Stripped.Sequence"
)
```

## Arguments

- pr_matrix:

  The report.pr_matrix file

- id_column:

  The column in pr_matrix that is used as identifier in the output file

- peptide_level:

  The column in pr_matrix that is used to aggregate peptide data

## Value

a matrix with median peptide intensities per sample per proteinGroup

## Examples

``` r
mpi <- get_median_intensities_prMatrix(report.pr_matrix)
```
