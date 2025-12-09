# Calculate protein intensities from peptide information

Calculate protein intensities from peptide information

## Usage

``` r
get_intensities_prMatrix(
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

  he column specifying at which level intensities are aggregated. See
  [summarize_peptide_intensities](https://dijkjel.github.io/DIANNmv/reference/summarize_peptide_intensities.md)

## Value

A matrix with the summed peptides intenties of razor/unique peptides

## Examples

``` r
mat <- get_intensities_prMatrix(report.pr_matrix)
```
