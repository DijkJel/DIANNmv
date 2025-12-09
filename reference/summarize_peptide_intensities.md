# Combine different peptide variants into a single intensity

Peptide intensities can be split out over different entries based on
different charge states or modifications. This function combines these
variants into one intensity per peptide sequence.

## Usage

``` r
summarize_peptide_intensities(pr_matrix, peptide_level = "Stripped.Sequence")
```

## Arguments

- pr_matrix:

  The report.pr_matrix file

- peptide_level:

  The column in pr_matrix that is used to aggregate.

## Value

a matrix with a single peptide intensity per sample for each entry in
the peptide_level column

## Examples

``` r
mat <- summarize_peptide_intensities(report.pr_matrix)
```
