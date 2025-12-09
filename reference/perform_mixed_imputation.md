# Perform mixed imputation over a matrix containing data of a single condition.

Perform mixed imputation over a matrix containing data of a single
condition.

## Usage

``` r
perform_mixed_imputation(matrix, matrix_masks)
```

## Arguments

- matrix:

  A matrix with intensities. Should contain all replicates of a single
  condition.

- matrix_masks:

  Boolean masks returning from
  [create_imputation_mask](https://dijkjel.github.io/DIANNmv/reference/create_imputation_mask.md)

## Value

A matrix with all complete observations after imputation.

## Examples

``` r
if (FALSE) { # \dontrun{
masks <- create_imputation_mask(data_matrix)
data <- perform_mixed_imputation(data_matrix, masks)
} # }
```
