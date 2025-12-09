# Returns boolean masks specifying MAR/MNAR imputation

Returns boolean masks specifying MAR/MNAR imputation

## Usage

``` r
create_imputation_mask(matrix, cutoff = "empirically")
```

## Arguments

- matrix:

  A matrix with log-transformed intensities

- cutoff:

  The cutoff that specifies MAR vs MNAR. See details.

## Value

A list with two matrices: one for MAR, one for MNAR.

## Details

A value is considered as MNAR when the mean value of the replicates of a
condition is below a threshold. This can be a fixed value between 0-1.
E.g. '0.1' specifies that this threshold is the 10th percentile of all
values over the replicates of a condition.

When set to 'empirically', this threshold is determined based on the
data. In this case, all intensities that are the only non-missing value
within replicates are collected, and the median value of this set is
used as cutoff.

## Examples

``` r
if (FALSE) { # \dontrun{
masks <- create_imputation_mask(data_matrix) # Default option.
masks <- create_imputation_mask(data_matrix, cutoff = 0.1) # Sets the
#MAR/MNAR cutoff at the 10th percentile of observed values
} # }
```
