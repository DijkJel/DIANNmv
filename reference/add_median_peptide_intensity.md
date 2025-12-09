# Adds median peptide intensities to summerizedExperiment object

Adds median peptide intensities to summerizedExperiment object

## Usage

``` r
add_median_peptide_intensity(se)
```

## Arguments

- se:

  SummerizedExperiment object returned from DEP::make_se()

## Value

a summerizedExperiment object with peptide intensities as extra assay

## Examples

``` r
se <- prepare_se(report.pg_matrix, expDesign, report.pr_matrix)

#> Imputing along margin 2 (samples/columns).
#> [1] 0.3058978
#> Imputing along margin 1 (features/rows).
#> Warning: 36 rows with more than 50 % entries missing;
#>  mean imputation used for these rows
#> Cluster size 5511 broken into 3577 1934 
#> Cluster size 3577 broken into 2319 1258 
#> Cluster size 2319 broken into 3 2316 
#> Done cluster 3 
#> Cluster size 2316 broken into 1144 1172 
#> Done cluster 1144 
#> Done cluster 1172 
#> Done cluster 2316 
#> Done cluster 2319 
#> Done cluster 1258 
#> Done cluster 3577 
#> Cluster size 1934 broken into 1298 636 
#> Done cluster 1298 
#> Done cluster 636 
#> Done cluster 1934 
se <- add_median_peptide_intensity(se)
```
