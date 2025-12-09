# Make bar plot showing number of differential proteins between conditions.

Make bar plot showing number of differential proteins between
conditions.

## Usage

``` r
plot_DEP_barplot(res, comparisons = "all", names = NULL)
```

## Arguments

- res:

  A data frame with results from get_DEPresults()

- comparisons:

  A character vector specifying which comparisons to include. The
  default 'all' considers all comparisons present in res. The bars are
  plotted in the order of specified comparisons.

- names:

  (optional). Specifies the labels for the different bars. Defaults to
  the comparison names in colnames(res)

## Value

A ggplot bar plot

## Examples

``` r
se <- prepare_se(report.pg_matrix, expDesign)

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
res <- get_DEPresults(se, type  = 'all')
#> Tested contrasts: neg_ctrl_vs_motif1, neg_ctrl_vs_motif2, motif1_vs_motif2

plot_DEP_barplot(res)

plot_DEP_barplot(res, comparisons = c('neg_ctrl_vs_motif1',
                                       'neg_ctrl_vs_motif2'))


```
