# Prepare data for volcano plots

Prepare data for volcano plots

## Usage

``` r
prepare_volcano_data(
  res,
  pval_cutoff = 0.05,
  fc_cutoff = 1,
  label = "sig",
  top_n = NULL,
  up_color = "red3",
  down_color = "dodgerblue",
  ns_color = "grey70",
  remove_overimputed = F
)
```

## Arguments

- res:

  Data frame with results from get_DEPresults.

- pval_cutoff:

  The p.adj significance cutoff.

- fc_cutoff:

  The log2 fold change significance cutoff.

- label:

  Specifies which points to label. The default is 'sig', labeling all
  significant points. Entering a value for top_n limits the labeling to
  the top_n up- and down-regulated proteins based on the p.adj. When
  providing a vector with protein names, onlythose points are labeled.

- top_n:

  Specifies how many significant points to label.

- up_color:

  Specifies the color of the significantly up-regulated proteins.

- down_color:

  Specifies the color of the significantly down-regulated proteins.

- ns_color:

  Specifies the color of non-significant proteins.

- remove_overimputed:

  Boolean value. If set to TRUE, proteins with too many imputed values
  in 1-vs-1 comparisons are removed. Requires 'missing_thr' set in
  [get_DEPresults](https://dijkjel.github.io/DIANNmv/reference/get_DEPresults.md).

## Value

A list with a data frame for each comparison present in res.

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
res <- get_DEPresults(se, type = 'all')
#> Tested contrasts: neg_ctrl_vs_motif1, neg_ctrl_vs_motif2, motif1_vs_motif2
data <- prepare_volcano_data(res)
```
