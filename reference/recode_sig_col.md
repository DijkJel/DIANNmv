# Relabel the significance column in the DEP results data frame

Relabel the significance column in the DEP results data frame

## Usage

``` r
recode_sig_col(res, pval_cutoff = 0.05, fc_cutoff = 1)
```

## Arguments

- res:

  Results data frame after DEP analysis.

- pval_cutoff:

  The p.adj significance cutoff

- fc_cutoff:

  The absolute log2 fold change significance cutoff

## Value

The results data frame with adjusted significance column

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
res <- recode_sig_col(res, pval_cutoff = 0.01, fc_cutoff = 1)
```
