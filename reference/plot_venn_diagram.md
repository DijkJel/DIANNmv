# Plots a Venn diagram showing overlapping significant proteins

Plots a Venn diagram showing overlapping significant proteins

## Usage

``` r
plot_venn_diagram(
  res,
  ...,
  comparisons = "all",
  colors = c("#b3e2cd", "#fdcdac", "#cbd5e8", "#f4cae4", "#e6f5c9"),
  names = NULL
)
```

## Arguments

- res:

  A data frame with results from get_DEPresults()

- ...:

  Additional results data frames (optional)

- comparisons:

  A character vector specifying which comparisons to include. The
  default 'all' considers all comparisons present in res.

- colors:

  The colors used for the Venn diagram.

- names:

  (optional). Specifies the labels for the different subsets. Defaults
  to the comparison names in colnames(res)

## Value

A Venn diagram plot

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
plot_venn_diagram(res) # Default option
#> INFO [2026-01-21 15:23:50] [[1]]
#> INFO [2026-01-21 15:23:50] venn_list
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $filename
#> INFO [2026-01-21 15:23:50] NULL
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $disable.logging
#> INFO [2026-01-21 15:23:50] T
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $fill
#> INFO [2026-01-21 15:23:50] colors
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $fontfamily
#> INFO [2026-01-21 15:23:50] [1] "sans"
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $fontface
#> INFO [2026-01-21 15:23:50] [1] "bold"
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $cat.fontfamily
#> INFO [2026-01-21 15:23:50] [1] "sans"
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $lty
#> INFO [2026-01-21 15:23:50] [1] 0
#> INFO [2026-01-21 15:23:50] 

plot_venn_diagram(res, colors = c('green', 'blue', 'red')) # Uses non-default
#> INFO [2026-01-21 15:23:50] [[1]]
#> INFO [2026-01-21 15:23:50] venn_list
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $filename
#> INFO [2026-01-21 15:23:50] NULL
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $disable.logging
#> INFO [2026-01-21 15:23:50] T
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $fill
#> INFO [2026-01-21 15:23:50] colors
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $fontfamily
#> INFO [2026-01-21 15:23:50] [1] "sans"
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $fontface
#> INFO [2026-01-21 15:23:50] [1] "bold"
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $cat.fontfamily
#> INFO [2026-01-21 15:23:50] [1] "sans"
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $lty
#> INFO [2026-01-21 15:23:50] [1] 0
#> INFO [2026-01-21 15:23:50] 

# colors
plot_venn_diagram(res, comparisons = c('neg_ctrl_vs_motif1',
                                       'neg_ctrl_vs_motif2')) # Only
#> INFO [2026-01-21 15:23:50] [[1]]
#> INFO [2026-01-21 15:23:50] venn_list
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $filename
#> INFO [2026-01-21 15:23:50] NULL
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $disable.logging
#> INFO [2026-01-21 15:23:50] T
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $fill
#> INFO [2026-01-21 15:23:50] colors
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $fontfamily
#> INFO [2026-01-21 15:23:50] [1] "sans"
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $fontface
#> INFO [2026-01-21 15:23:50] [1] "bold"
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $cat.fontfamily
#> INFO [2026-01-21 15:23:50] [1] "sans"
#> INFO [2026-01-21 15:23:50] 
#> INFO [2026-01-21 15:23:50] $lty
#> INFO [2026-01-21 15:23:50] [1] 0
#> INFO [2026-01-21 15:23:50] 

                                       # includes two comparisons

```
