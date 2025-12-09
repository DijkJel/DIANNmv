# Perform GSEA based on DEP results

Perform GSEA based on DEP results

## Usage

``` r
perform_GSEA(res, genesets)
```

## Arguments

- res:

  Results from get_DEPresults(). See details.

- genesets:

  MSigDB gene sets retrieved with get_genesets()

## Value

A data frame with GSEA results

## Details

The res object should be obtained by running get_DEPresults() with a
one-vs-one (manual) comparison, otherwise it does not know which FC
column to use to rank the proteins.

## Examples

``` r
if (FALSE) { # \dontrun{
se <- prepare_se(report.pg_matrix, expDesign)
res <- get_DEPresults(se, condition_1 = 'motif1', condition2 = 'neg_ctrl', type = 'manual')
GOBP <- get_genesets(db, 'c5', 'GOBP')
gsea <- perform_GSEA(res, GOBP)
} # }
```
