# Rank proteins/genes based on fold change

Rank proteins/genes based on fold change

## Usage

``` r
get_ranked_genes(res)
```

## Arguments

- res:

  Results from get_DEPresults(). See details.

## Value

A ranked, named vector with FC values and gene names.

## Details

The res object should be obtained by running get_DEPresults() with a
one-vs-one (manual) comparison, otherwise it does not know which FC
column to use to rank the proteins.

## Examples

``` r
if (FALSE) { # \dontrun{
se <- prepare_se(report.pg_matrix, expDesign)
res <- get_DEPresults(se, condition1 = 'motif1', condition2 = 'neg_ctrl', type = 'manual')
GOBP <- get_genesets(db, collection = 'c5', subcollection = 'GOBP')
ranked_genes <- get_ranked_genes(res)

} # }
```
