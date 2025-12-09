# Creates individual GSEA enrichment plots for selected pathwways.

Creates individual GSEA enrichment plots for selected pathwways.

## Usage

``` r
plot_gsea_enrichment_plot(
  genesets,
  pathways,
  res,
  line_color = "green",
  ES_color = "red",
  linewidth = 1,
  ticksSize = 0.2
)
```

## Arguments

- genesets:

  A list with pathways/genesets created with
  [get_genesets](https://dijkjel.github.io/DIANNmv/reference/get_genesets.md).

- pathways:

  A character vector specifying which pathways need to be plotted.

- res:

  Results from
  [get_DEPresults](https://dijkjel.github.io/DIANNmv/reference/get_DEPresults.md).
  Requires a 1-vs-1 comparison.

- line_color:

  Character specifying which color the solid line in the plot is.

- ES_color:

  Character specifying which colro the dashed line specifying ES-score
  in the plot is.

- linewidth:

  Numeric value specifying the thickness of the plotted lines.

- ticksSize:

  Numeric value specifying the thickness of the ticks in lower part
  graph.

## Value

An single or list of ggplot2 objects.

## Examples

``` r
if (FALSE) { # \dontrun{
db <- load_msigdb()
genesets <- get_genesets(db, 'h')
res <- get_DEPresults(se, 'motif1', 'neg_ctrl')
pathways <- c('HALLMARK_MYC_TARGETS_V1', 'HALLMARK_TGF_BETA_SIGNALING')
enrichment_plots <- plot_gsea_enrichment_plot(genesets, pathways, res)
} # }
```
