# Plots GSEA data as volcano plot

Plots GSEA data as volcano plot

## Usage

``` r
plot_gsea_volcano(
  gsea,
  padj_cutoff = 0.05,
  label = "sig",
  top_n = NULL,
  up_color = "red3",
  down_color = "dodgerblue",
  ns_color = "grey70"
)
```

## Arguments

- gsea:

  gsea Output from
  [perform_GSEA](https://dijkjel.github.io/DIANNmv/reference/perform_GSEA.md).

- padj_cutoff:

  padj_cutoff The maximum p.adjust value allowed for inclusion of the
  pathway.

- label:

  Specifies which points to label. The default is 'sig', labeling all
  significant points. Entering a value for top_n limits the labeling to
  the top_n up- and top_n down-regulated proteins based on the p.adj.
  When providing a vector with protein names, only those points are
  labeled.

- top_n:

  Specifies how many significant points to label.

- up_color:

  Character string specifying the color for significant points with
  positive NES.

- down_color:

  Character string specifying the color for significant points with
  positive NES.

- ns_color:

  Character string specifying the color for non-significant points.

## Value

A ggplot2 object

## Examples

``` r
if (FALSE) { # \dontrun{
gsea <- perform_GSEA(res, genesets)
volcano <- plot_gsea_volcano(gsea, top_n = 10)
} # }
```
