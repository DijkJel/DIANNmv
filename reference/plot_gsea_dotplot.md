# Plots GSEA data as a faceted dot plot showing NES, padj, and set size.

Plots GSEA data as a faceted dot plot showing NES, padj, and set size.

## Usage

``` r
plot_gsea_dotplot(
  gsea,
  padj_cutoff = 0.05,
  top_n = Inf,
  remove_prefix = F,
  max_name_length = Inf,
  break_names = T
)
```

## Arguments

- gsea:

  gsea Output from
  [perform_GSEA](https://dijkjel.github.io/DIANNmv/reference/perform_GSEA.md).

- padj_cutoff:

  padj_cutoff The maximum p.adjust value allowed for inclusion of the
  pathway.

- top_n:

  The maximum number of pathways included. Takes the top_n pathways with
  the lowest p.adj values.

- remove_prefix:

  Boolean specifying to remove the prefix from pathway names.

- max_name_length:

  Numeric value specifying the max length of pathway names.

- break_names:

  Boolean value that indicates if pathway names should be shown on two
  lines when the length exceedds the max_name_length parameter

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
gsea <- perform_GSEA(res, genesets)
dotplot <- plot_gsea_dotplot(gsea, top_n = 10)
} # }
```
