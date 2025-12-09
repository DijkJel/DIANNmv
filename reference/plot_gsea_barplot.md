# Plots GSEA results as bar plot.

Plots GSEA results as bar plot.

## Usage

``` r
plot_gsea_barplot(
  gsea,
  pos_color = "gold1",
  neg_color = "darkblue",
  padj_cutoff = 0.05,
  top_n = Inf,
  remove_prefix = F,
  max_name_length = Inf,
  break_names = T
)
```

## Arguments

- gsea:

  Output from
  [perform_GSEA](https://dijkjel.github.io/DIANNmv/reference/perform_GSEA.md).

- pos_color:

  Character vector specifying which color to use for activated pathways.

- neg_color:

  Character vector specifying which color to use for repressed pathways.

- padj_cutoff:

  The maximum p.adjust value allowed for inclusion of the pathway.

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

A ggplot bar plot object

## Examples

``` r
if (FALSE) { # \dontrun{
gsea <- perform_GSEA(res, genesets)

# Bar plot of 10 pathways with lowest padj values.
barplot <- plot_gsea_barplot(gsea, top_n = 10)

#Bar plot of 10 pathways with lowest values and different colors for bars.
barplot <- plot_gsea_barplot(gsea, pos_color = 'red3', neg_color = 'dodgerblue', top_n = 10)
} # }
```
