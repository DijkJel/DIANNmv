# Prepares data for the different GSEA visualization options

Prepares data for the different GSEA visualization options

## Usage

``` r
prepare_gsea_data(
  gsea,
  padj_cutoff = 0.05,
  top_n = Inf,
  remove_prefix = F,
  break_names = T,
  max_name_length = Inf
)
```

## Arguments

- gsea:

  Output from
  [perform_GSEA](https://dijkjel.github.io/DIANNmv/reference/perform_GSEA.md).

- padj_cutoff:

  The maximum p.adjust value allowed for inclusion of the pathway.

- top_n:

  The maximum number of pathways included. Takes the top_n pathways with
  the lowest p.adj values.

- remove_prefix:

  Boolean specifying to remove the prefix from pathway names.

- break_names:

  Boolean value that indicates if pathway names should be shown on two
  lines when the length exceedds the max_name_length parameter

- max_name_length:

  Numeric value indicating the maximum characters of the pathway names
  before they are shortened.

## Value

A data frame with filtered GSEA data.

## Examples
