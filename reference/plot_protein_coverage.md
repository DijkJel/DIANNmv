# Plots the coverage of proteins.

Plots the coverage of proteins.

## Usage

``` r
plot_protein_coverage(
  se,
  genes,
  positions = NULL,
  zoom = NULL,
  fasta = NULL,
  organism = "hs",
  combine_overlap = T,
  dodge_labels = T,
  scaling = "centered",
  condition_order = NULL
)
```

## Arguments

- se:

  The SummarizedExperiment object from
  [prepare_se](https://dijkjel.github.io/DIANNmv/reference/prepare_se.md).

- genes:

  The gene names of the proteins that you would like to plot.

- positions:

  A numeric vector indiciting the amino acid positions that you want to
  highlight.

- zoom:

  A numeric with two values (start and stop) indicating the specific
  part of the protein that you want to plot.

- fasta:

  Data frame with protein information when not using a default option.
  See
  [get_ibaq_peptides](https://dijkjel.github.io/DIANNmv/reference/get_ibaq_peptides.md)

- organism:

  Specifies which default protein database to use. 'hs' for human, or
  'mm' for mouse.

- combine_overlap:

  Boolean specifying whether to combine miscleaved peptides.

- dodge_labels:

  Boolean specifying wheter position labels should be plotted on
  different heights so that they do not overlap in the plot.

- scaling:

  Boolean value indicating whether the peptide intensities should be
  centered over the different replicates.

- condition_order:

  Character vector that specifies the order of the conditions on the
  y-axis.

## Value

A ggplot2 object showing the found peptides for a protein.

## Examples

``` r
se <- prepare_se(report.pg_matrix, expDesign, report.pr_matrix)

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
smad4 <- plot_protein_coverage(se, 'SMAD4', positions = c(100, 150))
#> Warning: Duplicated aesthetics after name standardisation: ymin
#> Coordinate system already present.
#> ℹ Adding new coordinate system, which will replace the existing one.

# Zoom in on first 200 AA of protein
smad4 <- plot_protein_coverage(se, 'SMAD4', positions = c(100, 150), zoom = c(1, 200))
#> Warning: Duplicated aesthetics after name standardisation: ymin
#> Coordinate system already present.
#> ℹ Adding new coordinate system, which will replace the existing one.
```
