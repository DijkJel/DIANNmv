# Perform in silico tryptic digest on uniprot fasta file

Perform in silico tryptic digest on uniprot fasta file

## Usage

``` r
get_ibaq_peptides(fasta_location, miscleaves = 0, enzyme = "trypsin")
```

## Arguments

- fasta_location:

  Path to fasta file

- miscleaves:

  Allowed number of miscleavages. Defaults to 0 (standard for iBAQ)

- enzyme:

  The enzyme used for digestion. Defaults to 'trypsin'. See 'cleaver'
  package vignette for more options.

## Value

A data frame with the number of theoretically observable peptides for
each entry in the fasta file.

## Examples

``` r
if (FALSE) { # \dontrun{
no_ibaq_peptides <- get_ibaq_peptides('path/to/fasta.fasta')
} # }
```
