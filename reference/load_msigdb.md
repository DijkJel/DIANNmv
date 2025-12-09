# Downloads MSigDB

Downloads MSigDB

## Usage

``` r
load_msigdb(organism = "hs", id = "SYM")
```

## Arguments

- organism:

  Specifies for which organism to download data. 'hs' for human, 'mm'
  for mouse.

- id:

  a character, representing the ID type to use ("SYM" for gene symbols
  and "EZID" for Entrez IDs).

## Value

A GSEAbase object with genesets

## Examples

``` r
if (FALSE) { # \dontrun{
db <- load_msigdb('mm') #downloads data for mouse.
db <- load_msigdb('hs', 'EZID') #downloads data for human with EnsembleID
} # }
```
