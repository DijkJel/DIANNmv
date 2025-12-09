# Extract a collection of gene sets from the MSigDB object

Extract a collection of gene sets from the MSigDB object

## Usage

``` r
get_genesets(
  db,
  collection = c("c2", "c5", "h"),
  subcollection = c("GOBP", "GOCC", "GOMF", "KEGG")
)
```

## Arguments

- db:

  An MSigDB object downloaded with load_msigdb()

- collection:

  Specifies which collection of gene sets to use. Options are: c2 (e.g.
  Reactome, KEGG), c5 (GO sets), or h (cancer hallmarks)

- subcollection:

  Specifies which subcollection to use. E.g. 'Reactome' or 'KEGG' for
  c2, or 'GOBP/GOCC/GOMF' for 'c5'

## Value

A list with gene sets

## Examples

``` r
if (FALSE) { # \dontrun{
Reactome <- get_genesets(db, collection = 'c2', subcollection = 'Reactome')
} # }
```
