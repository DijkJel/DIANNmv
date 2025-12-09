# Adds a Potential.contaminant column to pg_matrix based on MaxQuant contaminants.txt

Adds a Potential.contaminant column to pg_matrix based on MaxQuant
contaminants.txt

## Usage

``` r
add_contaminants(pg_matrix)
```

## Arguments

- pg_matrix:

  The report.pg_matrix output file from DIANN

## Value

A data frame with a added Potential.contaminant column.

## Examples

``` r
pg_matrix <- add_contaminants(report.pg_matrix)
```
