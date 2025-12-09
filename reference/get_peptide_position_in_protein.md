# Find the amino acid locations of peptides within a protein sequence

Find the amino acid locations of peptides within a protein sequence

## Usage

``` r
get_peptide_position_in_protein(peptide_sequences, full_sequence)
```

## Arguments

- peptide_sequences:

  A character vector with the peptide sequences.

- full_sequence:

  A character vector with the sequence of the parent protein.

## Value

A data frame with the peptide sequences, their start position and end
position in the parent protein.

## Examples

``` r
pep_sequence <- 'AARTGGGPLR'
full_sequence <- 'ISIFQKRSRWCAARTGGGPLRDPEIYLWWW'
positions <- get_peptide_position_in_protein(pep_sequence, full_sequence)
```
