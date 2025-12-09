# Prepare summarizedExperiment object from diann report.pg_matrix file

Prepare summarizedExperiment object from diann report.pg_matrix file

## Usage

``` r
prepare_se(
  pg_matrix,
  expDesign,
  pr_matrix = NULL,
  missing_thr = 0,
  min_peptides = 1,
  impute = "mixed",
  mixed_cutoff = "empirically",
  remove_contaminants = TRUE
)
```

## Arguments

- pg_matrix:

  the report.pg_matrix file from DIANN

- expDesign:

  A data frame with the experimental design. Should contain at least
  'label', 'condition', and 'replicate' columns.

- pr_matrix:

  Optional argument. If the report.pr_matrix file from DIANN is
  provided, peptide information will be added to output.

- missing_thr:

  Integer specifying which proteinGroups are filtered out based on
  missing values.

- min_peptides:

  An integer specifing the cutoff for razor/unique peptides. The default
  is 0.

- impute:

  Specifies which imputatation method to use (default: knn). No
  imputation is done when entering 'none'. See details for options.

- mixed_cutoff:

  Either 'empirally' or a value between 0-1. For details, see
  [mixed_imputation](https://dijkjel.github.io/DIANNmv/reference/mixed_imputation.md)

- remove_contaminants:

  A logical value specifying if potential contaminants should be removed
  from the pg_matrix.

## Value

A summarized Experiment object

## Details

For standard imputation options, see ?DEP::impute. For mixed imputation,
see
[mixed_imputation](https://dijkjel.github.io/DIANNmv/reference/mixed_imputation.md)

## Examples

``` r
se <- prepare_se(report.pg_matrix,
                expDesign, missing_thr = 1,
                impute = 'knn') # creates se with missing values imputed

#> Imputing along margin 1 (features/rows).
#> Warning: 173 rows with more than 50 % entries missing;
#>  mean imputation used for these rows
#> Cluster size 5598 broken into 2470 3128 
#> Cluster size 2470 broken into 763 1707 
#> Done cluster 763 
#> Cluster size 1707 broken into 715 992 
#> Done cluster 715 
#> Done cluster 992 
#> Done cluster 1707 
#> Done cluster 2470 
#> Cluster size 3128 broken into 2008 1120 
#> Cluster size 2008 broken into 953 1055 
#> Done cluster 953 
#> Done cluster 1055 
#> Done cluster 2008 
#> Done cluster 1120 
#> Done cluster 3128 

# creates se without imputing missing values.
se <- prepare_se(report.pg_matrix,
                expDesign,
                 missing_thr = 1,
                 impute = 'none')
```
