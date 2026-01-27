# Changelog

## DIANNmv 1.3.2

### Bug Fix:

- The bug with overimputation in v1.3.1 is fixed.

## DIANNmv 1.3.1

### New Features:

- [`get_DEPresults()`](https://dijkjel.github.io/DIANNmv/reference/get_DEPresults.md)
  includes an extra argument `missing_thr` that will add an extra column
  to the output data frame for each 1-vs-1 comparison when set to other
  than NA. This column indicates whether too many (depending on value of
  `missing_thr`) values were imputed for this specific comparison.

- [`plotVolcano()`](https://dijkjel.github.io/DIANNmv/reference/plotVolcano.md)
  includes an extra argument `remove_overimputed`. When set to TRUE,
  only the values that are not overimputed (see above) are being
  plotted.

- Initial CRAN submission.
