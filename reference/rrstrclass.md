# Data class rrstr

data class rrstr

## Fields

- `estimate`: vector with point and interval estimate

- `rnd`: how many digits to round display

- `alpha`: complement of c.i.

- `estimator`: either `"PF"` or `"RR"`

- `hom`: list of homogeneity statistic, p-value, and degrees of freedom.
  If `Phi == 0 | Phi == 1`, homogeneity test is not possible and error
  message displays

- `Y`: data.frame of restructured input

- `compare`: groups compared

## See also

rrstr

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)
