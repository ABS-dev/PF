# Rao-Scott weights.

Rao-Scott weights.

## Usage

``` r
rsb(y = NULL, n = NULL, formula = NULL, data = NULL, id = NULL)
```

## Arguments

- y:

  vector of number positive.

- n:

  vector of total number.

- formula:

  Formula of the form cbind(y, n) ~ id, where y is the number positive,
  n is the total number, id is a factor for estimating the weights by
  subset.

- data:

  data.frame containing variables of formula.

- id:

  vector of factor for estimating the weights by subset.

## Value

A list with the following elements.

- `w`: vector of weights

- `d`: vector of \\{d}\_{i}\\ estimates

## Details

Estimates the cluster design effect \\{d}\_{i}\\ as the variance
inflation due to clustering by the method of Rao and Scott. `rsb`
estimates the \\{d}\_{i}\\ for use by `rsbWt` or other functions.

## References

Rao JNK, Scott AJ, 1992. A simple method for the analysis of clustered
binary data. *Biometrics* 48:577-585.

## See also

[rsbWt](https://abs-dev.github.io/PF/reference/rsbWt.md).

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
# Weil's rat data (Table 1 of Rao and Scott)
rsb(rat$y, rat$n, id = rat$group)$d
#>  control  treated 
#> 1.232495 3.952861 
#  control  treated
# 1.232495 3.952861
rsb(data = rat, formula = cbind(y, n) ~ group)$d
#>  control  treated 
#> 1.232495 3.952861 
#  control  treated
# 1.232495 3.952861
```
