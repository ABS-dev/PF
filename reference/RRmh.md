# Mantel-Haenszel method, CI for common RR over strata or clusters with sparse data.

Estimates confidence intervals for the risk ratio or prevented fraction
from clustered or stratified data, using a Mantel-Haenszel estimator for
sparse data.

## Usage

``` r
RRmh(
  formula = NULL,
  data = NULL,
  compare = c("vac", "con"),
  Y,
  alpha = 0.05,
  pf = TRUE,
  rnd = 3
)
```

## Arguments

- formula:

  Formula of the form `cbind(y, n) ~ x + cluster(w)`, where `Y` is the
  number positive, `n` is the group size, `x` is a factor with two
  levels of treatment, and `w` is a factor indicating the clusters.

- data:

  `data.frame` containing variables for formula

- compare:

  Text vector stating the factor levels: `compare[1]` is the vaccinate
  group to which `compare[2]` (control or reference) is compared.

- Y:

  Matrix of data, \\K \times 4\\. Each row is a stratum or cluster. The
  columns are \\y1, n1, y2, n2\\, where the y's are the number of
  positive in each group, and the n is the total in each group. Group 1
  corresponds to vaccinates and group 2 are controls or reference. If
  data entered by formula and dataframe, `Y` is generated automatically.

- alpha:

  Complement of the confidence level.

- pf:

  Estimate *RR* or its complement *PF*?

- rnd:

  Number of digits for rounding. Affects display only, not estimates.

## Value

An object of class
[rr1](https://abs-dev.github.io/PF/reference/rr1-class.md) with the
following fields.

- `estimate`: vector of point and interval estimates: point estimate,
  lower confidence limit, upper confidence limit

- `estimator`: either `"PF"` or `"RR"`

- `y`: data.frame of restructured input

- `rnd`: how many digits to round the display

- `alpha`: complement of confidence level

## Details

Based on the Mantel-Haenszel (1959) procedure for sparse data developed
by Greenland and Robins (1985). The confidence limits are based on
asymptotic normality of the log(risk ratio). Agresti and Hartzel (2000)
favor this procedure for small, sparse data sets, but they warn that it
is less efficient than maximum likelihood for large data sets.

## Note

If either all y1's or all y2's are zero, a division by zero may occur,
and a NaN returned for some values.

Vignette *Examples for Stratified Designs* forthcoming with more
examples.

Call to this function may be one of two formats: (1) specify `data` and
`formula` or (2) as a matrix `Y`

`RRmh(formula, data, compare = c("b", "a"), pf = TRUE, alpha = 0.05, rnd = 3)`

`RRmh(Y, pf = TRUE, alpha = 0.05, rnd = 3)`

## References

Mantel N, Haenszel W, 1959. Statistical aspects of the analysis of data
from retrospective studies of disease. *Journal of the National Cancer
Institute* 22:719-748.

Greenland S, Robins JM, 1985. Estimation of a common effect parameter
from sparse follow-up data. *Biometrics* 41: 55-68. Errata, 45:
1323-1324.

Agresti A, Hartzel J, 2000. Strategies for comparing treatments on a
binary response with multi-centre data. *Statistics in Medicine* 19:
1115-1139.

Lachin JM, 2000. *Biostatistical Methods: The Assessment of Relative
Risks* (Wiley, New York), Sec. 4.3.1.

## See also

[rr1](https://abs-dev.github.io/PF/reference/rr1-class.md)

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
## Table 1 from Gart (1985)
##  as data frame

# tx group "b" is control
RRmh(cbind(y, n) ~ tx + cluster(clus),
     Table6,
     compare = c("a", "b"), pf = FALSE)
#> 
#> RR 
#> 95% interval estimates
#> 
#>   RR   LL   UL 
#> 2.67 1.37 5.23 
#> 

# RR
# 95% interval estimates
#
#   RR   LL   UL
# 2.67 1.37 5.23
#

## or as matrix
RRmh(Y = table6, pf = FALSE)
#> 
#> RR 
#> 95% interval estimates
#> 
#>   RR   LL   UL 
#> 2.67 1.37 5.23 
#> 

# RR
# 95% interval estimates
#
#   RR   LL   UL
# 2.67 1.37 5.23
```
