# Gart-Nam method, CI for common RR over strata or clusters.

Estimates confidence intervals for the risk ratio or prevented fraction
from clustered or stratified data.

## Usage

``` r
RRstr(
  formula = NULL,
  data = NULL,
  compare = c("vac", "con"),
  Y,
  alpha = 0.05,
  pf = TRUE,
  trace.it = FALSE,
  iter.max = 24,
  converge = 1e-06,
  rnd = 3,
  multiplier = 0.7,
  divider = 1.1
)
```

## Arguments

- formula:

  Formula of the form `cbind(y, n) ~ x + cluster(w)`, where y is the
  number positive, n is the group size, x is a factor with two levels of
  treatment, and w is a factor indicating the clusters.

- data:

  data.frame containing variables of formula

- compare:

  Text vector stating the factor levels: `compare[1]` is the control or
  reference group to which `compare[2]` is compared

- Y:

  Matrix of data. Each row is a stratum or cluster. The columns are y2,
  n2, y1, n1. If data entered by formula and dataframe, Y is generated
  automatically.

- alpha:

  Size of the homogeneity test and complement of the confidence level.

- pf:

  Estimate *RR* or its complement *PF*?

- trace.it:

  verbose tracking of the iterations?

- iter.max:

  Maximum number of iterations

- converge:

  Convergence criterion

- rnd:

  Number of digits for rounding. Affects display only, not estimates.

- multiplier:

  internal control parameter for algorithm

- divider:

  internal control parameter for algorithm

## Value

A [rrstr](https://abs-dev.github.io/PF/reference/rrstrclass.md) object
with the following fields:

- `estimate`: matrix of point and interval estimates - starting value,
  MLE, and skewness corrected

- `hom`: list of homogeneity statistic, p-value, and degrees of freedom,
  or error message if appropriate.

- `estimator`: either `"PF"` or `"RR"`

- `y`: `data.frame` of restructured input

- `compare`: groups compared

- `rnd`: how many digits to round the display

- `alpha`: size of test; complement of confidence level

## Details

Uses the DUD algorithm to estimate confidence intervals by the method of
Gart.

## Note

Vignette *Examples for Stratified Designs* forthcoming with more
examples.

Call to this function may be one of two formats: (1) specify `data` and
`formula` or (2) as a matrix `Y`

`RRstr(formula, data, compare = c("b", "a"), pf = TRUE, alpha = 0.05, trace.it = FALSE, iter.max = 24, converge = 1e-6, rnd = 3, multiplier = 0.7, divider = 1.1)`

`RRstr(Y, compare = c("b", "a"), pf = TRUE, alpha = 0.05, trace.it = FALSE, iter.max = 24, converge = 1e-6, rnd = 3, multiplier = 0.7, divider = 1.1)`

## References

Gart JJ, 1985. Approximate tests and interval estimation of the common
relative risk in the combination of \\2 x 2\\ tables. *Biometrika*
72:673-677.

Gart JJ, Nam J, 1988. Approximate interval estimation of the ratio of
binomial parameters: a review and corrections for skewness. *Biometrics*
44:323-338.

Ralston ML, Jennrich RI, 1978. DUD, A Derivative-Free Algorithm for
Nonlinear Least Squares. *Technometrics* 20:7-14.

## See also

[rrstr](https://abs-dev.github.io/PF/reference/rrstrclass.md)

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
## Table 1 from Gart (1985)
##  as data frame
## "b" is control group
RRstr(cbind(y, n) ~ tx + cluster(clus),
      Table6,
      compare = c("a", "b"), pf = FALSE)
#> 
#> Test of homogeneity across clusters
#> 
#> stat  0.954 
#> df    3 
#> p     0.812 
#> 
#> RR 
#> 95% interval estimates
#> 
#>             RR   LL   UL
#> starting  2.66 1.37 5.18
#> mle       2.65 1.39 5.03
#> skew corr 2.65 1.31 5.08
#> 

# Test of homogeneity across clusters

# stat     0.954
# df       3
# p        0.812

#  RR estimates

#             RR   LL   UL
# starting   2.66 1.37 5.18
# mle        2.65 1.39 5.03
# skew corr  2.65 1.31 5.08

## or as matrix
RRstr(Y = table6, pf = FALSE)
#> 
#> Test of homogeneity across clusters
#> 
#> stat  0.954 
#> df    3 
#> p     0.812 
#> 
#> RR 
#> 95% interval estimates
#> 
#>             RR   LL   UL
#> starting  2.66 1.37 5.18
#> mle       2.65 1.39 5.03
#> skew corr 2.65 1.31 5.08
#> 

tst <- data.frame(y = c(0, 2, 0, 4, 0, 3, 0, 7),
n = rep(10, 8),
tx = rep(c("a", "b"), 4),
clus = rep(paste("Row", 1:4, sep = ""), each = 2))
```
