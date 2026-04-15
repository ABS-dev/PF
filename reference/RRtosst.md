# RR exact CI, TOSST method.

Estimates confidence interval for the risk ratio or prevented fraction;
exact method based on the score statistic (inverts two one-sided tests).

## Usage

``` r
RRtosst(
  y = NULL,
  formula = NULL,
  data = NULL,
  compare = c("vac", "con"),
  alpha = 0.05,
  pf = TRUE,
  stepstart = 0.1,
  iter.max = 36,
  converge = 1e-06,
  rnd = 3,
  trace.it = FALSE,
  nuisance.points = 120,
  gamma = 1e-06
)
```

## Arguments

- y:

  Data vector c(y1, n1, y2, n2) where y are the positives, n are the
  total, and group 1 is compared to group 2 (control or reference
  group).

- formula:

  Formula of the form `cbind(y, n) ~ x`, where y is the number positive,
  n is the group size, x is a factor with two levels of treatment.

- data:

  data.frame containing variables of formula.

- compare:

  Text vector stating the factor levels: `compare[1]` is the vaccinate
  group to which `compare[2]` (control or reference) is compared.

- alpha:

  Complement of the confidence level.

- pf:

  Estimate *RR* or its complement *PF*?

- stepstart:

  starting interval for step search

- iter.max:

  Maximum number of iterations

- converge:

  Convergence criterion

- rnd:

  Number of digits for rounding. Affects display only, not estimates.

- trace.it:

  Verbose tracking of the iterations?

- nuisance.points:

  number of points over which to evaluate nuisance parameter

- gamma:

  parameter for Berger-Boos correction (restricts range of nuisance
  parameter evaluation)

## Value

A [rr1](https://abs-dev.github.io/PF/reference/rr1-class.md) object with
the following fields.

- `estimate`: vector with point and interval estimate

- `estimator`: either `"PF"` or `"RR"`

- `y`: data.frame with "y1", "n1", "y2", "n2" values.

- `rnd`: how many digits to round the display

- `alpha`: complement of confidence level

## Details

Estimates confidence intervals based on the score statistic that are
'exact' in the sense of accounting for discreteness. Inverts two
one-sided score tests. The score statistic is used to select tail area
tables, and the binomial probability is estimated over the tail area by
taking the maximum over the nuisance parameter. Algorithm is a simple
step search.

The data may also be a matrix. In that case `Y` would be entered as

`matrix(c(y1, n1-y1, y2, n2-y2), 2, 2, byrow = TRUE)`.

## References

Koopman PAR, 1984. Confidence intervals for the ratio of two binomial
proportions. *Biometrics* 40:513-517.

Agresti A, Min Y, 2001. On small-sample confidence intervals for
parameters in discrete distribution. *Biometrics* 57: 963-971.

Berger RL, Boos DD, 1994. P values maximized over a confidence set for
the nuisance parameter. *Journal of the American Statistical
Association* 89:214-220.

## See also

[RRotsst](https://abs-dev.github.io/PF/reference/RRotsst.md),
[rr1](https://abs-dev.github.io/PF/reference/rr1-class.md)

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
# Both examples represent the same observation, with data entry by vector
# and matrix notation.

y_vector <- c(4, 24, 12, 28)
RRtosst(y_vector)
#> 
#> PF 
#> 95% interval estimates
#> 
#>    PF    LL    UL 
#> 0.611 0.012 0.902 
#> 

# PF
# 95% interval estimates

#    PF    LL    UL
# 0.611 0.012 0.902

y_matrix <- matrix(c(4, 20, 12, 16), 2, 2, byrow = TRUE)
#      [, 1] [, 2]
# [1, ]    4   20
# [2, ]   12   16

RRtosst(y_matrix)
#> 
#> PF 
#> 95% interval estimates
#> 
#>    PF    LL    UL 
#> 0.611 0.012 0.902 
#> 

# PF
# 95% interval estimates

#    PF    LL    UL
# 0.611 0.012 0.902

require(dplyr)
data1 <- data.frame(group = rep(c("treated", "control"), each = 2),
  y = c(1, 3, 7, 5),
  n = c(12, 12, 14, 14),
  cage = rep(paste("cage", 1:2), 2))
data2 <- data1 |>
  group_by(group) |>
  summarize(sum_y = sum(y),
    sum_n = sum(n))
RRtosst(data = data2, formula =  cbind(sum_y, sum_n) ~ group,
  compare = c("treated", "control"))
#> 
#> PF 
#> 95% interval estimates
#> 
#>    PF    LL    UL 
#> 0.611 0.012 0.902 
#> 

# PF
# 95% interval estimates

#    PF    LL    UL
# 0.611 0.012 0.902
```
