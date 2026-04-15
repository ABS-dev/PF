# RR likelihood support interval.

likelihood support interval for the risk ratio or prevented fraction by
the likelihood profile.

## Usage

``` r
RRlsi(
  y = NULL,
  formula = NULL,
  data = NULL,
  compare = c("vac", "con"),
  alpha = 0.05,
  k = 8,
  use.alpha = FALSE,
  pf = TRUE,
  iter.max = 50,
  converge = 1e-06,
  rnd = 3,
  start = NULL,
  track = FALSE,
  full.track = FALSE
)
```

## Arguments

- y:

  Data vector `c(y1, n1, y2, n2)` where `y` are the positives, `n` are
  the total, and group 1 is compared to group 2 (control or reference
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

  Complement of the confidence level (see details).

- k:

  Likelihood ratio criterion.

- use.alpha:

  Base choice of k on its relationship to alpha?

- pf:

  Estimate *RR* or its complement *PF*?

- iter.max:

  Maximum number of iterations

- converge:

  Convergence criterion

- rnd:

  Number of digits for rounding. Affects display onlyRR, not estimates.

- start:

  Optional starting value.

- track:

  Verbose tracking of the iterations?

- full.track:

  Verbose tracking of the iterations?

## Value

An object of class
[rrsi](https://abs-dev.github.io/PF/reference/rrsi.md) with the
following fields: `estimate`: matrix of point and interval estimates -
see details `estimator`: either `"PF"` or `"RR"` `y`: data.frame with
"y1", "n1", "y2", "n2" values. `rnd`: how many digits to round the
display `k`: likelihood ratio criterion `alpha`: complement of
confidence level

## Details

Estimates a likelihood support interval for *RR* or *PF* by the profile
likelihood method using the DUD algorithm.

Likelihood support intervals are usually formed based on the desired
likelihood ratio, often 1 / 8 or 1 / 32. Under some conditions the log
likelihood ratio may follow the chi square distribution. If so, then
\\\alpha = 1 - F(2log(k), 1)\\, where \\F\\ is a chi-square CDF. if
`use.alpha = TRUE`, `RRlsi()` will make the conversion from \\\alpha\\
to \\k\\

The data may also be a matrix. In that case `Y` would be entered as

`matrix(c(y1, n1-y1, y2, n2-y2), 2, 2, byrow = TRUE)`.

## References

Royall R. *Statistical Evidence: A Likelihood Paradigm*. Chapman & Hall,
Boca Raton, 1997. Section 7.6

Ralston ML, Jennrich RI, 1978. DUD, A Derivative-Free Algorithm for
Nonlinear Least Squares. *Technometrics* 20:7-14.

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
# All examples represent the same observation, with data entry by vector,
# matrix, and formula+data notation.

y_vector <- c(4, 24, 12, 28)
RRlsi(y_vector)
#> 
#> 1/8 likelihood support interval for PF 
#> 
#> corresponds to 95.858% confidence
#>   (under certain assumptions)
#> 
#> PF 
#>     PF     LL     UL 
#> 0.6111 0.0168 0.8859 
#> 

# 1 / 8 likelihood support interval for PF

# corresponds to 95.858% confidence
#   (under certain assumptions)

# PF
#     PF     LL     UL
# 0.6111 0.0168 0.8859

y_matrix <- matrix(c(4, 20, 12, 16), 2, 2, byrow = TRUE)
y_matrix
#>      [,1] [,2]
#> [1,]    4   20
#> [2,]   12   16
#      [, 1] [, 2]
# [1, ]    4   20
# [2, ]   12   16

RRlsi(y_matrix)
#> 
#> 1/8 likelihood support interval for PF 
#> 
#> corresponds to 95.858% confidence
#>   (under certain assumptions)
#> 
#> PF 
#>     PF     LL     UL 
#> 0.6111 0.0168 0.8859 
#> 

# 1 / 8 likelihood support interval for PF

# corresponds to 95.858% confidence
#   (under certain assumptions)

# PF
#     PF     LL     UL
# 0.6111 0.0168 0.8859

require(dplyr)
data1 <- data.frame(group = rep(c("treated", "control"), each = 2),
  y = c(1, 3, 7, 5),
  n = c(12, 12, 14, 14),
  cage = rep(paste("cage", 1:2), 2))

data2 <- data1 |>
  group_by(group) |>
  summarize(sum_y = sum(y),
    sum_n = sum(n))
RRlsi(data = data2, formula =  cbind(sum_y, sum_n) ~ group,
   compare = c("treated", "control"))
#> 
#> 1/8 likelihood support interval for PF 
#> 
#> corresponds to 95.858% confidence
#>   (under certain assumptions)
#> 
#> PF 
#>     PF     LL     UL 
#> 0.6111 0.0168 0.8859 
#> 

# 1 / 8 likelihood support interval for PF
#
# corresponds to 95.858% confidence
# (under certain assumptions)
#
# PF
# PF     LL     UL
# 0.6111 0.0168 0.8859
```
