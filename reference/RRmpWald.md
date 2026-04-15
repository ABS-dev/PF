# Wald confidence intervals for RR from matched pairs

Estimates confidence intervals for the risk ratio or prevented fraction
from matched pairs.

## Usage

``` r
RRmpWald(
  formula = NULL,
  data = NULL,
  compare = c("vac", "con"),
  affected = 1,
  x,
  alpha = 0.05,
  pf = TRUE,
  tdist = TRUE,
  df = NULL,
  rnd = 3
)
```

## Arguments

- formula:

  Formula of the form `y ~ x + cluster(w)`, where `y` is the indicator
  for an individual's positive response, `x` is a factor with two levels
  of treatment, and `w` identifies the pairs.

- data:

  `data.frame` containing variables in formula

- compare:

  Text vector stating the factor levels: `compare[1]` is the vaccinate
  group to which `compare[2]` (control or reference) is compared.

- affected:

  Indicator for positive response

- x:

  Alternative data input. Instead of formula and data frame, data may be
  input as frequency vector. See example for how to order this vector.

- alpha:

  Complement of the confidence level

- pf:

  Estimate *RR* or its complement *PF*?

- tdist:

  Use t distribution?

- df:

  Degrees of freedom. When NULL, the function will default to \`df = N

  - 2\`, where N is the total number of pairs.

- rnd:

  Number of digits for rounding. Affects display only, not estimates.

## Value

A [rrmp](https://abs-dev.github.io/PF/reference/rrmp.md) object with the
following fields:

- `estimate`: vector of point and interval estimates - see details

- `estimator`: either `"PF"` or `"RR"`

- `compare`: text vector, same as input

- `alpha`: complement of confidence level

- `rnd`: how many digits to round the display

- `multvec`: data frame showing the multinomial representation of the
  data

## Details

Estimates confidence intervals for the risk ratio or prevented fraction
from matched pairs. The response is the tetranomial vector
`c(11, 12, 21, 22)`, where the first index is the row and the the second
index is the column when displayed as a 2x2 table. Wald type confidence
intervals are found by applying the delta method to the multinomial
variance. This method fails when there are no responders in one of the
treatment groups.

Alternative forms of data entry are illustrated by the output, say `Y`,
where `c(Y$xtable) = Y$freqvec = Y$multvec$Freq`.

If `RR = 0` (`PF = 1`), the function will return degenerate interval.

## Note

Experimental functions for estimating profile likelihood intervals are
in the CVBmisc package.

Call to this function may be one of two formats: (1) specify `data` and
`formula` or (2) as a vector `x`

`RRmpWald(formula, data, compare = c("vac", "con"), affected = 1, alpha = 0.05, pf = TRUE, tdist = TRUE, df = NULL, rnd = 3)`

`RRmpWald(x, compare = c("vac", "con"), affected = 1, alpha = 0, 05, pf = TRUE, tdist = TRUE, df = NULL, rnd = 3)`

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
RRmpWald(pos ~ tx + cluster(cage), New, compare = c("vac", "con"))
#> 
#> PF 
#> 95% interval estimates
#> 
#>    PF    LL    UL 
#> 0.550 0.183 0.752 
#> 

thistable <- New |>
  tidyr::spread(tx, pos) |>
  tidyr::drop_na() |>
  dplyr::mutate(vac = factor(vac, levels = 1:0),
    con = factor(con, levels = 1:0)) |>
  with(table(vac, con))
thistable
#>    con
#> vac  1  0
#>   1  7  2
#>   0 13  4
as.vector(thistable)
#> [1]  7 13  2  4

RRmpWald(x = as.vector(thistable))
#> 
#> PF 
#> 95% interval estimates
#> 
#>    PF    LL    UL 
#> 0.550 0.183 0.752 
#> 
```
