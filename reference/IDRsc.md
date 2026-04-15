# IDR confidence interval.

Estimates confidence interval for the incidence density ratio or
prevented fraction based on it.

## Usage

``` r
IDRsc(
  y = NULL,
  data = NULL,
  formula = NULL,
  compare = c("con", "vac"),
  alpha = 0.05,
  pf = TRUE,
  rnd = 3
)
```

## Arguments

- y:

  Data vector c(y1, n1, y2, n2) where y are the positives, n are the
  total, and group 1 is compared to group 2 (control or reference).

- data:

  data.frame containing variables of formula.

- formula:

  Formula of the form cbind(y, n) ~ x, where y is the number positive, n
  is the group size, x is a factor with two levels of treatment.

- compare:

  Text vector stating the factor levels: `compare[1]` is the vaccinate
  group to which `compare[2]` (control or reference) is compared.

- alpha:

  Complement of the confidence level.

- pf:

  Estimate *IDR*, or its complement *PF*?

- rnd:

  Number of digits for rounding. Affects display only, not estimates.

## Value

A [rr1](https://abs-dev.github.io/PF/reference/rr1-class.md) object with
the following elements.

- `estimate`: vector with point and interval estimate

- `estimator`: either *PF* or *IDR*

- `y`: data vector

- `rnd`: how many digits to round the display

- `alpha`: complement of confidence level

## Details

The incidence density is the number of cases per subject-time; its
distribution is assumed Poisson. `IDRsc` estimates a confidence interval
for the incidence density ratio using Siev's formula based on the
Poisson score statistic. \\ IDR = \widehat{IDR}\left\\ 1 + \left(
\frac{1}{{{y}\_{1}}} + \frac{1}{{{y}\_{2}}} \right)\frac{z\_{\alpha /
2}^{2}}{2}\\ \\ \pm \\ \\ \frac{z\_{\alpha /
2}^{2}}{2{{y}\_{1}}{{y}\_{2}}}\sqrt{{{y}\_{\bullet }}\left(
{{y}\_{\bullet }}z\_{\alpha / 2}^{2} + 4{{y}\_{1}}{{y}\_{2}} \right)}
\right\\ \\

The data may also be a matrix. In that case `y` would be entered as

`matrix(c(y1, n1 - y1, y2, n2 - y2), 2, 2, byrow = TRUE)`.

## References

Siev D, 1994. Estimating vaccine efficacy in prospective studies.
*Preventive Veterinary Medicine* 20:279-296, Appendix 1.

Graham PL, Mengersen K, Morton AP, 2003. Confidence limits for the ratio
of two rates based on likelihood scores:non-iterative method *Statistics
in Medicine* 22:2071-2083.

Siev D, 2004. Letter to the editor. *Statistics in Medicine* 23:693.
(Typographical error in formula: replace the two final minus signs with
subscript dots.)

## See also

[IDRlsi](https://abs-dev.github.io/PF/reference/IDRlsi.md)

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
# All examples represent the same observation, with data entry by vector,
# matrix, and formula+data notation.

y_vector <- c(26, 204, 10, 205)
IDRsc(y_vector, pf = FALSE)
#> 
#> IDR 
#> 95% interval estimates
#> 
#>  IDR   LL   UL 
#> 2.61 1.28 5.34 
#> 

# IDR
# 95% interval estimates

#  IDR   LL   UL
# 2.61 1.28 5.34

y_matrix <- matrix(c(26, 178, 10, 195), 2, 2, byrow = TRUE)
y_matrix
#>      [,1] [,2]
#> [1,]   26  178
#> [2,]   10  195
#      [, 1] [, 2]
# [1, ]   26  178
# [2, ]   10  195

IDRsc(y_matrix, pf = FALSE)
#> 
#> IDR 
#> 95% interval estimates
#> 
#>  IDR   LL   UL 
#> 2.61 1.28 5.34 
#> 

# IDR
# 95% interval estimates

#  IDR   LL   UL
# 2.61 1.28 5.34

require(dplyr)
data1 <- data.frame(group = rep(c("treated", "control"), each = 5),
            n = c(rep(41, 4), 40, rep(41, 5)),
            y = c(4, 5, 7, 6, 4, 1, 3, 3, 2, 1),
            cage = rep(paste("cage", 1:5), 2))
data2 <- data1 |>
  group_by(group) |>
  summarize(sum_y = sum(y),
  sum_n = sum(n))
IDRsc(data = data2, formula =  cbind(sum_y, sum_n) ~ group,
    compare = c("treated", "control"), pf = FALSE)
#> 
#> IDR 
#> 95% interval estimates
#> 
#>  IDR   LL   UL 
#> 2.61 1.28 5.34 
#> 

# IDR
# 95% interval estimates

#  IDR   LL   UL
# 2.61 1.28 5.34
```
