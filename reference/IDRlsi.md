# IDR likelihood support interval.

Estimates likelihood support interval for the incidence density ratio or
prevented fraction based on it.

## Usage

``` r
IDRlsi(
  y = NULL,
  formula = NULL,
  data = NULL,
  alpha = 0.05,
  k = 8,
  use.alpha = FALSE,
  pf = TRUE,
  converge = 1e-08,
  rnd = 3,
  start = NULL,
  trace.it = FALSE,
  iter.max = 24,
  compare = c("con", "vac")
)
```

## Arguments

- y:

  Data vector c(y1, n1, y2, n2) where y are the positives, n are the
  total, and group 1 is compared to group 2 (control or reference).

- formula:

  Formula of the form cbind(y, n) ~ x, where y is the number positive, n
  is the group size, x is a factor with two levels of treatment.

- data:

  data.frame containing variables of the formula.

- alpha:

  Complement of the confidence level.

- k:

  Likelihood ratio criterion.

- use.alpha:

  Base choice of k on its relationship to alpha?

- pf:

  Estimate *IDR* or its complement *PF*?

- converge:

  Convergence criterion

- rnd:

  Number of digits for rounding. Affects display only, not estimates.

- start:

  describe here.

- trace.it:

  Verbose tracking of the iterations?

- iter.max:

  Maximum number of iterations

- compare:

  Text vector stating the factor levels: `compare[1]` is the vaccinate
  group to which `compare[2]` (control or reference) is compared.

## Value

A [rrsi](https://abs-dev.github.io/PF/reference/rrsi.md) object with the
following elements.

- `estimate`: vector with point and interval estimate

- `estimator`: either *PF* or *IDR*

- `y`: data.frame with "y1", "n1", "y2", "n2" values.

- `k`: Likelihood ratio criterion

- `rnd`: how many digits to round the display

- `alpha`: complement of confidence level

## Details

Estimates likelihood support interval for the incidence density ratio
based on orthogonal factoring of reparameterized likelihood. The
incidence density is the number of cases per subject-time; its
distribution is assumed Poisson.

Likelihood support intervals are usually formed based on the desired
likelihood ratio, often 1 / 8 or 1 / 32. Under some conditions the log
likelihood ratio may follow the chi square distribution. If so, then
\\\alpha = 1 - F(2log(k), 1)\\, where \\F\\ is a chi-square CDF. if
``` use.alpha = TRUE``RRsc() ``` will make the conversion from
\\\alpha\\ to \\k\\ .

The data may also be a matrix, in which case `y` would be entered as
`matrix(c(y1, n1 - y1, y2, n2 - y2), 2, 2, byrow = TRUE)`.

## References

Royall R. *Statistical Evidence: A Likelihood Paradigm*. Chapman & Hall,
Boca Raton, 1997. Section 7.2.

## See also

[IDRsc](https://abs-dev.github.io/PF/reference/IDRsc.md)

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
# Both examples represent the same observation, with data entry by vector
# and matrix notation.

y_vector <- c(26, 204, 10, 205)
IDRlsi(y_vector, pf = FALSE)
#> 
#> 1/8 likelihood support interval for IDR 
#> 
#> corresponds to 95.858% confidence
#>   (under certain assumptions)
#> 
#> IDR 
#>  IDR   LL   UL 
#> 2.61 1.26 5.88 
#> 

# 1 / 8 likelihood support interval for IDR

# corresponds to 95.858% confidence
#   (under certain assumptions)

# IDR
#  IDR   LL   UL
# 2.61 1.26 5.88

y_matrix <- matrix(c(26, 178, 10, 195), 2, 2, byrow = TRUE)
y_matrix
#>      [,1] [,2]
#> [1,]   26  178
#> [2,]   10  195
#      [, 1] [, 2]
# [1, ]   26  178
# [2, ]   10  195

IDRlsi(y_matrix, pf = FALSE)
#> 
#> 1/8 likelihood support interval for IDR 
#> 
#> corresponds to 95.858% confidence
#>   (under certain assumptions)
#> 
#> IDR 
#>  IDR   LL   UL 
#> 2.61 1.26 5.88 
#> 

# 1 / 8 likelihood support interval for IDR

# corresponds to 95.858% confidence
#   (under certain assumptions)

# IDR
#  IDR   LL   UL
# 2.61 1.26 5.88

data1 <- data.frame(group = rep(c("treated", "control"), each = 5),
             n = c(rep(41, 4), 40, rep(41, 5)),
             y = c(4, 5, 7, 6, 4, 1, 3, 3, 2, 1),
             cage = rep(paste("cage", 1:5), 2))
IDRlsi(data = data1, formula = cbind(y, n) ~ group,
               compare = c("treated", "control"), pf = FALSE)
#> 'y' and 'n' values will be summed by treated and control designation in datacolumn group
#> 
#> 1/8 likelihood support interval for IDR 
#> 
#> corresponds to 95.858% confidence
#>   (under certain assumptions)
#> 
#> IDR 
#>  IDR   LL   UL 
#> 2.61 1.26 5.88 
#> 

# 1 / 8 likelihood support interval for IDR

# corresponds to 95.858% confidence
#   (under certain assumptions)

# IDR
#  IDR   LL   UL
# 2.61 1.26 5.88

require(dplyr)
#> Loading required package: dplyr
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
data2 <- data1 |>
  group_by(group) |>
  summarize(sum_y = sum(y),
    sum_n = sum(n))

IDRlsi(data = data2, formula = cbind(sum_y, sum_n) ~ group,
               compare = c("treated", "control"), pf = FALSE)
#> 
#> 1/8 likelihood support interval for IDR 
#> 
#> corresponds to 95.858% confidence
#>   (under certain assumptions)
#> 
#> IDR 
#>  IDR   LL   UL 
#> 2.61 1.26 5.88 
#> 

# 1 / 8 likelihood support interval for IDR

# corresponds to 95.858% confidence
#   (under certain assumptions)

# IDR
#  IDR   LL   UL
# 2.61 1.26 5.88
```
