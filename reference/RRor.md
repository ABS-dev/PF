# RR estimate from logistic regression.

Model based interval estimate of the risk ratio or prevented fraction
from a logistic regression model.

## Usage

``` r
RRor(
  fit = NULL,
  beta.hat = NULL,
  var.beta.hat = NULL,
  degf = NULL,
  which = c(1, 2),
  pf = TRUE,
  norm = FALSE,
  alpha = 0.05,
  rnd = 3
)
```

## Arguments

- fit:

  A [glm](https://rdrr.io/r/stats/glm.html) object.

- beta.hat:

  Parameters estimates from a logistic regression with no intercept.

- var.beta.hat:

  Variance-covariance matrix from a logistic regression with no
  intercept.

- degf:

  Degrees of freedom.

- which:

  Numeric vector indicating which parameters to compare, so that
  `RR = compare[2] / compare[1]`

- pf:

  Estimate *RR* or its complement *PF*?

- norm:

  Estimate confidence interval using quantiles of Guassian rather

  than t distribution quantiles?

- alpha:

  Complement of the confidence level.

- rnd:

  Number of digits for rounding. Affects display only, not estimates.

## Value

A [rror](https://abs-dev.github.io/PF/reference/rrorclass.md) object
with the following fields.

- `estimate`: vector with point and interval estimate

- `estimator`: either *PF* or *RR*

- `mu`: matrix with rows giving probability estimates for each of the
  groups

- `rnd`: how many digits to round the display

- `alpha`: complement of confidence level

- `norm`: logical indicating Gaussian or t-interval

- `degf`: degrees of freedom

## Details

Estimates confidence intervals using the delta method on parameters from
a generalized linear model with logit link.

\\RR = {{{\mu}}\_{2}} / {{{\mu}}\_{1}}\\, where \\{\mu}\_{i}\\ are the
estimated probabilities from the model.

## Note

Call to this function may be one of two formats: (1) specify `fit` or
(2) `beta.hat`, `var.beta.hat`, `degf`

`RRor(fit, degf = NULL, pf = TRUE, alpha = 0.05, which = c(1, 2), norm = TRUE, rnd = 3)`

`RRor(beta.hat, var.beta.hat, degf, pf = TRUE, alpha = 0.05, which = c(1, 2), norm = TRUE, rnd = 3)`

## See also

[rror](https://abs-dev.github.io/PF/reference/rrorclass.md),
[phiWt](https://abs-dev.github.io/PF/reference/phiWt.md),
[tauWt](https://abs-dev.github.io/PF/reference/tauWt.md)
[StatWI007](https://www.aphis.usda.gov/animal_health/vet_biologics/publications/STATWI0007.pdf)
for more examples

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
bird.fit <- glm(cbind(y, n - y) ~ tx - 1, binomial, bird)
RRor(tauWt(bird.fit))
#> 
#> 95% t intervals on 4 df
#> 
#> PF 
#>     PF     LL     UL 
#>  0.500 -0.583  0.842 
#> 
#>       mu.hat    LL     UL
#> txcon  0.733 0.943 0.3121
#> txvac  0.367 0.752 0.0997

# 95% t intervals on 4 df
#
# PF
#     PF     LL     UL
#  0.500 -0.583  0.842
#
#       mu.hat    LL     UL
# txcon  0.733 0.943 0.3121
# txvac  0.367 0.752 0.0997

RRor(phiWt(bird.fit))
#> 
#> 95% t intervals on 4 df
#> 
#> PF 
#>     PF     LL     UL 
#>  0.500 -0.583  0.842 
#> 
#>       mu.hat    LL     UL
#> txcon  0.733 0.943 0.3121
#> txvac  0.367 0.752 0.0997
# 95% t intervals on 4 df
#
# PF
#     PF     LL     UL
#  0.500 -0.583  0.842
#
#       mu.hat    LL     UL
# txcon  0.733 0.943 0.3121
# txvac  0.367 0.752 0.0997

```
