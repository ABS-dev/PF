# Binomial dispersion parameter.

MME estimate of dispersion parameter phi.

## Usage

``` r
phiWt(fit, subset.factor = NULL, fit.only = TRUE, show.warns = FALSE)
```

## Arguments

- fit:

  A [glm](https://rdrr.io/r/stats/glm.html) object.

- subset.factor:

  Factor for estimating phi by subset. Will be converted to a factor if
  it is not a factor.

- fit.only:

  Return only the new fit? If FALSE, also returns the weights and phi
  estimates.

- show.warns:

  Show warnings

## Value

A list with the following elements. `fit`: the new model fit, updated by
the estimated weights `weights`: vector of weights `phi`: vector of phi
estimates

## Details

Estimates binomial dispersion parameter \\\phi\\ by the method of
moments. Refits the model, weighting the observations by \\1/\phi\\.
Uses `quasibinomial` family in
[`glm()`](https://rdrr.io/r/stats/glm.html).

## References

Wedderburn RWM, 1974. Quasi-likelihood functions, generalized linear
models, and the Gauss-Newton method. *Biometrika* 61:439-447.

## See also

[tauWt](https://abs-dev.github.io/PF/reference/tauWt.md),
[RRor](https://abs-dev.github.io/PF/reference/RRor.md).

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
birdm.fit <- glm(cbind(y, n - y) ~ tx-1, binomial, birdm)
RRor(phiWt(birdm.fit))
#> 
#> 95% t intervals on 4 df
#> 
#> PF 
#>     PF     LL     UL 
#>  0.479 -0.537  0.823 
#> 
#>       mu.hat   LL    UL
#> txcon  0.768 0.95 0.367
#> txvac  0.400 0.78 0.111
#
# 95% t intervals on 4 df
#
# PF
#     PF     LL     UL
#  0.479 -0.537  0.823
#
#       mu.hat   LL    UL
# txcon  0.768 0.95 0.367
# txvac  0.400 0.78 0.111
#
```
