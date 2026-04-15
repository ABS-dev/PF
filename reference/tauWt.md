# Binomial dispersion: intra-cluster correlation parameter.

MME estimates of binomial dispersion parameter tau (intra-cluster
correlation).

## Usage

``` r
tauWt(
  fit,
  subset.factor = NULL,
  fit.only = TRUE,
  iter.max = 12,
  converge = 1e-06,
  trace.it = FALSE
)
```

## Arguments

- fit:

  A [glm](https://rdrr.io/r/stats/glm.html) object.

- subset.factor:

  Factor for estimating phi by subset. Will be converted to a factor if
  it is not a factor.

- fit.only:

  Return only the final fit? If FALSE, also returns the weights and tau
  estimates.

- iter.max:

  Maximum number of iterations.

- converge:

  Convergence criterion: difference between model degrees of freedom and
  Pearson's chi-square. Default 1e-6.

- trace.it:

  Display print statements indicating progress

## Value

A list with the following elements. `fit`: the new model fit, updated by
the estimated weights `weights`: vector of weights `phi`: vector of phi
estimates

## Details

Estimates binomial dispersion parameter \\\tau\\ by the method of
moments. Iteratively refits the model by the Williams procedure,
weighting the observations by \\1 / \phi\_{ij}\\, where \\\phi\_{ij} =
1 + \tau_j(n\_{ij} - 1)\\, \\j\\ indexes the subsets, and \\i\\ indexes
the observations.

## References

Williams DA, 1982. Extra-binomial variation in logistic linear models.
*Applied Statistics* 31:144-148.

Wedderburn RWM, 1974. Quasi-likelihood functions, generalized linear
models, and the Gauss-Newton method. *Biometrika* 61:439-447.

## See also

[phiWt](https://abs-dev.github.io/PF/reference/phiWt.md),
[RRor](https://abs-dev.github.io/PF/reference/RRor.md).

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
birdm.fit <- glm(cbind(y, n - y) ~ tx - 1, binomial, birdm)
RRor(tauWt(birdm.fit))
#> 
#> 95% t intervals on 4 df
#> 
#> PF 
#>     PF     LL     UL 
#>  0.489 -0.578  0.835 
#> 
#>       mu.hat    LL    UL
#> txcon  0.737 0.944 0.320
#> txvac  0.376 0.758 0.104

# 95% t intervals on 4 df
#
# PF
#     PF     LL     UL
#  0.489 -0.578  0.835
#
#       mu.hat    LL    UL
# txcon  0.737 0.944 0.320
# txvac  0.376 0.758 0.104
#
# binomial family only
# any link
```
