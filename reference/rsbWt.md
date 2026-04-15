# Rao-Scott weighting.

Rao-Scott weighting of clustered binomial observations.

## Usage

``` r
rsbWt(fit = NULL, subset.factor = NULL, fit.only = TRUE)
```

## Arguments

- fit:

  A [stats::glm](https://rdrr.io/r/stats/glm.html) object.

- subset.factor:

  Factor for estimating phi by subset. Will be converted to a factor if
  it is not a factor.

- fit.only:

  Return only the new fit? If FALSE, also returns the weights and phi
  estimates.

## Value

A list with the following elements.

- `fit`: the new model fit, updated by the estimated weights

- `weights`: vector of weights

- `d`: vector of \\{d}\_{i}\\ estimates

## Details

Estimates the cluster design effect \\{d}\_{i}\\ as the variance
inflation due to clustering by the method of Rao and Scott. Observations
are then weighted by the inverse of the \\{d}\_{i}\\.

## References

Rao JNK, Scott AJ, 1992. A simple method for the analysis of clustered
binary data. *Biometrics* 48:577-585.

## See also

[RRor](https://abs-dev.github.io/PF/reference/RRor.md),
[rsb](https://abs-dev.github.io/PF/reference/rsb.md).

## Author

[PF-package](https://abs-dev.github.io/PF/reference/PF-package.md)

## Examples

``` r
birdm.fit <- glm(cbind(y, n - y) ~ tx-1, binomial, birdm)
RRor(rsbWt(birdm.fit))
#> 
#> 95% t intervals on 4 df
#> 
#> PF 
#>     PF     LL     UL 
#>  0.479 -1.061  0.868 
#> 
#>       mu.hat    LL     UL
#> txcon  0.768 0.968 0.2659
#> txvac  0.400 0.848 0.0737
#
# 95% t intervals on 4 df
#
# PF
#     PF     LL     UL
#  0.479 -1.061  0.868
#
#       mu.hat    LL     UL
# txcon  0.768 0.968 0.2659
# txvac  0.400 0.848 0.0737
#
```
