# Internal function.

Internal function.

## Usage

``` r
.rr.score.asymp(y, alpha = 0.05, iter.max = 18, converge = 1e-04, mn = FALSE)
```

## Arguments

- y:

  data

- alpha:

  alpha

- iter.max:

  maximum number of iterations

- converge:

  convergence criterion

- mn:

  boolean whether to calculate MN or use default value of 1.0

## Examples

``` r
.rr.score.asymp(c(0, 18, 16, 19), mn = FALSE)
#>     point        LL        UL 
#> 0.0000000 0.0000000 0.2105332 
.rr.score.asymp(c(0, 18, 16, 19), mn = TRUE)
#>     point        LL        UL 
#> 0.0000000 0.0000000 0.2154228 
```
