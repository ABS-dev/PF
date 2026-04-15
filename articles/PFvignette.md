# PF Vignette

## 1 Introduction

The `PF` package is a collection of functions related to estimating
prevented fraction, \\PF=1-RR\\ , where
\\RR={{{\pi}\_{2}}}/{{{\pi}\_{1}}}\\\\ .

### 1.1 Technical notes

*Optimization*. Unless otherwise stated, optimization is by the DUD
algorithm (Ralston and Jennrich 1978).

*Level tested*. The help files indicate the level of testing undergone
by each function. In some cases that is a subjective judgement, since
most of these functions were originally tested in SPlus and have been
ported to R more recently.

## 2 Score based methods

### 2.1 The score statistic

Confidence intervals for the risk ratio may be based on the score
statistic Miettinen and Nurminen (1985), \\\frac{{{{\hat{\pi
}}}\_{2}}-{{\rho }\_{0}}{{{\hat{\pi }}}\_{1}}}{\sqrt{\left(
{{{{\tilde{\pi }}}\_{2}}(1-{{{\tilde{\pi }}}\_{2}})}/{{{n}\_{2}}}\\
\right)+\rho \_{0}^{2}\left( {{{{\tilde{\pi }}}\_{1}}(1-{{{\tilde{\pi
}}}\_{1}})}/{{{n}\_{1}}}\\ \right)}}\\ where hat indicates the MLE and
tilde indicates the MLE under the restriction that \\\rho ={{\rho
}\_{0}}\\.

### 2.2 Asymptotic intervals

[`RRsc()`](https://abs-dev.github.io/PF/reference/RRsc.md) estimates
asymptotic confidence intervals for the risk ratio or prevented fraction
based on the score statistic. Interval estimates are returned for three
estimators. The score method was originally introduced by (Koopman
1984). Gart and Nam’s modification includes a skewness correction (Gart
and Nam 1988). The method of (Miettinen and Nurminen 1985) is a version
made slightly more conservative than Koopman’s by including a factor of
\\{(N-1)}/{N}\\\\.

``` r
require(PF)
```

    ## Loading required package: PF

``` r
RRsc(c(4, 24, 12, 28))
```

    ## 
    ## PF 
    ## 95% interval estimates
    ## 
    ##                 PF     LL    UL
    ## MN method    0.611 0.0251 0.857
    ## score method 0.611 0.0328 0.855
    ## skew corr    0.611 0.0380 0.876

Starting estimates for the algorithm are obtained by the modified Katz
method (log method with 0.5 added to each cell). Both forms of the Katz
estimate may be retrieved from the returned object using
`RRsc()\$estimate`.

### 2.3 Exact intervals

These methods give intervals that are \`exact’ in the sense that they
are based on the actual sampling distribution rather than an
approximation to it. The score statistic is used to select the \\2
\times 2\\ tables that would fall in the tail area, and the binomial
probability is estimated over the tail area by taking the maximum over
the nuisance parameter. The search over the tail area is made more
efficient by the Berger-Boos correction (Berger and Boos 1994).

[`RRtosst()`](https://abs-dev.github.io/PF/reference/RRtosst.md) gives
intervals obtained by inverting two one-sided score tests;
[`RRotsst()`](https://abs-dev.github.io/PF/reference/RRotsst.md) gives
intervals obtained by inverting one two-sided score test.
[`RRtosst()`](https://abs-dev.github.io/PF/reference/RRtosst.md) is thus
more conservative, preserving at least \\{\alpha }/{2}\\\\ in each tail.
(Agresti and Min 2001) discuss the properties and relative benefits of
the two approaches. The price of exactnesss is conservatism, due to the
discreteness of the binomial distribution (Agresti 2001). This means
that the actual coverage of the confidence interval does not exactly
conform to the nominal coverage, but it will not be less than it. (See
also (Agresti 2003).) Both functions use a simple step search algorithm.

``` r
RRotsst(c(4, 24, 12, 28))
```

    ## 
    ## PF 
    ## 95% interval estimates
    ## 
    ##     PF     LL     UL 
    ## 0.6111 0.0148 0.8519

``` r
RRtosst(c(4, 24, 12, 28))
```

    ## 
    ## PF 
    ## 95% interval estimates
    ## 
    ##    PF    LL    UL 
    ## 0.611 0.012 0.902

## 3 Stratified designs

Methods for estimating a common *RR* from stratified or clustered
designs depend on homogeneity with respect to the common parameter.

### 3.1 Gart-Nam score method

(Gart 1985) and (Gart and Nam 1988) derived a score statistic for a
common estimate of *RR* from designs with multiple independent strata,
and they showed that it is identical to one proposed by (Radhakrishnan
1965) from a different perspective.

[`RRstr()`](https://abs-dev.github.io/PF/reference/RRstr.md) provides
confidence intervals and a homogeneity test based on Gart’s statistic.

Data may be input two ways, either using a formula and data frame, or as
a matrix.

``` r
RRstr(cbind(y, n) ~ tx + cluster(clus), Table6,
      compare = c("a", "b"), pf = FALSE)
```

    ## 
    ## Test of homogeneity across clusters
    ## 
    ## stat  0.954 
    ## df    3 
    ## p     0.812 
    ## 
    ## RR 
    ## 95% interval estimates
    ## 
    ##             RR   LL   UL
    ## starting  2.66 1.37 5.18
    ## mle       2.65 1.39 5.03
    ## skew corr 2.65 1.31 5.08

``` r
# Data matrix input: RRstr(Y = table6, pf = F)
```

### 3.2 Mantel-Haenszel estimator

A widely-used heuristic method for sparse frequency tables is the
weighted average approach of (Mantel and Haenszel 1959)[^1]. MH interval
estimates are based on the asymptotic normality of the log of the risk
ratio. [`RRmh()`](https://abs-dev.github.io/PF/reference/RRmh.md)
utilizes the variance estimator given by (Greenland and Robins 1985) for
sparse strata. The resulting asymptotic estimator is consistent for both
the case of sparse strata where the number of strata is assumed
increasing, and the case of limited number of strata where the stratum
size is assumed increasing. In the latter case, however, it is less
efficient than maximum likelihood Greenland and Robins (1985).
Additional discussion may be found in Section 4.3.1 (Lachin 2000),
(Landis et al. 2005), and (Somes and O’Brien 2006) [^2].

``` r
RRmh(cbind(y, n) ~ tx + cluster(clus), Table6,
     compare = c("a", "b"), pf = FALSE)
```

    ## 
    ## RR 
    ## 95% interval estimates
    ## 
    ##   RR   LL   UL 
    ## 2.67 1.37 5.23

``` r
# Data matrix input: RRmh(Y = table6, pf = F)
```

### 3.3 Examples

For a fuller set of examples, see the vignette *Examples for Stratified
Designs*.

## 4 Model based intervals

### 4.1 Logistic regression estimates

Intervals may be estimated from logistic regression models with
[`RRor()`](https://abs-dev.github.io/PF/reference/RRor.md). It takes the
fit of a [`glm()`](https://rdrr.io/r/stats/glm.html) object and
estimates the intervals by the delta method.

``` r
bird.fit <- glm(cbind(y, n - y) ~ tx - 1, binomial, bird)
RRor(bird.fit)
```

    ## 
    ## 95% t intervals on 4 df
    ## 
    ## PF 
    ##      PF      LL      UL 
    ##  0.5000 -0.0406  0.7598 
    ## 
    ##       mu.hat    LL    UL
    ## txcon  0.733 0.896 0.466
    ## txvac  0.367 0.624 0.168

### 4.2 Estimating the dispersion parameter

The binomial GLM weights are \\\frac{\hat{\pi }(1-\hat{\pi
})}{a(\hat{\varphi })/n\\}\\ where \\a(\hat{\varphi })\\ is a function
of the dispersion parameter.

#### 4.2.1 Dispersion parameter \\\varphi\\

A simple estimator of the dispersion parameter, \\\varphi\\, may be
estimated by the method of moments \[Wed74\]. It is given by
[`phiWt()`](https://abs-dev.github.io/PF/reference/phiWt.md). This form
of the dispersion parameter has \\a(\varphi)=\varphi\\, and \\\varphi\\
is estimated by \\{X^2}/{df}\\\\, the Pearson statistic divided by the
degrees of freedom.

Note that \\\varphi\\ is the same estimator as may be obtained by the
`quasibinomial` family in [`glm()`](https://rdrr.io/r/stats/glm.html)
which is, in fact, what is used by
[`phiWt()`](https://abs-dev.github.io/PF/reference/phiWt.md) to reweight
the original fit:

``` r
phiWt(bird.fit, fit.only = FALSE)$phi
```

    ##      all 
    ## 2.471592

``` r
summary(update(bird.fit, family = quasibinomial))$disp
```

    ## [1] 2.471592

[`phiWt()`](https://abs-dev.github.io/PF/reference/phiWt.md) makes it
easy to estimate *PF* intervals with a single command.

``` r
# model weighted by phi hat
RRor(phiWt(bird.fit))
```

    ## 
    ## 95% t intervals on 4 df
    ## 
    ## PF 
    ##     PF     LL     UL 
    ##  0.500 -0.583  0.842 
    ## 
    ##       mu.hat    LL     UL
    ## txcon  0.733 0.943 0.3121
    ## txvac  0.367 0.752 0.0997

It also allows different estimates of \\{\hat{\varphi}}\\ for specified
subsets of the data.

``` r
# model with separate phi for vaccinates and controls
RRor(phiWt(bird.fit, subset.factor = bird$tx))
```

    ## 
    ## 95% t intervals on 4 df
    ## 
    ## PF 
    ##     PF     LL     UL 
    ##  0.500 -0.645  0.848 
    ## 
    ##       mu.hat    LL     UL
    ## txcon  0.733 0.938 0.3330
    ## txvac  0.367 0.767 0.0925

If you want to subtract a degree of freedom for each additional
parameter, you can do that by entering the degrees of freedom as an
argument to [`RRor()`](https://abs-dev.github.io/PF/reference/RRor.md).

``` r
# subtract 2 degrees of freedom
RRor(phiWt(bird.fit, subset.factor = bird$tx), degf = 2)
```

    ## 
    ## 95% t intervals on 2 df
    ## 
    ## PF 
    ##     PF     LL     UL 
    ##  0.500 -2.164  0.921 
    ## 
    ##       mu.hat    LL     UL
    ## txcon  0.733 0.975 0.1635
    ## txvac  0.367 0.895 0.0377

#### 4.2.2 Dispersion parameter \\\tau\\

When overdispersion is due to intra-cluster correlation, it may make
sense to estimate the dispersion as a function of the intra-cluster
correlation parameter \\\tau\\. In other words, \\{a({\varphi
}\_{ij})}=1+{{\tau }\_{j}}({{n}\_{ij}}-1)\\.
[`tauWt()`](https://abs-dev.github.io/PF/reference/tauWt.md) does this
using the Williams procedure (Williams 1982).

``` r
# model weighted using tau hat
RRor(tauWt(bird.fit, subset.factor = bird$tx))
```

    ## 
    ## 95% t intervals on 4 df
    ## 
    ## PF 
    ##     PF     LL     UL 
    ##  0.500 -0.645  0.848 
    ## 
    ##       mu.hat    LL     UL
    ## txcon  0.733 0.938 0.3330
    ## txvac  0.367 0.767 0.0925

In this example the
[`tauWt()`](https://abs-dev.github.io/PF/reference/tauWt.md) estimates
are the same as the
[`phiWt()`](https://abs-dev.github.io/PF/reference/phiWt.md) estimates.
That is because the cluster sizes are all the same. Let’s see what
happens if we modify the `bird` data set. The `birdm` data set has the
same cluster fractions but differing cluster sizes.

``` r
# different cluster sizes, same cluster fractions
birdm.fit <- glm(cbind(y, n - y) ~ tx - 1, binomial, birdm)
RRor(tauWt(birdm.fit, subset.factor = birdm$tx))
```

    ## 
    ## 95% t intervals on 4 df
    ## 
    ## PF 
    ##     PF     LL     UL 
    ##  0.490 -0.605  0.838 
    ## 
    ##       mu.hat    LL    UL
    ## txcon  0.737 0.942 0.328
    ## txvac  0.376 0.764 0.101

Note that increasing cluster size can make things worse when there is
intra-cluster correlation.

Now let’s compare the weights from
[`phiWt()`](https://abs-dev.github.io/PF/reference/phiWt.md) and
[`tauWt()`](https://abs-dev.github.io/PF/reference/tauWt.md) with
unequal cluster sizes. In the output below, `w` represents
\\1/a(\hat{\varphi })\\\\ and `nw` is \\n/a(\hat{\varphi })\\\\

``` r
# Compare phi and tau weights
#
phi.wts <- phiWt(birdm.fit, fit.only = FALSE, subset.factor = birdm$tx)$weights
tau.wts <- tauWt(birdm.fit, fit.only = FALSE, subset.factor = birdm$tx)$weights
w <- cbind(w.phi  = phi.wts,
           w.tau  = tau.wts,
           nw.phi = phi.wts * birdm$n,
           nw.tau = tau.wts * birdm$n)
print(cbind(birdm[, c(3, 1, 2)], round(w, 2)), row.names = FALSE)
```

    ##   tx  y  n w.phi w.tau nw.phi nw.tau
    ##  vac  1 10  0.32  0.35   3.20   3.55
    ##  vac  8 20  0.32  0.21   6.40   4.13
    ##  vac  9 15  0.32  0.26   4.80   3.92
    ##  con  8 16  0.21  0.27   3.39   4.33
    ##  con  8 10  0.21  0.38   2.12   3.82
    ##  con 27 30  0.21  0.16   6.36   4.84

Look at the last two rows. Note that the `nw.phi` are directly
proportional to `n` within treatment group, while the `nw.tau` are not.
With intra-cluster correlation, increasing cluster size does not give a
corresponding increase in information.

### 4.3 Rao-Scott weights

(Rao and Scott 1992) give a method of weighting clustered binomial
observations based on the variance inflation due to clustering. They
relate their approach to the concepts of design effect and effective
sample size familiar in survey sampling, and they illustrate its use in
a variety of contexts.
[`rsbWt()`](https://abs-dev.github.io/PF/reference/rsbWt.md) implements
it in the same manner as
[`phiWt()`](https://abs-dev.github.io/PF/reference/phiWt.md) and
[`tauWt()`](https://abs-dev.github.io/PF/reference/tauWt.md). For more
general use, the function
[`rsb()`](https://abs-dev.github.io/PF/reference/rsb.md) just returns
the design effect estimates and the weights.

``` r
# model weighted with Rao-Scott weights
RRor(rsbWt(birdm.fit, subset.factor = birdm$tx))
```

    ## 
    ## 95% t intervals on 4 df
    ## 
    ## PF 
    ##     PF     LL     UL 
    ##  0.479 -0.314  0.793 
    ## 
    ##       mu.hat    LL    UL
    ## txcon  0.768 0.960 0.311
    ## txvac  0.400 0.717 0.149

``` r
# just the design effect estimates
rsb(birdm$y, birdm$n)$d
```

    ##      all 
    ## 6.334514

## 5 Likelihood based intervals

The [`RRlsi()`](https://abs-dev.github.io/PF/reference/RRlsi.md)
function estimates likelihood support intervals for *RR* by the profile
likelihood Section 7.6 (Royall 1997).

Likelihood support intervals are usually formed based on the desired
likelihood ratio, \\{1}/{k}\\\\, often \\{1}/{8}\\\\ or \\{1}/{32}\\\\.
Under some conditions the log likelihood ratio may follow the chi-square
distribution. If so, then \\\alpha =1-{{F}\_{{{\chi }^{2}}}}\left( 2\log
(k),1 \right)\\.
[`RRsc()`](https://abs-dev.github.io/PF/reference/RRsc.md) will make the
conversion from \\\alpha\\ to \\k\\ with the argument `use.alpha = T`.

``` r
RRlsi(c(4, 24, 12, 28))
```

    ## 
    ## 1/8 likelihood support interval for PF 
    ## 
    ## corresponds to 95.858% confidence
    ##   (under certain assumptions)
    ## 
    ## PF 
    ##     PF     LL     UL 
    ## 0.6111 0.0168 0.8859

``` r
RRlsi(c(4, 24, 12, 28), use.alpha = TRUE)
```

    ## 
    ## 1/6.826 likelihood support interval for PF 
    ## 
    ## corresponds to 95% confidence
    ##   (under certain assumptions)
    ## 
    ## PF 
    ##     PF     LL     UL 
    ## 0.6111 0.0495 0.8792

## 6 Incidence ratio

The incidence is the number of cases per subject-time. Its distribution
is assumed Poisson. Under certain designs, the incidence ratio (\\IR\\)
is used as a measure of treatment effect. Correspondingly,
\\{{PF}\_{IR}}=1-IR\\ would be used as a measure of effect for an
intervention that is preventive, such as vaccination. \\IR\\ is also
called incidence density ratio (\\IDR\\), and that is the term used in
the following functions.

### 6.1 Score method

[`IDRsc()`](https://abs-dev.github.io/PF/reference/IDRsc.md) estimates a
confidence interval for the incidence density ratio using Siev’s formula
Appendix 1 (Siev 1994) based on the Poisson score statistic[^3].
\\IDR=\widehat{IDR}\left\\ 1+\left(
\frac{1}{{{y}\_{1}}}+\frac{1}{{{y}\_{2}}} \right)\frac{z\_{\alpha
/2}^{2}}{2}\\ \\ \pm \\ \\ \frac{z\_{\alpha
/2}^{2}}{2{{y}\_{1}}{{y}\_{2}}}\sqrt{{{y}\_{\bullet }}\left(
{{y}\_{\bullet }}z\_{\alpha /2}^{2}+4{{y}\_{1}}{{y}\_{2}} \right)}
\right\\\\

``` r
IDRsc(c(26, 204, 10, 205), pf = FALSE)
```

    ## 
    ## IDR 
    ## 95% interval estimates
    ## 
    ##  IDR   LL   UL 
    ## 2.61 1.28 5.34

### 6.2 Likelihood method

A likelihood support interval for \\IDR\\ may be estimated based on
orthogonal factoring of the reparameterized likelihood. Section 7.2
(Royall 1997)
[`IDRlsi()`](https://abs-dev.github.io/PF/reference/IDRlsi.md)
implements this method.

Likelihood support intervals are usually formed based on the desired
likelihood ratio, \\{1}/{k}\\\\, often \\{1}/{8}\\\\ or \\{1}/{32}\\\\.
Under some conditions the log likelihood ratio may follow the chi square
distribution. If so, then \\\alpha =1-{{F}\_{{{\chi }^{2}}}}\left( 2\log
(k),1 \right)\\.
[`IDRlsi()`](https://abs-dev.github.io/PF/reference/IDRlsi.md) will make
the conversion from \\\alpha\\ tp \\k\\ with the argument
`use.alpha = T`.

``` r
IDRlsi(c(26, 204, 10, 205), pf = FALSE)
```

    ## 
    ## 1/8 likelihood support interval for IDR 
    ## 
    ## corresponds to 95.858% confidence
    ##   (under certain assumptions)
    ## 
    ## IDR 
    ##  IDR   LL   UL 
    ## 2.61 1.26 5.88

``` r
IDRlsi(c(26, 204, 10, 205), pf = FALSE, use.alpha = TRUE)
```

    ## 
    ## 1/6.826 likelihood support interval for IDR 
    ## 
    ## corresponds to 95% confidence
    ##   (under certain assumptions)
    ## 
    ## IDR 
    ##  IDR   LL   UL 
    ## 2.61 1.30 5.69

## References

Agresti, A. 2001. “Exact Inference for Categorical Data: Recent Advances
and Continuing Controversies.” *Statistics in Medicine* 20: 2709–22.

———. 2003. “Dealing with Discreteness: Making ‘Exact’ Confidence
Intervals for Proportions, Differences of Proportions, and Odds Ratios
More Exact.” *Statistical Methods in Medical Research* 12: 3–21.

Agresti, A., and J. Hartzel. 2000. “Strategies for Comparing Treatments
on a Binary Response with Multi-Centre Data.” *Statistics in Medicine*
19: 1115–39.

Agresti, A., and Y. Min. 2001. “On Small-Sample Confidence Intervals for
Parameters in Discrete Distributions.” *Biometrics* 57: 963–71.

Berger, R. L., and D. D. Boos. 1994. “P Values Maximized over a
Confidence Set for the Nuisance Parameter.” *Journal of the American
Statistical Association* 89: 214–20.

Cochran, W. G. 1954. “Some Methods for Strengthening the Common
Chi-Square Tests.” *Biometrics* 10: 417–51.

Gart, J. J. 1985. “Approximate Tests and Interval Estimation of the
Common Relative Risk in the Combination of \\2 \times 2\\ Tables.”
*Biometrika* 72: 673–77.

Gart, J. J., and J. Nam. 1988. “Approximate Interval Estimation of the
Ratio of Binomial Parameters: A Review and Corrections for Skewness.”
*Biometrics* 44: 323–38.

Graham, P. L., K. Mengersen, and A. P. Morton. 2003. “Confidence Limits
for the Ratio of Two Rates Based on Likelihood Scores: Non-Iterative
Method.” *Statistics in Medicine* 22: 2071–83.

Greenland, S., and J. M. Robins. 1985. “Estimation of a Common Effect
Parameter from Sparse Follow-up Data.” *Biometrics* 41: 55–68.

Koopman, P. A. R. 1984. “Confidence Intervals for the Ratio of Two
Binomial Proportions.” *Biometrics* 40: 513–17.

Kuritz, S. J., J. R. Landis, and G. G. Koch. 1988. “A General Overview
of Mantel-Haenszel Methods: Applications and Recent Developments.”
*Annual Review of Public Health* 9: 123–60.

Lachin, J. M. 2000. *Biostatistical Methods: The Assessment of Relative
Risks*. New York: Wiley.

Landis, J. R., T. J. Sharp, S. J. Kuritz, and G. G. Koch. 2005.
“Mantel-Haenszel Methods.” In *Encyclopedia of Biostatistics*, edited by
P. Armitage and T. Colton, Second, 4:2937–50. Chichester: Wiley.

Mantel, N., and W. Haenszel. 1959. “Statistical Aspects of the Analysis
of Data from Retrospective Studies of Disease.” *Journal of the National
Cancer Institute* 22: 719–48.

Miettinen, O., and M. Nurminen. 1985. “Comparative Analysis of Two
Rates.” *Statistics in Medicine* 4: 213–26.

Radhakrishnan, S. 1965. “Combination of Results from Several \\2 \times
2\\ Contingency Tables.” *Biometrics* 21: 86–98.

Ralston, M. L., and R. I. Jennrich. 1978. “DUD, a Derivative-Free
Algorithm for Nonlinear Least Squares.” *Technometrics* 20: 7–14.

Rao, J. N. K., and A. J. Scott. 1992. “A Simple Method for the Analysis
of Clustered Binary Data.” *Biometrics* 48: 577–85.

Royall, R. 1997. *Statistical Evidence: A Likelihood Paradigm*. Boca
Raton: Chapman; Hall.

Sato, T. 1988. “Confidence Intervals for Effect Parameters from Cohort
Studies Based on Efficient Scores.” *Japanese Journal of Applied
Statistics* 17: 43–54.

Siev, D. 1994. “Estimating Vaccine Efficacy in Prospective Studies.”
*Preventive Veterinary Medicine* 20: 279–96.

———. 2004. “Letter to the Editor: Confidence Limits for the Ratio of Two
Rates Based on Likelihood Scores, Non-Iterative Method.” *Statistics in
Medicine* 23: 693.

Somes, G. W., and K. F. O’Brien. 2006. “Mantel-Haenszel Statistics.” In
*Encyclopedia of Statistical Sciences*, edited by S. Kotz, N.
Balakrishnan, C. B. Read, and B. Vidakovic, Second, 7:4483–86. Hoboken:
Wiley.

Williams, D. A. 1982. “Extra-Binomial Variation in Logistic Linear
Models.” *Applied Statistics* 31: 144–48.

[^1]: (Kuritz, Landis, and Koch 1988) review the Mantel-Haenzel approach
    and point out its relationship to a method proposed by (Cochran
    1954), which was the basis of Rhadakrishnan’s method (Radhakrishnan
    1965), alluded to in Section 3.1.

[^2]: SAS Proc FREQ provides MH interval estimates of *RR*. The other
    estimator calculated by Proc FREQ, which it calls \`\`logit,’’ is
    actually a weighted least squares estimator (Lachin 2000) that has a
    demonstrable and severe bias for sparse data (Greenland and Robins
    1985). It should be avoided.

[^3]: This formula was published in a Japanese journal (Sato 1988)
    several years before Siev. See also (Graham, Mengersen, and Morton
    2003) and (Siev 2004).
