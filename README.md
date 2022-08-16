GET: Global envelopes
=====================

https://cran.r-project.org/package=GET

The `R` package `GET` provides global envelopes which can be used for central regions of functional or multivariate data (e.g. outlier detection, functional boxplot), for graphical Monte Carlo and permutation tests where the test statistic is a multivariate vector or function (e.g. goodness-of-fit testing for point patterns and random sets, functional ANOVA, functional GLM, n-sample test of correspondence of distribution functions), and for global confidence and prediction bands (e.g. confidence band in polynomial regression, Bayesian posterior prediction).

## The development version

The github repository holds a copy of the current development version of the contributed R package `GET`.

This development version is as or more recent than the official release of `GET` on the Comprehensive R Archive Network (CRAN) at https://cran.r-project.org/package=GET

## Where is the official release?

For the most recent official release of `GET`, see https://cran.r-project.org/package=GET

## Installation

### Installing the official release

To install the official release of `GET` from CRAN, start `R` and type

```R
install.packages('GET')
```

### Installing the development version

The easiest way to install the `GET` library from github is through the `remotes` package. Start `R` and type:

```R
require(remotes)
install_github('myllym/GET')
```
If you do not have the R library `remotes` installed, install it first by running

```R
install.packages("remotes")
```

After installation, in order to start using `GET`, load it to R and see
the main help page, which describes the functions of the library:
```R
require(GET)
help('GET-package')
```

If you want to have also vignettes working, you should also install packages from the 'suggests' field (`fda` and `fda.usc`),
have MiKTeX on your computer, and install the library with
```R
install_github('myllym/GET', build_vignettes = TRUE)
```

## Vignettes

The package contains two vignettes. The GET vignette describes the package in general. It is available by starting `R` and typing
```R
library("GET")
vignette("GET")
```
This vignette corresponds to Myllymäki and Mrkvička (2020).

The package provides also a vignette for global envelopes for point pattern analyses, which is available by starting `R` and typing
```R
library("GET")
vignette("pointpatterns")
```

## Branches

Currently two branches are provided in the development version. The main branch of GET is called `master`.

The other branch is called `FDR` and it includes also the experimental FDR envelopes tested in
Mrkvička and Myllymäki (2022, False discovery rate envelopes. arXiv:2008.10108 [stat.ME]).
The main branch includes the FDR envelopes which were found to have good performance in
Mrkvička and Myllymäki (2022).


## References

To cite GET in publications use

Myllymäki, M. and Mrkvička, T. (2020). GET: Global envelopes in R. arXiv:1911.06583 [stat.ME] URL: https://arxiv.org/abs/1911.06583

Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017).
Global envelope tests for spatial processes. Journal of the Royal Statistical Society:
Series B (Statistical Methodology) 79: 381-404. doi: 10.1111/rssb.12172 http://dx.doi.org/10.1111/rssb.12172
(You can find the preprint version of the article here: http://arxiv.org/abs/1307.0239v4)

and a suitable selection of:

Mrkvička, T., Myllymäki, M. and Hahn, U. (2017).
Multiple Monte Carlo testing, with applications in spatial point processes.
Statistics and Computing 27 (5): 1239-1255. https://doi.org/10.1007/s11222-016-9683-9

Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2020).
A one-way ANOVA test for functional data with graphical interpretation.
Kybernetika 56 (3), 432-458. http://doi.org/10.14736/kyb-2020-3-0432

Mrkvička, T., Myllymäki, M., Kuronen, M. and Narisetty, N. N. (2022).
New methods for multiple testing in permutation inference for the general linear model.
Statistics in Medicine 41(2), 276-297. https://doi.org/10.1002/sim.9236

Mrkvička and Myllymäki (2022).
False discovery rate envelopes. arXiv:2008.10108 [stat.ME]

Mrkvička, T., Roskovec, T. and Rost, M. (2021).
A nonparametric graphical tests of significance in functional GLM.
Methodology and Computing in Applied Probability 23, 593-612. https://doi.org/10.1007/s11009-019-09756-y

Mrkvička, T., Soubeyrand, S., Myllymäki, M., Grabarnik, P., and Hahn, U. (2016).
Monte Carlo testing in spatial statistics, with applications to spatial residuals.
Spatial Statistics 18, Part A: 40--53. https://doi.org/10.1016/j.spasta.2016.04.005

Myllymäki, M., Grabarnik, P., Seijo, H., and Stoyan, D. (2015).
Deviation test construction and power comparison for marked spatial point
patterns. Spatial Statistics 11: 19-34. https://doi.org/10.1016/j.spasta.2014.11.004
(You can find the preprint version of the article here: http://arxiv.org/abs/1306.1028)

Myllymäki, M., Kuronen, M. and Mrkvička, T. (2020).
Testing global and local dependence of point patterns on covariates in parametric models.
Spatial Statistics 42, 100436, https://doi.org/10.1016/j.spasta.2020.100436

Dai, W., Athanasiadis, S. and Mrkvička, T. (2022).
A new functional clustering method with combined dissimilarity sources and graphical interpretation.
Intech open. https://doi.org/10.5772/intechopen.100124
