# GET: Global envelopes in R

The R package provides global envelopes which can be used for central regions of functional or multivariate data (e.g. outlier detection, functional boxplot), for graphical Monte Carlo and permutation tests where the test statistic is a multivariate vector or function (e.g. goodness-of-fit testing for point patterns and random sets, functional ANOVA, functional GLM, n-sample test of correspondence of distribution functions), and for global confidence and prediction bands (e.g. confidence band in polynomial regression, Bayesian posterior prediction).

## Installation

You can install the `GET` library from github through `remotes` package with the following two R commands:

```R
library(remotes)
install_github('myllym/GET')
```
If you do not have the R library `remotes` installed, install it first by

```R
install.packages("remotes")
```

After installation, in order to start using `GET`, load it to R and see
the main help page, which describes the functions of the library:
```R
library(GET)
help('GET-package')
```

## Branches

The branch for public use is called `master`. There are no other public branches at the moment.
The branch 'no_fastdepth' of the library `spptest` was taken as the master branch of `GET` September 21, 2016.
The library `spptest` has become `GET`, which is developed further.

## References

To cite GET in publications use (a selection of):

Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017).
Global envelope tests for spatial processes. Journal of the Royal Statistical Society:
Series B (Statistical Methodology) 79: 381-404. doi: 10.1111/rssb.12172 http://dx.doi.org/10.1111/rssb.12172
(You can find the preprint of the article here: http://arxiv.org/abs/1307.0239v4)

Mrkvička, T., Myllymäki, M. and Hahn, U. (2017).
Multiple Monte Carlo testing, with applications in spatial point processes.
Statistics and Computing 27 (5): 1239-1255. https://doi.org/10.1007/s11222-016-9683-9

Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2018).
A one-way ANOVA test for functional data with graphical interpretation.
arXiv:1612.03608 [stat.ME] (http://arxiv.org/abs/1612.03608)

Mrkvička, T., Myllymäki, M. and Narisetty, N. N. (2019).
New methods for multiple testing in permutation inference for the general linear model. arXiv:1906.09004 [stat.ME]

Mrkvička, T., Roskovec, T. and Rost, M. (2019).
A nonparametric graphical tests of significance in functional GLM. arXiv:1902.04926 [stat.ME]

Myllymäki, M., Grabarnik, P., Seijo, H., and Stoyan, D. (2015).
Deviation test construction and power comparison for marked spatial point
patterns. Spatial Statistics 11: 19-34. https://doi.org/10.1016/j.spasta.2014.11.004
(You can find the preprint of the article here: http://arxiv.org/abs/1306.1028)

Mrkvička, T., Soubeyrand, S., Myllymäki, M., Grabarnik, P., and Hahn, U. (2016).
Monte Carlo testing in spatial statistics, with applications to spatial residuals.
Spatial Statistics 18, Part A: 40--53. https://doi.org/10.1016/j.spasta.2016.04.005
