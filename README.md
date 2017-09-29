# GET

An R library for Global Envelope Tests

## Installation

If you do not have the R library `devtools` installed, install it first by

```R
install.packages("devtools")
```
Then install the `GET` library with the following two R commands:

```R
library(devtools)
install_github('myllym/GET')
```

After installation, in order to start using `GET`, load it to R and see
the main help page, which describes the usage of the functions of the library:
```R
library(GET)
help(GET) # or help(GET, help="html")
```

In order to use the function random_labelling, the R library `marksummary` is
needed. It is available at https://github.com/myllym/marksummary.

## Branches

The branch for public use is called `master`. There are no other public branches at the moment.
The branch 'no_fastdepth' of the library `spptest` was taken as the master branch of `GET` September 21, 2016.
The library `spptest` has become `GET`, which is developed further.

## References

####To cite GET in publications use:

Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017).
Global envelope tests for spatial processes. Journal of the Royal Statistical Society:
Series B (Statistical Methodology) 79: 381-404. doi: 10.1111/rssb.12172 http://dx.doi.org/10.1111/rssb.12172
(You can find the preprint of the article here: http://arxiv.org/abs/1307.0239v4)

Mrkvička, T., Myllymäki, M. and Hahn, U. (2017).
Multiple Monte Carlo testing, with applications in spatial point processes.
Statistics and Computing 27 (5): 1239-1255. https://doi.org/10.1007/s11222-016-9683-9

####If you use deviation tests, please also cite:

Myllymäki, M., Grabarnik, P., Seijo, H., and Stoyan, D. (2015).
Deviation test construction and power comparison for marked spatial point
patterns. Spatial Statistics 11: 19-34. https://doi.org/10.1016/j.spasta.2014.11.004
(You can find the preprint of the article here: http://arxiv.org/abs/1306.1028)

####If you use apply tests to spatial residuals, please also cite:

Mrkvička, T., Soubeyrand, S., Myllymäki, M., Grabarnik, P., and Hahn, U. (2016).
Monte Carlo testing in spatial statistics, with applications to spatial residuals.
Spatial Statistics 18, Part A: 40--53. https://doi.org/10.1016/j.spasta.2016.04.005

####If you use functional ANOVA, please also cite:

Mrkvička, T., Hahn, U. and Myllymäki, M. (2016).
A one-way ANOVA test for functional data with graphical interpretation.
arXiv:1612.03608 [stat.ME] (http://arxiv.org/abs/1612.03608)