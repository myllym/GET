# GET

An R library for Global Envelope Tests

## Installation

Install the library to R with the following two R commands:

```R
library(devtools)
install_github('myllym/GET')
```

If you do not have the library ´devtools´ installed, install it first by

```R
install.packages("devtools")
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


## References

Myllymäki, M., Grabarnik, P., Seijo, H., and Stoyan, D. (2015).
Deviation test construction and power comparison for marked spatial point
patterns. Spatial Statistics 11, 19-34.
(Preprint of the article: http://arxiv.org/abs/1306.1028)

Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2016).
Global envelope tests for spatial processes. Journal of the Royal Statistical Society:
Series B (Statistical Methodology). doi: 10.1111/rssb.12172
(Preprint of the article: http://arxiv.org/abs/1307.0239v4)

Mrkvička, T., Myllymäki, M. and Hahn, U. (2016).
Multiple Monte Carlo testing, with applications in spatial point processes.
Statistics and Computing. doi: 10.1007/s11222-016-9683-9

Mrkvička, T., Soubeyrand, S., Myllymäki, M., Grabarnik, P., and Hahn, U. (2016).
Monte Carlo testing in spatial statistics, with applications to spatial residuals.
Spatial Statistics. doi: http://dx.doi.org/10.1016/j.spasta.2016.04.005