# spptest

An R library for spatial point pattern testing

spptest provides envelope and deviation tests for spatial point processes.
The main motivation for a new package are the scalings for deviation tests
and global envelope tests.

## Installation

Install the library to R with the following two R commands:

```R
library(devtools)
install_github('spptest', username = 'myllym', ref = 'no_fastdepth')
```

If you do not have the library ´devtools´ installed, install it first by

```R
install.packages("devtools")
```

After installation, in order to start using `spptest`, load it to R and see
the main help page, which describes the usage of the functions of the library:
```R
library(spptest)
help(spptest) # or help(spptest, help="html")
```

In order to use the function random_labelling, the R library `marksummary` is
needed. It is currently available by request.

## Branches

The branch for current public use is called `no_fastdepth`. It includes all the
same features as the `master` branch. The only difference is that the `master` 
branch includes references to functional depth measures provided by the R 
library fastdepth which is not yet publicly available.


