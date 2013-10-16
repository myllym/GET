\name{random_labelling_test}
\alias{random_labelling_test}
\title{Random labeling test}
\usage{
  random_labelling_test(pattern, mtf_name = "m",
    nsim = 999L, r_min = NULL, r_max = NULL, r_vec = NULL,
    measure = "max", scaling = "qdir", use_L = TRUE,
    edge_correction = "translate", use_theo = TRUE,
    method = "permute", use_biased_lambda2 = FALSE,
    ties = "random")
}
\arguments{
  \item{pattern}{A \code{\link[spatstat]{ppp}} object as
  the simple marked point pattern to be analysed. The marks
  need to be in the form of a numeric vector. The window
  has to have the type "rectangle".}

  \item{mtf_name}{A vector of mark test function names. "1"
  stands for the unmarked K-function. Accepted values are
  '1', 'm', 'mm', 'gamma', 'gammaAbs' and 'morAbs'.}

  \item{r_max}{A positive scalar value representing the
  maximum radius that should be considered. r_vec overrides
  r_max. By default, r_max is NULL and will get a sensible
  default.}

  \item{r_vec}{A monotonically increasing vector of
  non-negative r-values to act as the endpoints of the bins
  for the K_f-functions. r_vec overrides r_max. The bins
  are exclusive on the left and inclusive on the right. If
  the first vector element has value zero, it will be
  regarded as the collapsed bin [0, 0], and the next bin
  will start from and exclude 0.}

  \item{edge_correction}{The name of the edge correction to
  be used. Options are 'translate', 'none' and 'periodic'.}

  \item{method}{The name of the method to create
  simulations under the null hypothesis. 'permute' results
  in permutations of the marks. Using 'sample' will sample
  the marks from the empirical mark distribution with
  replacement. 'permute' is the default.}

  \item{use_biased_lambda2}{A logical scalar on whether to
  use the biased or the unbiased (in the Poisson case)
  estimate of the intensity squared.}

  \item{r_min}{The minimum radius to include.}

  \item{use_theo}{Whether to use the theoretical summary
  function or the mean of the simulations.}

  \item{scaling}{The name of the scaling to use. Options
  include 'none', 'q', 'qdir' and 'st'. 'qdir' is default.}

  \item{measure}{The deviation measure to use. Default is
  'max'. Must be one of the following: 'max', 'int' or
  'int2'.}

  \item{ties}{The method to treat tied values with. If one
  or more of the elements of sim_vec are equal to obs, how
  should the rank of obs be determined? For 'conservative'
  the resulting p-value will be the highest possible. For
  'liberal' the p-value will be the lowest possible. For
  'random' the rank of the obs within the tied values is
  uniformly sampled so that the resulting p-value is at
  most the conservative option and at least the liberal
  option. For 'midrank' the mid-rank within the tied values
  is taken. 'midrank' is the default.}
}
\description{
  Do a random labelling deviation test for a pattern using
  K_f- or L_f-functions.
}
\details{
  Given a marked point pattern, the function performs
  simulations under the random labelling hypothesis,
  estimates K_f- or L_f-functions for the pattern and
  simulations and makes a deviation test.
}
\examples{
require(spatstat)
mpp <- spruces
# T(r) = \\hat{L}_m(r), an estimator of the L_m(r) function
p <- random_labelling_test(mpp, mtf_name = 'm', nsim=999, r_min=0, r_max=9.5)
p
}

