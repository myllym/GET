\name{deviation_test}
\alias{deviation_test}
\title{Deviation test}
\usage{
  deviation_test(curve_set, r_min = NULL, r_max = NULL,
    use_theo = TRUE, scaling = "qdir", measure = "max",
    savedevs = FALSE)
}
\arguments{
  \item{scaling}{The name of the scaling to use. Options
  include 'none', 'q', 'qdir' and 'st'. 'qdir' is default.}

  \item{savedevs}{Logical. Should the global rank values
  k_i, i=1,...,nsim+1 be returned? Default: FALSE.}

  \item{curve_set}{A curve_set or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}

  \item{r_min}{The minimum radius to include.}

  \item{r_max}{The maximum radius to include.}

  \item{use_theo}{Whether to use the theoretical summary
  function or the mean of the simulations.}

  \item{measure}{The deviation measure to use. Default is
  'max'. Must be one of the following: 'max', 'int',
  'int2', 'MBD', 'MHRD', 'MMHRD', 'ISD' or 'MISD'.}
}
\value{
  If 'savedevs=FALSE' (default), the p-value is returned.
  If 'savedevs=TRUE', then a list containing the p-value
  and calculated deviation measures u_i, i=1,...,nsim+1
  (where u_1 corresponds to the data pattern) is returned.
}
\description{
  Crop the curve set to the interval of distances [r_min,
  r_max], calculate residuals, scale the residuals and
  perform a deviation test with a chosen deviation measure.
}
\details{
  The deviation test is based on a test function T(r) and
  it works as follows:

  1) The test function estimated for the data, T_1(r), and
  for nsim simulations from the null model, T_2(r), ....,
  T_{nsim+1}(r), must be saved in 'curve_set' and given to
  the deviation_test function.

  2) The deviation_test function then \itemize{ \item Crops
  the functions to the chosen range of distances [r_min,
  r_max].  \item If the curve_set is not consist of
  residuals, i.e. curve_set$is_residual is FALSE, then the
  residuals d_i(r) = T_i(r) - T_0(r) are calculated, where
  T_0(r) is the expectation of T(r) under the null
  hypothesis.  If use_theo = TRUE, as T_0(r) the
  theoretical value given in the curve_set$theo is used, if
  it is given. Otherwise, T_0(r) is estimated by the mean
  of T_j(r), j=2,...,nsim+1.  \item Scales the residuals.
  Options are \itemize{ \item 'none' No scaling. Nothing
  done.  \item 'q' Quantile scaling.  \item 'qdir'
  Directional quantile scaling.  \item 'st' Studentised
  scaling.  } \item Calculates the global deviation measure
  u_i, i=1,...,nsim+1, see options for 'measure'.  \item
  Calculates the p-value. }
}
\examples{
## Testing complete spatial randomness (CSR)
#-------------------------------------------
require(spatstat)
pp <- unmark(spruces)
Generate nsim simulations under CSR, calculate L-function for the data and simulations
env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
# The deviation test using the integral deviation measure
res <- deviation_test(env, measure = 'int')
res
# or
res <- deviation_test(env, r_min=0, r_max=7, measure='int')

## Random labeling test
#----------------------
mpp <- spruces
# T(r) = \\hat{L}_m(r), an estimator of the L_m(r) function
curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=4999, r_min=0, r_max=9.5)
res <- deviation_test(curve_set, measure='int')
res
}
\references{
  Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D.
  (2013). Deviation test construction and power comparison
  for marked spatial point patterns. arXiv:1306.1028
}

