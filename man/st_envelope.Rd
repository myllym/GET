\name{st_envelope}
\alias{st_envelope}
\title{Studentised envelope test}
\usage{
  st_envelope(curve_set, alpha = 0.05, savedevs = FALSE,
    ...)
}
\arguments{
  \item{savedevs}{Logical. Should the global deviation
  values u_i, i=1,...,nsim+1 be returned? Default: FALSE.}

  \item{curve_set}{A curve_set (see
  \code{\link{create_curve_set}}) or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}

  \item{alpha}{The significance level. Simultaneous
  100(1-alpha) percent envelopes will be calculated.}

  \item{...}{Additional parameters passed to
  \code{\link{estimate_p_value}} to obtain a point estimate
  for the p-value.}
}
\value{
  An "envelope_test" object containing the following
  fields: \itemize{ \item r Distances for which the test
  was made.  \item method The envelope method.  \item p A
  point estimate for the p-value.  \item u_alpha \item u
  Global deviation values. Returned only if savedevs =
  TRUE.  \item central_curve The mean test function
  (median) calculated from simulations.  \item data_curve
  The test function for the data.  \item lower The lower
  envelope.  \item upper The upper envelope.  \item call
  The call of the function. }
}
\description{
  The studentised envelope test, which takes into account
  the unequal variances of the test function T(r) for
  different distances r.
}
\examples{
## Testing complete spatial randomness (CSR)
#-------------------------------------------
library(spatstat)
pp <- spruces
## Test for complete spatial randomness (CSR)
# Generate nsim simulations under CSR, calculate L-function for the data and simulations
env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
# The studentised envelope test
res <- st_envelope(env)
plot(res)
# or (requires R library ggplot2)
plot(res, use_ggplot2=T)

## Advanced use:
# Create a curve set, choosing the interval of distances [r_min, r_max]
curve_set <- crop_curves(env, r_min = 1, r_max = 8)
# For better visualisation, take the L(r)-r function
curve_set <- residual(curve_set, use_theo = TRUE)
# The studentised envelope test
res <- st_envelope(curve_set); plot(res, use_ggplot2=TRUE)

## Random labeling test
#----------------------
mpp <- spruces
# T(r) = \\hat{L}_m(r), an estimator of the L_m(r) function
curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=4999, r_min=0, r_max=9.5)
res <- st_envelope(curve_set)
plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))
}
\references{
  Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D.
  (2013). Deviation test construction and power comparison
  for marked spatial point patterns. arXiv:1306.1028

  Myllymäki, M., Mrkvička, T., Seijo, H. and Grabarnik, P.
  (2013). Global envelope tests for spatial point patterns.
}

