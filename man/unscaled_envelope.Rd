\name{unscaled_envelope}
\alias{unscaled_envelope}
\title{Unscaled envelope test}
\usage{
  unscaled_envelope(curve_set, alpha = 0.05,
    savedevs = FALSE, ...)
}
\arguments{
  \item{curve_set}{A curve_set (see
  \code{\link{create_curve_set}}) or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}

  \item{alpha}{The significance level. Simultaneous
  100(1-alpha)\% envelopes will be calculated.}

  \item{savedevs}{Logical. Should the deviation values u_i,
  i=1,...,nsim+1 be returned? Default: FALSE.}

  \item{...}{Additional parameters passed to
  \code{\link{estimate_p_value}} to obtain a point estimate
  for the p-value. The default point estimate is the
  mid-rank p-value. The choice should not affect the
  result, since no ties are expected to occur.}
}
\value{
  An "envelope_test" object containing the following
  fields: \itemize{ \item r = Distances for which the test
  was made.  \item method = The name of the envelope test.
  \item p = A point estimate for the p-value (default is
  the mid-rank p-value).  \item u_alpha = The value of u
  corresponding to the 100(1-alpha)\% simultaneous
  envelope.  \item u = Deviation values (u[1] is the value
  for the data pattern). Returned only if savedevs = TRUE.
  \item central_curve = If the curve_set (or envelope
  object) contains a component 'theo', then this function
  is used as the central curve and returned in this
  component.  Otherwise, the central_curve is the mean of
  the test functions T_i(r), i=2, ..., s+1.  \item
  data_curve = The test function for the data.  \item lower
  = The lower envelope.  \item upper = The upper envelope.
  \item call = The call of the function. }
}
\description{
  The unscaled envelope test, which leads to envelopes with
  constant width over the distances r. It corresponds to
  the classical maximum deviation test without scaling.
}
\details{
  This test suffers from unequal variance of T(r) over the
  distances r and from the asymmetry of distribution of
  T(r). We recommend to use the rank_envelope (if number of
  simulations close to 5000 can be afforded) or
  st_envelope/qdir_envelope (if large number of simulations
  cannot be afforded) instead.
}
\examples{
## Testing complete spatial randomness (CSR)
#-------------------------------------------
require(spatstat)
pp <- spruces
## Test for complete spatial randomness (CSR)
# Generate nsim simulations under CSR, calculate L-function for the data and simulations
env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
# The studentised envelope test
res <- unscaled_envelope(env)
plot(res)
# or (requires R library ggplot2)
plot(res, use_ggplot2=TRUE)

## Advanced use:
# Create a curve set, choosing the interval of distances [r_min, r_max]
curve_set <- crop_curves(env, r_min = 1, r_max = 8)
# For better visualisation, take the L(r)-r function
curve_set <- residual(curve_set, use_theo = TRUE)
# The studentised envelope test
res <- unscaled_envelope(curve_set); plot(res, use_ggplot2=TRUE)

## Random labeling test
#----------------------
# requires library 'marksummary'
mpp <- spruces
# Use the test function T(r) = \\hat{L}_m(r), an estimator of the L_m(r) function
curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=4999, r_min=1.5, r_max=9.5)
res <- unscaled_envelope(curve_set)
plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))
}
\references{
  Ripley, B.D. (1981). Spatial statistics. Wiley, New
  Jersey.
}

