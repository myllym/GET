\name{rank_envelope}
\alias{rank_envelope}
\title{The rank envelope test}
\usage{
  rank_envelope(curve_set, alpha = 0.05, savedevs = FALSE,
    alternative = "two.sided", lexo = FALSE, ...)
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

  \item{savedevs}{Logical. Should the global rank values
  k_i, i=1,...,nsim+1 be returned? Default: FALSE.}

  \item{alternative}{A character string specifying the
  alternative hypothesis. Must be one of the following:
  "two.sided" (default), "less" or "greater".}

  \item{lexo}{Logical, whether or not to use lexical
  ordering when ranking. See details.}

  \item{...}{Additional parameters passed to
  \code{\link{estimate_p_value}} to obtain a point estimate
  for the p-value. The default point estimate is the
  mid-rank p-value.}
}
\value{
  An "envelope_test" object containing the following
  fields: \itemize{ \item r = Distances for which the test
  was made.  \item method = The name of the envelope test.
  \item alternative = The alternative specified in the
  function call.  \item p = A point estimate for the
  p-value (default is the mid-rank p-value).  \item lexo =
  As the argument \code{lexo}.  \item p_interval = The
  p-value interval [p_liberal, p_conservative].  \item
  k_alpha = The value of k corresponding to the
  100(1-alpha)\% simultaneous envelope.  \item k = Global
  rank values (k[1] is the value for the data pattern).
  Returned only if savedevs = TRUE.  \item central_curve =
  If the curve_set (or envelope object) contains a
  component 'theo', then this function is used as the
  central curve and returned in this component.  Otherwise,
  the central_curve is the mean of the test functions
  T_i(r), i=2, ..., s+1.  Used for visualization only.
  \item data_curve = The test function for the data.  \item
  lower = The lower envelope.  \item upper = The upper
  envelope.  \item call = The call of the function. }
}
\description{
  The rank envelope test, p-value and simultaneous envelope
}
\details{
  The rank envelope test is a completely non-parametric
  test, which provides a p-value interval given by the most
  liberal and the most conservative p-value estimate and
  the simultaneous 100(1-alpha)\% envelope for the chosen
  test function T(r) on the chosen interval of distances.

  Given a curve_set object, the test is carried out as
  follows.

  For each curve in the curve_set, both the data curve and
  the simulations, the global rank measure k is determined.
  If savedevs = TRUE, then the global rank values k_1, k_2,
  ..., k_(s+1) are returned in the component 'k', where
  k[1] is the value for the data.

  Based on k_i, i=1, ..., s+1, the p-interval is
  calculated. This interval is by default plotted for the
  object returned by the rank_envelope function. Also a
  single p-value is calculated and returned in component
  'p'. By default this p-value is the mid-rank p-value, but
  another option can be used by passing additional
  parameters to \code{\link{estimate_p_value}}.

  The simultaneous 100(1-alpha)\% envelope is given by the
  'k_alpha'th lower and upper envelope. For details see
  Myllymäki et al. (2013).

  The above holds if \code{lexo == FALSE} and then the test
  corresponds to the rank envelope test by Myllymaki et. al
  (2013). If \code{lexo == TRUE}, then all the pointwise
  ranks are used to rank the curves, by so called lexical
  ordering. The test is no longer usable as a graphical
  test. If \code{lexo} was set to a value \code{TRUE}, a
  warning is issued when plotting using
  \code{\link{plot.envelope_test}}).
}
\examples{
## Testing complete spatial randomness (CSR)
#-------------------------------------------
require(spatstat)
pp <- unmark(spruces)
# Generate nsim simulations under CSR, calculate L-function for the data and simulations
env <- envelope(pp, fun="Lest", nsim=4999, savefuns=TRUE, correction="translate")
# The rank envelope test
res <- rank_envelope(env)
# Plot the result.
# - The central curve is now obtained from env[['theo']], which is the
# value of the L-function under the null hypothesis (L(r) = r).
plot(res)
# or (requires R library ggplot2)
plot(res, use_ggplot2=TRUE)

## Advanced use:
# Create a curve set, choosing the interval of distances [r_min, r_max]
curve_set <- crop_curves(env, r_min = 1, r_max = 7)
# For better visualisation, take the L(r)-r function
curve_set <- residual(curve_set, use_theo = TRUE)
# Do the rank envelope test
res <- rank_envelope(curve_set); plot(res, use_ggplot2=TRUE)

## Random labeling test
#----------------------
# requires library 'marksummary'
mpp <- spruces
# 1) Perform simulations under the random labelling hypothesis and calculate
# the test function T(r) for the data pattern (mpp) and each simulation.
# The command below specifies that the test function is T(r) = \\hat{L}_m(r),
# which is an estimator of the mark-weighted L function, L_m(r),
# with translational edge correction (default).
# The random_labelling function returns the centred functions \\hat{L}_m(r)-T_0(r),
# where T_0(r) = \\hat{L}(r) is the unmarked L function.
curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=4999, r_min=1.5, r_max=9.5)
# 2) Do the rank envelope test
res <- rank_envelope(curve_set)
# 3) Plot the test result
plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))

# Make the test using instead the test function T(r) = \\hat{L}_mm(r);
# which is an estimator of the mark-weighted L function, L_mm(r),
# with translational edge correction (default).
curve_set <- random_labelling(mpp, mtf_name = 'mm', nsim=4999, r_min=1.5, r_max=9.5)
res <- rank_envelope(curve_set)
plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[mm](r)-L(r))))

## Goodness-of-fit test
#----------------------
pp <- unmark(spruces)
# Minimum distance between points in the pattern
min(nndist(pp))
# Fit a model
fittedmodel <- ppm(pp, interaction=Hardcore(hc=1)) # Hardcore process

\dontrun{
# Simulating Gibbs process by 'envelope' is slow, because it uses the MCMC algorithm
#env <- envelope(fittedmodel, fun="Jest", nsim=999, savefuns=TRUE, correction="none", r=seq(0, 4, length=500))

# Using direct algorihm can be faster, because the perfect simulation is used here.
simulations <- NULL
for(j in 1:4999) {
   simulations[[j]] <- rHardcore(beta=exp(fittedmodel$coef[1]), R = fittedmodel$interaction$par$hc, W = pp$window);
   if(j\%\%10==0) cat(j, "...", sep="")
}
env <- envelope(pp, simulate=simulations, fun="Jest", nsim=length(simulations), savefuns=TRUE, correction="none", r=seq(0, 4, length=500))
curve_set <- crop_curves(env, r_min = 1, r_max = 3.5)
res <- rank_envelope(curve_set); plot(res, use_ggplot2=TRUE)
}
}
\references{
  Myllymäki, M., Mrkvička, T., Seijo, H., Grabarnik, P.
  (2013). Global envelope tests for spatial point patterns.
  arXiv:1307.0239 [stat.ME]
}
\seealso{
  \code{\link{random_labelling}}
}

