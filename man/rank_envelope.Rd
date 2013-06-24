\name{rank_envelope}
\alias{rank_envelope}
\title{The rank envelope test}
\usage{
  rank_envelope(curve_set, alpha = 0.05, savedevs = FALSE,
    ...)
}
\arguments{
  \item{curve_set}{A curve_set (see
  \code{\link{create_curve_set}}) or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}

  \item{alpha}{The significance level. Simultaneous
  100(1-alpha) percent envelopes will be calculated.}

  \item{savedevs}{Logical. Should the global rank values
  k_i, i=1,...,nsim+1 be returned? Default: FALSE.}

  \item{...}{Additional parameters passed to
  \code{\link{estimate_p_value}} to obtain a point estimate
  for the p-value.}
}
\value{
  An "envelope_test" object containing the following
  fields: \itemize{ \item r Distances for which the test
  was made.  \item method The envelope method.  \item p A
  point estimate for the p-value.  \item p_interval The
  p-value interval [p_liberal, p_conservative].  \item
  k_alpha \item k Global rank values. Returned only if
  savedevs = TRUE.  \item central_curve The mean test
  function (median) calculated from simulations.  \item
  data_curve The test function for the data.  \item lower
  The lower envelope.  \item upper The upper envelope.
  \item call The call of the function. }
}
\description{
  The rank envelope test, p-value and simultaneous
  envelopes
}
\details{
  The rank envelope test is a completely non-parametric
  test, which provides a p-value interval given by the most
  liberal and the most conservative p-value estimate and
  the simultaneous 100(1-alpha) percent envelopes for the
  chosen test function T(r) on the chosen interval of
  distances.
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
plot(res)
# or (requires R library ggplot2)
plot(res, use_ggplot2=TRUE)

## Advanced use:
# Create a curve set, choosing the interval of distances [r_min, r_max]
curve_set <- crop_curves(env, r_min = 1, r_max = 7)
# For better visualisation, take the L(r)-r function
curve_set <- residual(curve_set, use_theo = TRUE)
# The rank envelope test
res <- rank_envelope(curve_set); plot(res, use_ggplot2=TRUE)

## Random labeling test
#----------------------
mpp <- spruces
# T(r) = \\hat{L}_m(r), an estimator of the L_m(r) function
curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=4999, r_min=0, r_max=9.5)
res <- rank_envelope(curve_set)
plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))
# T(r) = \\hat{L}_mm(r), an estimator of the L_mm(r) function
curve_set <- random_labelling(mpp, mtf_name = 'mm', nsim=4999, r_min=0, r_max=9.5)
res <- rank_envelope(curve_set)
plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[mm](r)-L(r))))

## Goodness-of-fit test
#----------------------
# FIXME What is a reasonable model to be fitted to this data?
pp <- unmark(spruces)
# Minimum distance between points in the pattern
min(nndist(pp))
# Fit a model
## fittedmodel <- ppm(pp, interaction=Strauss(r=4)) # Strauss process
fittedmodel <- ppm(pp, interaction=Hardcore(hc=1)) # Hardcore process

\dontrun{
## Simulating Strauss process by 'envelope' is slow
#env <- envelope(fittedmodel, fun="Lest", nsim=999, savefuns=TRUE, correction="none", r=seq(0, 9.5, length=500))

simulations <- NULL
for(j in 1:999) {
   ##simulations[[j]] <- rStrauss(beta=exp(fittedmodel$coef[1]), gamma=exp(fittedmodel$coef[2]), R=fittedmodel$interaction$par$r, W=pp$window);
   simulations[[j]] <- rHardcore(beta=exp(fittedmodel$coef[1]), R = fittedmodel$interaction$par$hc, W = pp$window);
   if(j\%\%10==0) cat(j, "...", sep="")
}
env <- envelope(pp, simulate=simulations, fun="Lest", nsim=length(simulations), savefuns=TRUE, correction="none", r=seq(0, 9.5, length=500))
res <- rank_envelope(env); plot(res, use_ggplot2=TRUE)
}
}
\references{
  Myllymäki, M., Mrkvička, T., Seijo, H., Grabarnik, P.
  (2013). Global envelope tests for spatial point patterns.
}
\seealso{
  \code{\link{random_labelling}}
}

