\name{sd_envelope}
\alias{sd_envelope}
\title{Variance stabilized envelope test}
\usage{
  sd_envelope(curve_set, alpha = 0.05, ...)
}
\arguments{
  \item{curve_set}{A curve_set or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}

  \item{alpha}{}

  \item{...}{Additional parameters passed to
  \code{\link{estimate_p_value}}.}
}
\description{
  Variance stabilized envelope test
}
\examples{
library(spatstat)
pp <- spruces
env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE)
curve_set <- spptest:::crop_curves(env, r_min = 0, r_max = 8)
curve_set <- spptest:::residual(curve_set, use_theo = TRUE)
res <- sd_envelope(curve_set)
plot(res)
plot(res, xlab=expression(italic(r)), ylab=expression(italic(L(r)-r)))
}

