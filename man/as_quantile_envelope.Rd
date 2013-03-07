\name{as_quantile_envelope}
\alias{as_quantile_envelope}
\title{Asymmetric quantile envelope test}
\usage{
  as_quantile_envelope(curve_set, alpha = 0.05, ...)
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
  Asymmetric quantile envelope test
}
\examples{
library(spatstat)
pp <- spruces
env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE)
curve_set <- crop_curves(env, r_min = 0, r_max = 8)
curve_set <- residual(curve_set, use_theo = TRUE)
res <- as_quantile_envelope(curve_set)
plot(res)
}

