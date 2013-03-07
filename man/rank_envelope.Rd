\name{rank_envelope}
\alias{rank_envelope}
\title{The rank envelope test}
\usage{
  rank_envelope(curve_set, alpha = 0.05, ...)
}
\arguments{
  \item{alpha}{}

  \item{...}{Additional parameters passed to
  \code{\link{estimate_p_value}}.}

  \item{curve_set}{A curve_set or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}
}
\description{
  The rank envelope test
}
\examples{
library(spatstat)
pp <- spruces
env <- envelope(pp, fun="Lest", nsim=3999, savefuns=TRUE)
curve_set <- spptest:::crop_curves(env, r_min = 1, r_max = 7)
curve_set <- spptest:::residual(curve_set, use_theo = TRUE)
res <- rank_envelope(curve_set)
plot(res)
}

