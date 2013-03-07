\name{normal_envelope}
\alias{normal_envelope}
\title{Approximative normal envelope test}
\usage{
  normal_envelope(curve_set, alpha = 0.05, n_norm = 2e+05,
    ...)
}
\arguments{
  \item{n_norm}{Number of simulations drawn from the
  multivariate normal distribution (dimension = number of
  distances r).}

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
  Approximative normal envelope test
}
\examples{
library(spatstat)
pp <- spruces
env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE)
curve_set <- spptest:::crop_curves(env, r_min = 0, r_max = 8)
curve_set <- spptest:::residual(curve_set, use_theo = TRUE)
res <- normal_envelope(curve_set, n_norm=20000)
plot(res)
}

