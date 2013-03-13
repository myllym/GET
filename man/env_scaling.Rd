\name{env_scaling}
\alias{env_scaling}
\title{Env scaling.}
\usage{
  env_scaling(curve_set, probs = c(0.025, 0.975), ...)
}
\arguments{
  \item{probs}{A two-element vector containing the lower
  and upper quantiles for the envelope, in that order and
  on the interval [0, 1].  The default values are 0.025 and
  0.975 as in the article by Møller and Berthelsen.}

  \item{...}{Further arguments passed to quantile.}

  \item{curve_set}{A curve_set or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}
}
\value{
  A scaled curve_set.
}
\description{
  Env scaling.
}
\references{
  J. Møller and K. K. Berthelsen, “Transforming spatial
  point processes into Poisson processes using random
  superposition,” Advances in Applied Probability, vol. 44,
  no. 1, pp. 42–62, 2012.
}

