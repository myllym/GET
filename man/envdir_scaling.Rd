\name{envdir_scaling}
\alias{envdir_scaling}
\title{EnvDir scaling.}
\usage{
  envdir_scaling(curve_set, probs = c(0.025, 0.975),
    type = 8, ...)
}
\arguments{
  \item{curve_set}{The curves to be scaled.}

  \item{probs}{The quantile percentiles to use. Should be
  within [0, 1].  See \code{\link[stats]{quantile}} for
  further details.}

  \item{type}{The quantile algorithm to use. See
  \code{\link[stats]{quantile}} for further details.}

  \item{...}{Further arguments passed to quantile.}
}
\description{
  EnvDir scaling.
}

