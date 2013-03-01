\name{scale_curves}
\alias{scale_curves}
\title{Scale summary functions.}
\usage{
  scale_curves(curve_set, scaling = "env", ...)
}
\arguments{
  \item{scaling}{The name of the scaling to use. Options
  include 'none', 'env', 'envdir' and 'scale'. 'env' is
  default.}

  \item{...}{Further arguments passed to the chosen scaling
  function.}

  \item{curve_set}{A curve_set or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}
}
\description{
  Scale summary functions.
}

