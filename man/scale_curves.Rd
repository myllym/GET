\name{scale_curves}
\alias{scale_curves}
\title{Scale summary functions.}
\usage{
  scale_curves(curve_set, scaling = "env", ...)
}
\arguments{
  \item{curve_set}{The estimated summary functions. Either
  a curve_set object or an \code{\link[spatstat]{envelope}}
  object.}

  \item{scaling}{The name of the scaling to use. Options
  include 'env', 'envdir', 'scale'. 'env' is default.}

  \item{...}{Further arguments passed to the chosen scaling
  function.}
}
\description{
  Scale summary functions.
}

