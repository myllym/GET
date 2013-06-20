\name{scale_curves}
\alias{scale_curves}
\title{Scale curves.}
\usage{
  scale_curves(curve_set, scaling = "qdir", ...)
}
\arguments{
  \item{scaling}{The name of the scaling to use. Options
  include 'none', 'q', 'qdir' and 'st'. 'qdir' is default.}

  \item{...}{Further arguments passed to the chosen scaling
  function.}

  \item{curve_set}{A curve_set or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}
}
\description{
  The most important use is: scale residuals of test
  functions.
}
\details{
  Given a set of curves in curve_set, the function
  scale_curves scales the curves by one of the following
  scalings: \itemize{ \item none No scaling. Nothing done.
  \item q Quantile scaling.  \item qdir Directional
  quantile scaling.  \item st Studentised scaling. } See
  for details Myllymäki et al. (2013).
}
\references{
  Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D.
  (2013). Deviation test construction and power comparison
  for marked spatial point patterns. arXiv:1306.1028
}

