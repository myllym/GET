\name{convert_envelope}
\alias{convert_envelope}
\title{Convert an envelope object to a curve_set object.}
\usage{
  convert_envelope(curve_set)
}
\arguments{
  \item{curve_set}{A curve_set or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}
}
\value{
  If an \code{\link[spatstat]{envelope}} object was given,
  return a corresponding curve_set object. If a curve_set
  object was given, return it unharmed.
}
\description{
  If given an envelope object, convert it into a curve_set
  object. If given a curve_set object, check its
  correctness and give it back.
}

