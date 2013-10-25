\name{get_T_0}
\alias{get_T_0}
\title{Define T_0 from a curve_set object}
\usage{
  get_T_0(curve_set)
}
\arguments{
  \item{curve_set}{A curve_set or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}
}
\description{
  Define T_0, the expectation of the test function under
  H_0, from a curve_set object.
}

