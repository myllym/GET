\name{deviation}
\alias{deviation}
\title{Deviation of curves.}
\usage{
  deviation(curve_set, measure = "max", ...)
}
\arguments{
  \item{curve_set}{A residual curve_set object.}

  \item{measure}{The deviation measure to use. Default is
  'max'. Must be one of the following: 'max', 'int',
  'int2', 'MBD', 'MHRD', 'MMHRD', 'ISD' or 'MISD'.}

  \item{...}{Further arguments delivered to the measure
  function.}
}
\value{
  A deviation_set object. The list has two elements: obs
  and sim.  obs is scalar while sim is a vector with at
  least one element.
}
\description{
  Deviation of curves.
}

