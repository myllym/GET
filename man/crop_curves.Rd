\name{crop_curves}
\alias{crop_curves}
\title{Crop the curves to a certain interval.}
\usage{
  crop_curves(curve_set, r_min = NULL, r_max = NULL)
}
\arguments{
  \item{curve_set}{The collection of curves. Must either
  have class \code{\link[spatstat]{envelope}} or class
  curve_set.}
}
\description{
  Crop the curves to a certain interval in preparation for
  a deviation test.
}

