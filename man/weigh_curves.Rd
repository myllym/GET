\name{weigh_curves}
\alias{weigh_curves}
\title{Multiply by a coefficient.}
\usage{
  weigh_curves(curve_set, coeff)
}
\arguments{
  \item{curve_set}{The curve_set to be weighted. Must be of
  class curve_set.}

  \item{coeff}{The coefficient vector, often of the length
  of one curve.}

  \item{zero_idx}{Something that can be used to index the
  curve vectors.  These values of the curves will be set to
  zero. If NULL, nothing will be zeroed unless coeff has
  zeros. Defaults to NULL.}
}
\description{
  Multiply by a coefficient.
}

