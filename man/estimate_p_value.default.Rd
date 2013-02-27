\name{estimate_p_value.default}
\alias{estimate_p_value.default}
\title{Estimate p-value.}
\usage{
  estimate_p_value.default(obs, sim_vec, ties = "random")
}
\arguments{
  \item{obs}{The data sample. A scalar real value. Must not
  be NULL.}

  \item{sim_vec}{The Monte Carlo samples. A vector of real
  values.  Must not be NULL.}

  \item{ties}{The method to treat tied values with. If one
  or more of the elements of sim_vec are equal to obs, how
  should the rank of obs be determined? For 'conservative'
  the resulting p-value will be the highest possible. For
  'liberal' the p-value will be the lowest possible. For
  'random' the rank of the obs within the tied values is
  uniformly sampled so that the resulting p-value is at
  most the conservative option and at least the liberal
  option. For 'midrank' we take the mid-rank within the
  tied values. 'random' is the default.}
}
\value{
  The p-value estimate. A scalar real value between 0 and
  1.
}
\description{
  Estimates the p-value of the given data sample for the
  given set of Monte Carlo samples. User can choose which
  method is used to treat tied values between the data
  sample and the Monte Carlo samples.
}
\references{
  Hajek, Sidak, Sen. Theory of Rank Tests. 1999. ff. 130.
}

