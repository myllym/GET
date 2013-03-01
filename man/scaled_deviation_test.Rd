\name{scaled_deviation_test}
\alias{scaled_deviation_test}
\title{Scale the curve_set and perform a deviation test.}
\usage{
  scaled_deviation_test(curve_set, r_min = NULL,
    r_max = NULL, use_theo = TRUE, scaling = "env",
    measure = "max")
}
\arguments{
  \item{curve_set}{A curve_set or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}

  \item{r_min}{The minimum radius to include.}

  \item{r_max}{The maximum radius to include.}
}
\description{
  Scale the curve_set and perform a deviation test.
}

