\name{residual}
\alias{residual}
\title{Subtract S_{H_0} from the summary functions.}
\usage{
  residual(curve_set, use_theo = TRUE)
}
\arguments{
  \item{use_theo}{Whether to use the theoretical summary
  function or the mean of the simulations.}

  \item{curve_set}{A curve_set or an
  \code{\link[spatstat]{envelope}} object. If an envelope
  object is given, it must contain the summary functions
  from the simulated patterns which can be achieved by
  setting savefuns = TRUE when calling envelope().}
}
\value{
  A curve set object containing residual summary functions.
  theo is no longer included.
}
\description{
  Subtract S_{H_0} from the summary functions.
}

