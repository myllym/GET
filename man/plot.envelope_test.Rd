\name{plot.envelope_test}
\alias{plot.envelope_test}
\title{Plot method for the class 'envelope_test'}
\usage{
  \method{plot}{envelope_test}(x, use_ggplot2=FALSE, main,
    ylim, xlab, ylab, ...)
}
\arguments{
  \item{x}{an 'envelope_test' object}

  \item{use_ggplot2}{TRUE/FALSE, If TRUE, then a plot with
  a coloured envelope ribbon is provided. Requires R
  library ggplot2.}

  \item{main}{See \code{\link{plot.default}}. A sensible
  default exists.}

  \item{ylim}{See \code{\link{plot.default}}. A sensible
  default exists.}

  \item{xlab}{See \code{\link{plot.default}}. A sensible
  default exists.}

  \item{ylab}{See \code{\link{plot.default}}. A sensible
  default exists.}

  \item{...}{Additional parameters to be passed to plot (if
  use_ggplot2=FALSE).}
}
\description{
  Plot method for the class 'envelope_test'
}
\seealso{
  \code{\link{rank_envelope}}, \code{\link{st_envelope}},
  \code{\link{qdir_envelope}}
}

