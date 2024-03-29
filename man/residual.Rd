% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residual.r
\name{residual}
\alias{residual}
\title{Residual form of the functions}
\usage{
residual(curve_set, use_theo = TRUE)
}
\arguments{
\item{curve_set}{A \code{\link{curve_set}} object, or
an \code{envelope} object of \pkg{spatstat}. If an envelope object is given,
it must contain the summary functions from the simulated patterns which can be
achieved by setting savefuns = TRUE when calling the \code{envelope} function.}

\item{use_theo}{Whether to use the theoretical summary function or the
mean of the functions in the curve_set.}
}
\value{
A curve set object containing residual summary functions. theo is
  no longer included.
}
\description{
Subtract the theoretical function \eqn{S_{H_0}}{S_{H_0}} or the mean of the functions
in the curve set. If the \code{curve_set} object contains already residuals
\eqn{T_i(r) - T_0(r)}{T_i(r) - T_0(r)}, \code{use_theo} ignored and the same object returned.
}
\details{
The mean of the functions in the \code{curve_set} is the mean of all functions.
If \code{use_theo = TRUE}, but the component \code{theo} does not exist in the
\code{curve_set}, the mean of the functions is used silently.
}
