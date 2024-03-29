% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/appl_spatialF.r
\name{GET.spatialF}
\alias{GET.spatialF}
\title{Testing global and local dependence of point patterns on covariates}
\usage{
GET.spatialF(
  X,
  formula.full,
  formula.reduced,
  fitfun,
  covariates,
  nsim,
  bw = spatstat.explore::bw.scott(X),
  bw.S = bw,
  dimyx = NULL,
  ...
)
}
\arguments{
\item{X}{A \code{ppp} object of \pkg{spatstat} representing the observed point pattern.}

\item{formula.full}{A formula for the trend of the full model.}

\item{formula.reduced}{A formula for the trend of the reduced model
that is a submodel of the full model.}

\item{fitfun}{A function of a point pattern, model formula and covariates,
giving a fitted model object that can be used with \code{\link[stats]{simulate}}.}

\item{covariates}{A list of covariates.}

\item{nsim}{The number of simulations.}

\item{bw}{The bandwidth for smoothed residuals.}

\item{bw.S}{The radius for the local S(u)-statistic.}

\item{dimyx}{Pixel array dimensions for smoothed residuals. See \code{as.mask} of \pkg{spatstat}.}

\item{...}{Additional arguments to be passed to \code{\link{global_envelope_test}}.}
}
\value{
list with three components
\itemize{
\item F = the global envelope test based on the F(u) statistic
\item S = the global envelope test based on the S(u) statistic
\item coef = the coefficients of the full model given by fitfun
}
}
\description{
Compute the spatial F- and S-statistics and perform the one-stage global envelope tests
proposed by Myllymäki et al. (2020).
}
\examples{
if(require("spatstat.model", quietly=TRUE)) {
  # Example of tropical rain forest trees
  data("bei")

  fullmodel <- ~ grad
  reducedmodel <- ~ 1
  fitppm <- function(X, model, covariates) {
    ppm(X, model, covariates=covariates)
  }
  \dontshow{res <- GET.spatialF(bei, fullmodel, reducedmodel, fitppm, bei.extra, 3, alpha=0.5, dimyx=32)}
  \donttest{
  nsim <- 19 # Increase nsim for serious analysis!
  res <- GET.spatialF(bei, fullmodel, reducedmodel, fitppm, bei.extra, nsim)
  }
  plot(res$F)
  plot(res$S)

  \donttest{
  # Example of forest fires
  data("clmfires")
  # Choose the locations of the lightnings in years 2004-2007:
  pp.lightning <- unmark(subset(clmfires, cause == "lightning" &
                   date >= "2004-01-01" & date < "2008-01-01"))

  covariates <- clmfires.extra$clmcov100
  covariates$forest <- covariates$landuse == "conifer" | covariates$landuse == "denseforest" |
                        covariates$landuse == "mixedforest"

  fullmodel <- ~ elevation + landuse
  reducedmodel <- ~ landuse
  nsim <- 19 # Increase nsim for serious analysis!
  res <- GET.spatialF(pp.lightning, fullmodel, reducedmodel, fitppm, covariates, nsim)
  plot(res$F)
  plot(res$S)

  # Examples of the fitfun functions for clustered and regular processes
  # fitfun for the log Gaussian Cox Process with exponential covariance function
  fitLGCPexp <- function(X, model, covariates) {
    kppm(X, model, clusters="LGCP", model="exponential", covariates=covariates)
  }
  # fitfun for the hardcore process with hardcore radius 0.01
  fitHardcore <- function(X, model, covariates) {
    ppm(X, model, interaction=Hardcore(0.01), covariates=covariates)
  }
  }
}
}
\references{
Myllymäki, M., Kuronen, M. and Mrkvička, T. (2020). Testing global and local dependence of point patterns on covariates in parametric models. Spatial Statistics 42, 100436. doi: 10.1016/j.spasta.2020.100436
}
