# Calculates the mean in the circle of size r, b(i,r), at each pixel of the image.
# At the edge, the mean is the mean of the image values in the image window.
# @param x An image.
# @param R The radius of the circle.
maverage2d <- function(x, R) {
  R <- mean(R) # Note: R can be a vector of length 2 corresponding to bandwidths in x and y directions
  X <- spatstat::fillNA(x, 0)
  # Calculates the sum of values of X in b(i,R) at each pixel i
  Y <- spatstat::second.moment.calc(X, what="smooth", kernel=function(x, y) x^2 + y^2 < R^2, sigma=0)
  # Calculate the number of pixels in each b(i,R)
  X[!is.na(x$v)] <- 1
  Yden <- spatstat::second.moment.calc(X, what="smooth", kernel=function(x, y) x^2 + y^2 < R^2, sigma=0)
  # Calculate the mean
  Z <- spatstat::eval.im(Y/Yden)
  Z[is.na(x$v)] <- NA
  Z
}

# Compute the spatial statistics F(u) and S(u), and the coefficients of the full model
# @param returncoef Logical. Whether or not to calculate and return also the coefficients
# of the full model.
compute_statistics <- function(X, formula.full, formula.reduced, fitfun, covariates, bw, bw.S, returncoef=FALSE, dimyx=NULL) {
  M.full <- fitfun(X, formula.full, covariates)
  M.reduced <- fitfun(X, formula.reduced, covariates)
  r0 <- residuals(M.reduced)
  r1 <- residuals(M.full)
  s0 <- spatstat::Smooth(r0, sigma=bw, dimyx=dimyx)
  s1 <- spatstat::Smooth(r1, sigma=bw, dimyx=dimyx)
  sm0 <- maverage2d(s0, bw.S)
  sq0 <- maverage2d(s0^2, bw.S)
  sm1 <- maverage2d(s1, bw.S)
  sq1 <- maverage2d(s1^2, bw.S)
  Fstat <- list(F=s0^2 - s1^2, S = (sq0-sm0^2)/(sq1-sm1^2))
  ret <- list(Fstat=Fstat)
  if(returncoef) ret$coef <- coef(summary(M.full))
  ret
}

# obs and sim_m are matrices that can include NAs, they are removed
curve_set_helper <- function(r, obs, sim_m) {
  ok <- !is.na(obs)
  xy <- expand.grid(x=r$x, y=r$y, KEEP.OUT.ATTRS = FALSE)[ok,]
  xy$width <- r$x[2] - r$x[1]
  xy$height <- r$y[2] - r$y[1]
  create_curve_set(list(r=xy, obs=obs[ok], sim_m = matrix(sim_m[ok], nrow=sum(ok))))
}

#' Testing global and local dependence of point patterns on covariates
#' 
#' Compute the spatial F- and S-statistics and perform the one-stage global envelope tests
#' proposed by Myllymäki et al. (2020).
#' 
#' @param X A \code{ppp} object of \pkg{spatstat} representing the observed point pattern.
#' @param formula.full A formula for the trend of the full model.
#' @param formula.reduced A formula for the trend of the reduced model
#' that is a submodel of the full model.
#' @param fitfun A function of a point pattern, model formula and covariates,
#'   giving a fitted model object that can be used with \code{\link[stats]{simulate}}.
#' @param covariates A list of covariates.
#' @param nsim The number of simulations.
#' @param bw The bandwidth for smoothed residuals.
#' @param bw.S The radius for the local S(u)-statistic.
#' @param dimyx Pixel array dimensions for smoothed residuals. See \code{as.mask} of \pkg{spatstat}.
#' @param ... Additional arguments to be passed to \code{\link{global_envelope_test}}.
#' @return list with three components
#' \itemize{
#' \item F = the global envelope test based on the F(u) statistic
#' \item S = the global envelope test based on the S(u) statistic
#' \item coef = the coefficients of the full model given by fitfun
#' }
#' @export
#' @references
#' Myllymäki, M., Kuronen, M. and Mrkvička, T. (2020). Testing global and local dependence of point patterns on covariates in parametric models. Spatial Statistics. doi: 10.1016/j.spasta.2020.100436
#' @importFrom stats simulate
#' @examples
#' if(require("spatstat", quietly=TRUE)) {
#'   # Example of tropical rain forest trees
#'   data("bei")
#'
#'   fullmodel <- ~ grad
#'   reducedmodel <- ~ 1
#'   fitppm <- function(X, model, covariates) {
#'     ppm(X, model, covariates=covariates)
#'   }
#'   \dontshow{res <- GET.spatialF(bei, fullmodel, reducedmodel, fitppm, bei.extra, 3, alpha=0.5, dimyx=32)}
#'   \donttest{
#'   nsim <- 19 # Increase nsim for serious analysis!
#'   res <- GET.spatialF(bei, fullmodel, reducedmodel, fitppm, bei.extra, nsim)
#'   }
#'   plot(res$F)
#'   plot(res$S)
#'
#'   \donttest{
#'   # Example of forest fires
#'   data("clmfires")
#'   # Choose the locations of the lightnings in years 2004-2007:
#'   pp.lightning <- unmark(subset(clmfires, cause == "lightning" &
#'                    date >= "2004-01-01" & date < "2008-01-01"))
#'
#'   covariates <- clmfires.extra$clmcov100
#'   covariates$forest <- covariates$landuse == "conifer" | covariates$landuse == "denseforest" |
#'                         covariates$landuse == "mixedforest"
#'
#'   fullmodel <- ~ elevation + landuse
#'   reducedmodel <- ~ landuse
#'   nsim <- 19 # Increase nsim for serious analysis!
#'   res <- GET.spatialF(pp.lightning, fullmodel, reducedmodel, fitppm, covariates, nsim)
#'   plot(res$F)
#'   plot(res$S)
#'
#'   # Examples of the fitfun functions for clustered and regular processes
#'   # fitfun for the log Gaussian Cox Process with exponential covariance function
#'   fitLGCPexp <- function(X, model, covariates) {
#'     kppm(X, model, clusters="LGCP", model="exponential", covariates=covariates)
#'   }
#'   # fitfun for the hardcore process with hardcore radius 0.01
#'   fitHardcore <- function(X, model, covariates) {
#'     ppm(X, model, interaction=Hardcore(0.01), covariates = covariates)
#'   }
#'   }
#' }
GET.spatialF <- function(X, formula.full, formula.reduced, fitfun, covariates, nsim,
                         bw=spatstat::bw.scott(X), bw.S=bw, dimyx=NULL, ...) {
  if(!spatstat::is.ppp(X)) stop("X should be a ppp object.")
  check_isnested(formula.full, formula.reduced)
  M.reduced <- fitfun(X, formula.reduced, covariates)
  sim <- simulate(M.reduced, nsim = nsim)
  fun <- function(data, returncoef=FALSE) {
    compute_statistics(data, formula.full, formula.reduced, fitfun, covariates, bw, bw.S, returncoef=returncoef, dimyx=dimyx)
  }
  obs <- fun(X, TRUE)
  r <- list(x=obs$Fstat[[1]]$xcol, y=obs$Fstat[[1]]$yrow)
  sim.funcs <- lapply(sim, fun)
  fun2 <- function(name) {
    obs.func <- t(as.matrix(obs$Fstat[[name]]))
    sim_m <- simplify2array(lapply(sim.funcs, function(x) t(as.matrix(x$Fstat[[name]]))))
    cset <- curve_set_helper(r, obs.func, sim_m)
    global_envelope_test(cset, alternative="greater", ...)
  }
  res <- lapply(names(obs$Fstat), fun2)
  names(res) <- names(obs$Fstat)
  res$coef <- obs$coef
  res
}
