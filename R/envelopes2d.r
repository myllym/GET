# Functionality for central_region_2d and global_envelope_test_2d
cr_or_GET_2d <- function(image_set, CR_or_GET = c("CR", "GET"), ...) {
  CR_or_GET <- match.arg(CR_or_GET)
  obs_d <- dim(obs)
  sim_d <- dim(sim)
  image_set <- check_image_set_dimensions(image_set)
  # Create curve_set transforming the 2d functions to vectors
  curve_set_v <- image_set_to_curve_set(image_set)

  # Prepare global envelope or global envelope test
  switch(CR_or_GET,
         CR = {
           res_v <- central_region(curve_set_v, ...)
         },
         GET = {
           res_v <- global_envelope_test(curve_set_v, ...)
         })

  # Transform the results to 2d
  if(is.vector(curve_set_v[['obs']]))
    res <- structure(list(rx=image_set$r[[1]], ry=image_set$r[[2]],
                          obs=obs, #matrix(res_v$obs, nrow=obs_d[1], ncol=obs_d[2]),
                          central=matrix(res_v$central, nrow=obs_d[1], ncol=obs_d[2]),
                          lo=matrix(res_v$lo, nrow=obs_d[1], ncol=obs_d[2]),
                          hi=matrix(res_v$hi, nrow=obs_d[1], ncol=obs_d[2])),
                     class = c("global_envelope_2d", "list"))
  else
    res <- structure(list(rx=image_set$r[[1]], ry=image_set$r[[2]],
                          central=matrix(res_v$central, nrow=obs_d[2], ncol=obs_d[3]),
                          lo=matrix(res_v$lo, nrow=obs_d[2], ncol=obs_d[3]),
                          hi=matrix(res_v$hi, nrow=obs_d[2], ncol=obs_d[3])),
                     class = c("global_envelope_2d", "list"))
  attr(res, "method") <- attr(res_v, "method")
  attr(res, "type") <- attr(res_v, "type")
  attr(res, "alternative") <- attr(res_v, "einfo")$alternative
  attr(res, "k_alpha") <- attr(res_v, "k_alpha")
  attr(res, "alpha") <- attr(res_v, "alpha")
  attr(res, "k") <- attr(res_v, "k")
  if(CR_or_GET == "GET") {
    attr(res, "p") <- attr(res_v, "p")
    if(attr(res_v, "type") == "rank") {
      attr(res, "p_interval") <- attr(res_v, "p_interval")
      attr(res, "ties") <- attr(res_v, "ties")
    }
  }
  attr(res, "call") <- match.call()
  res
}

#' 2D central region / global envelope
#'
#' Provides central regions or global envelopes or confidence bands in 2D
#'
#' @references
#' Mrkvička, T., Soubeyrand, S., Myllymäki, M., Grabarnik, P. and Hahn, U. (2016). Monte Carlo testing in spatial statistics, with applications to spatial residuals. Spatial Statistics 18, Part A: 40--53. doi: 10.1016/j.spasta.2016.04.005
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' Mrkvička, T., Hahn, U. and Myllymäki, M. (2018). A one-way ANOVA test for functional data with graphical interpretation. arXiv:1612.03608 [stat.ME]
#'
#' @param image_set An image set, i.e. a set of 2d functions. See \code{\link{create_image_set}}.
#' @param ... Additional parameters to be passed to \code{\link{central_region}}.
#' @return An object of class "global_envelope_2d" (and "list"),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item rx = the vector of values of the argument rx on x-coordinate at which the test was made
#' \item rx = the vector of values of the argument ry on y-coordinate at which the test was made
#' \item obs = the data function (matrix), if there is only one data function. Otherwise not existing.
#' \item lo = the lower envelope (matrix) based on the simulated functions
#' \item hi = the upper envelope (matrix) based on the simulated functions
#' \item central = the central curve, mean or median (default) of the test functions
#' T_i(r), i=2, ..., s+1. Used for visualization only.
#' }
#' Additionally, the return value has attributes \code{method}, \code{type}, \code{alternative},
#' \code{k_alpha}, \code{alpha}, \code{k}, and \code{call}, see more detailed description in
#' \code{\link{central_region}}.
#' @export
central_region_2d <- function(image_set, ...) {
  if(class(image_set)[1] != "image_set") stop("image_set should be an object of class image_set.\n")
  cr_or_GET_2d(image_set, CR_or_GET = "CR", ...)
}

#' 2D global envelope test
#'
#' Provides global envelope tests based on 2D functions
#'
#' @references
#' Mrkvička, T., Soubeyrand, S., Myllymäki, M., Grabarnik, P. and Hahn, U. (2016). Monte Carlo testing in spatial statistics, with applications to spatial residuals. Spatial Statistics 18, Part A: 40--53. doi: 10.1016/j.spasta.2016.04.005
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' Mrkvička, T., Hahn, U. and Myllymäki, M. (2018). A one-way ANOVA test for functional data with graphical interpretation. arXiv:1612.03608 [stat.ME]
#'
#' @inheritParams central_region_2d
#' @param ... Additional parameters to be passed to \code{\link{global_envelope_test}}.
#' @return An object of class "global_envelope_2d" (and "list"),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item rx = the vector of values of the argument rx on x-coordinate at which the test was made
#' \item rx = the vector of values of the argument ry on y-coordinate at which the test was made
#' \item obs = the data function (matrix), if there is only one data function. Otherwise not existing.
#' \item lo = the lower envelope (matrix) based on the simulated functions
#' \item hi = the upper envelope (matrix) based on the simulated functions
#' \item central = the central curve, mean (default) or median of the test functions
#' T_i(r), i=2, ..., s+1. Used for visualization only.
#' }
#' Additionally, the return value has attributes \code{method}, \code{type}, \code{alternative},
#' \code{k_alpha}, \code{alpha}, \code{k}, \code{p} (and \code{p_interval} and \code{ties} if
#' type is \code{'rank'}) and \code{call}, see more detailed description in
#' \code{\link{global_envelope_test}}.
#' @export
#' @examples
#' # Example of spatial point pattern residuals
#' #-------------------------------------------
#' if(require(spatstat, quietly=TRUE)) {
#'   data(cells)
#'   X <- cells
#'   # Fit the hard-core process
#'   model <- ppm(X, interaction=Hardcore())
#'   summary(model)
#'   HD <- 0.08168525 # Hard-core process
#'   # Choose a bandwitdh by Scott's rule of thumb
#'   ds <- bw.scott(X); ds
#'   # Calculate raw residuals of the fitted model
#'   u <- diagnose.ppm(model, type="raw", rbord = HD, which ="smooth",
#'                     sigma=ds, plot.it=FALSE)
#'   obs <- u$smooth$Z$v
#'   # Generate simulations from the hard-core null model
#'   nsim <- 1999
#'   simulations <- NULL
#'   ext.factor <- max(X$window$xrange[2]-X$window$xrange[1],
#'                     X$window$yrange[2]-X$window$yrange[1]) / 10
#'   win.extend <- owin(c(X$window$xrange[1]-ext.factor, X$window$xrange[2]+ext.factor),
#'                      c(X$window$yrange[1]-ext.factor, X$window$yrange[2]+ext.factor))
#'   mod02 <- list(cif="hardcore", par=list(beta=exp(model$fitin$coefs[1]),hc=HD), w=win.extend)
#'   # starting point pattern in an extended window
#'   x.start <- runifpoint(X$n, win=win.extend)
#'   # simulations
#'   for(sss in 1:nsim){
#'     uppp <- rmh(model=mod02, start=list(x.start=x.start), control=list(p=1,nrep=1e5,nverb=5000))
#'     f <- uppp$x > X$window$xrange[1] & uppp$x < X$window$xrange[2] &
#'          uppp$y > X$window$yrange[1] & uppp$y < X$window$yrange[2]
#'     simulations[[sss]] <- ppp(uppp$x[f], uppp$y[f], window=X$window)
#'   }
#'   # Calculate the raw residuals for simulations
#'   sim <- array(0, c(nsim, u$smooth$Z$dim))
#'   for(i in 1:length(simulations)) {
#'     model=ppm(simulations[[i]],interaction=Hardcore(HD));
#'     u_sim <- diagnose.ppm(model, type="raw", rbord = HD, which ="smooth", sigma=ds, plot.it=FALSE)
#'     sim[i,,] <- u_sim$smooth$Z$v
#'     if((i %% 100)==0) cat(i, ' ')
#'   }
#'   # Constract the global envelope test for the (2D) raw residuals
#'   iset <- create_image_set(list(obs=obs, sim_m=sim))
#'   res <- global_envelope_test_2d(iset)
#'   par(mfrow=c(2,3))
#'   plot(res)
#'   # The same colors for all plots:
#'   par(mfrow=c(2,3))
#'   plot(res, col=spatstat::colourmap(grDevices::gray(0:255/255), range=c(min(res$obs, res$lo), max(res$obs, res$hi))))
#' }
global_envelope_test_2d <- function(image_set, ...) {
  if(class(image_set)[1] != "image_set") stop("image_set should be an object of class image_set.\n")
  cr_or_GET_2d(image_set, CR_or_GET = "GET", ...)
}

#' Print method for the class 'global_envelope_2d'
#' @usage \method{print}{global_envelope_2d}(x, ...)
#'
#' @param x an 'global_envelope_2d' object
#' @param ... Ignored.
#'
#' @method print global_envelope_2d
#' @export
print.global_envelope_2d <- function(x, ...) {
  print.global_envelope(x, ...)
}

#' Plot method for the class 'global_envelope_2d'
#'
#' Plot method for the class 'global_envelope_2d'
#'
#'
#' Additional parameter \code{col} can be passed in \code{...} to \code{\link[spatstat]{plot.im}}.
#' If \code{col} not given, a \code{\link[spatstat]{colourmap}} of 255 grey values between the
#' minimum and maximum of the function values is used for each image separately.
#' If \code{col} is provided, the same specification will be used for each produced plot,
#' which may make it easier to compare the figures with each other.
#'
#' @usage \method{plot}{global_envelope_2d}(x, col, sign.col = c(255, 0, 0), transparency = 85, main, ...)
#'
#' @param x an 'global_envelope_2d' object
#' @param sign.col A vector of length 3 giving the color for the significant regions.
#' The elements should correspond to red, green, blue of \code{\link[grDevices]{rgb}}.
#' Default to red (255, 0, 0).
#' @param transparency A number between 0 and 255 (default 85, 33% transparency).
#' Corresponds to alpha of \code{\link[grDevices]{rgb}}. Used in plotting the significant regions.
#' @param main A vector of length 5 giving titles for the 4 (if there is an observed function) or
#' 5 plots. A default exists, but this argument allows modifying it.
#' The default is \code{c("Observed", "Lower envelope", "Upper envelope",
#' "Significance: below (red)", "Significance: above (red)")}.
#' @param ... Additional parameters to be passed to \code{\link[spatstat]{plot.im}}.
#'
#' @method plot global_envelope_2d
#' @export
#' @importFrom grDevices rgb
#' @importFrom grDevices gray
#' @importFrom spatstat as.im
#' @importFrom spatstat contour.im
plot.global_envelope_2d <- function(x, sign.col = c(255, 0, 0), transparency = 85, main, ...) {
  extraargs <- list(...)
  if(length(sign.col)!=3) stop("Unreasonable length of sign.col.\n")
  if(missing(main)) main <- c("Observed", "Lower envelope", "Upper envelope",
                              "Significance: below (red)",
                              "Significance: above (red)")
  if(length(main) != 5) {
    warning("Unreasonable main provided. Setting empty main(s).\n")
    main <- rep("", times=5)
  }
  # First plot
  if(!is.null(x[['obs']])) {
    obs.im <- spatstat::as.im(list(x=x$rx, y=x$ry, z= x$obs))
    if(!("col" %in% names(extraargs))) {
      col <- spatstat::colourmap(grDevices::gray(0:255/255), range=range(x$obs))
      spatstat::plot.im(obs.im, col=col, main=main[1], ...)
    }
    else spatstat::plot.im(obs.im, main=main[1], ...)
    spatstat::contour.im(obs.im, add=TRUE)
  }
  # Lower envelope
  if(attr(x, "alternative") != "greater") {
    lo.im <- spatstat::as.im(list(x=x$rx, y=x$ry, z= x$lo))
    if(!("col" %in% names(extraargs))) {
      if(max(x$lo)>min(x$lo))
        col <- spatstat::colourmap(grDevices::gray(0:255/255), range=range(x$lo))
      else col <- grDevices::gray(0)
      spatstat::plot.im(lo.im, col=col, main=main[2], ...)
    }
    else spatstat::plot.im(lo.im, main=main[2], ...)
    if(!is.character(col)) spatstat::contour.im(lo.im, add=TRUE)
  }
  # Upper envelope
  if(attr(res, "alternative") != "less") {
    hi.im <- spatstat::as.im(list(x=x$rx, y=x$ry, z= x$hi))
    if(!("col" %in% names(extraargs))) {
      if(max(x$hi)>min(x$hi))
        col <- spatstat::colourmap(grDevices::gray(0:255/255), range=range(x$hi))
      else col <- grDevices::gray(1)
      spatstat::plot.im(hi.im, col=col, main=main[3], ...)
    }
    else spatstat::plot.im(hi.im, main=main[3], ...)
    if(!is.character(col)) spatstat::contour.im(hi.im, add=TRUE)
  }
  # Significance
  transparent <- grDevices::rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")
  red <- grDevices::rgb(sign.col[1], sign.col[2], sign.col[3], max = 255, alpha = transparency, names = "red")
  # Below
  if(attr(x, "alternative") != "greater") {
    if(!("col" %in% names(extraargs))) {
      col <- spatstat::colourmap(grDevices::gray(0:255/255), range=range(x$obs))
      spatstat::plot.im(obs.im, col=col, main=main[4], ...)
    }
    else spatstat::plot.im(obs.im, main=main[4], ...)
    if(sum(x$obs < x$lo) > 0)
      spatstat::plot.im(spatstat::im(x$obs < x$lo, xcol=x$rx, yrow=x$ry),
                        col=c(transparent, red), add=TRUE)
  }
  # Above
  if(attr(res, "alternative") != "less") {
    if(!("col" %in% names(extraargs))) {
      col <- spatstat::colourmap(grDevices::gray(0:255/255), range=range(x$obs))
      spatstat::plot.im(obs.im, col=col, main=main[5], ...)
    }
    else spatstat::plot.im(obs.im, main=main[5], ...)
    if(sum(x$obs > x$hi) > 0)
      spatstat::plot.im(spatstat::im(x$obs > x$hi, xcol=x$rx, yrow=x$ry),
                        col=c(transparent, red), add=TRUE)
  }
}
