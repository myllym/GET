# res_v = global envelope test results for images (in image_sets)
# that have been transformed to vectors for testing
curve_set_results_to_image_results <- function(res_v, image_sets) {
  if(!(class(res_v)[1] %in% c("global_envelope", "combined_global_envelope"))) stop("Invalid res_v.\n")
  if(class(image_sets)[1] != "list") image_sets <- list(image_sets)
  obs_d <- dim(image_sets[[1]]$obs)
  switch(class(res_v)[1],
    global_envelope = {
      if(!is.null(res_v$obs)) {
        res <- structure(list(r=list(rx=image_sets[[1]]$r[[1]], ry=image_sets[[1]]$r[[2]]),
                              obs=matrix(res_v$obs, nrow=obs_d[1], ncol=obs_d[2]),
                              central=matrix(res_v$central, nrow=obs_d[1], ncol=obs_d[2]),
                              lo=matrix(res_v$lo, nrow=obs_d[1], ncol=obs_d[2]),
                              hi=matrix(res_v$hi, nrow=obs_d[1], ncol=obs_d[2])))
      }
      else {
        res <- structure(list(r=list(rx=image_sets[[1]]$r[[1]], ry=image_sets[[1]]$r[[2]]),
                              central=matrix(res_v$central, nrow=obs_d[1], ncol=obs_d[2]),
                              lo=matrix(res_v$lo, nrow=obs_d[1], ncol=obs_d[2]),
                              hi=matrix(res_v$hi, nrow=obs_d[1], ncol=obs_d[2])))
      }
      mostattributes(res) <- attributes(res_v)
      attr(res, "xlab") <- attr(res, "xexp") <- attr(res, "ylab") <- attr(res, "xexp") <- attr(res, "call") <- NULL
      class(res) <- c("global_envelope2d", "list")
    },
    combined_global_envelope = {
      res <- list()
      for(i in 1:length(res_v)) {
        res[[i]] <- structure(list(r=list(rx=image_sets[[1]]$r[[1]], ry=image_sets[[1]]$r[[2]]),
                                   obs=matrix(res_v[[i]]$obs, nrow=obs_d[1], ncol=obs_d[2]),
                                   central=matrix(res_v[[i]]$central, nrow=obs_d[1], ncol=obs_d[2]),
                                   lo=matrix(res_v[[i]]$lo, nrow=obs_d[1], ncol=obs_d[2]),
                                   hi=matrix(res_v[[i]]$hi, nrow=obs_d[1], ncol=obs_d[2])))
        mostattributes(res[[i]]) <- attributes(res_v[[i]])
        attr(res[[i]], "xlab") <- attr(res[[i]], "xexp") <- attr(res[[i]], "ylab") <- attr(res[[i]], "xexp") <- attr(res[[i]], "call") <- NULL
        class(res[[i]]) <- c("global_envelope2d", "list")
      }
      names(res) <- names(res_v)
      mostattributes(res) <- attributes(res_v)
      attr(res, "xlab") <- attr(res, "xexp") <- attr(res, "ylab") <- attr(res, "xexp") <- attr(res, "call") <- NULL
      class(res) <- c("combined_global_envelope2d", "list")
    }
  )
  res
}

# Functionality for central_region2d and global_envelope_test2d
cr_or_GET_2d <- function(image_sets, CR_or_GET = c("CR", "GET"), ...) {
  CR_or_GET <- match.arg(CR_or_GET)
  if(inherits(image_sets, "image_set")) image_sets <- list(image_sets)
  if(!all(sapply(image_sets, function(x) { inherits(x, "image_set") }))) stop("The list does not contain image_set objects.\n")
  obs_d <- lapply(image_sets, function(x) { dim(x$obs) })
  sim_d <- lapply(image_sets, function(x) { dim(x$sim_m) })
  # Check that dimensions of different image sets are the same
  if(!all(sapply(obs_d, FUN=identical, y=obs_d[[1]]))) stop("Dimensions of image sets (obs) do not match.\n")
  if(!all(sapply(sim_d, FUN=identical, y=sim_d[[1]]))) stop("Dimensions of image sets (sim_m) do not match.\n")
  # Check dimensions of each image set
  image_sets <- lapply(image_sets, check_image_set_dimensions)
  # Check equalities of the r values
  if(!all(sapply(image_sets, FUN = function(x) {identical(x$r, y=image_sets[[1]]$r)}))) stop("The r values of image sets do not match.\n")
  # Create curve sets transforming the 2d functions (matrices) to vectors
  curve_sets_v <- lapply(image_sets, image_set_to_curve_set)

  # Prepare global envelope or global envelope test
  switch(CR_or_GET,
         CR = {
           res_v <- central_region(curve_sets_v, ...)
         },
         GET = {
           res_v <- global_envelope_test(curve_sets_v, ...)
         })

  # Transform the results to 2d
  res <- curve_set_results_to_image_results(res_v, image_sets)
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
#' @param image_sets An image set, i.e. a set of 2d functions (or a list of them).
#' See \code{\link{create_image_set}}.
#' @param ... Additional parameters to be passed to \code{\link{central_region}}.
#' @return An object of class "global_envelope2d" (and "list"),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = a list of vectors of values of x- and y-coordinates at which the test was made
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
central_region2d <- function(image_sets, ...) {
  res <- cr_or_GET_2d(image_sets, CR_or_GET = "CR", ...)
  attr(res, "call") <- match.call()
  res
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
#' @inheritParams central_region2d
#' @param ... Additional parameters to be passed to \code{\link{global_envelope_test}}.
#' @return An object of class "global_envelope2d" (and "list"),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = a list of vectors of values of x- and y-coordinates at which the test was made
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
#' @aliases GET.2d
#' @examples
#' \donttest{
#' # Example of spatial point pattern residuals
#' #-------------------------------------------
#' if(require("spatstat", quietly=TRUE)) {
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
#'   nsim <- 499 # Number of simulations; increase for serious analysis!
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
#'   sim <- array(0, c(u$smooth$Z$dim, nsim))
#'   for(i in 1:length(simulations)) {
#'     model=ppm(simulations[[i]],interaction=Hardcore(HD));
#'     u_sim <- diagnose.ppm(model, type="raw", rbord = HD, which ="smooth", sigma=ds, plot.it=FALSE)
#'     sim[,,i] <- u_sim$smooth$Z$v
#'     if((i %% 100)==0) cat(i, ' ')
#'   }
#'   # Constract the global envelope test for the (2D) raw residuals
#'   iset <- create_image_set(list(obs=obs, sim_m=sim))
#'   res <- global_envelope_test2d(iset, type="area")
#'   plot(res)
#'   plot(res, contours=FALSE) + ggplot2::scale_fill_gradient(low="black", high="white")
#'   plot(res, fixedscales=FALSE)
#' }}
global_envelope_test2d <- function(image_sets, ...) {
  res <- cr_or_GET_2d(image_sets, CR_or_GET = "GET", ...)
  attr(res, "call") <- match.call()
  res
}

#' Print method for the class 'global_envelope2d'
#'
#' @param x an 'global_envelope2d' object
#' @param ... Ignored.
#' @export
print.global_envelope2d <- function(x, ...) {
  GEprinthelper(x, ...)
}

#' Print method for the class 'combined_global_envelope2d'
#'
#' @param x an 'combined_global_envelope2d' object
#' @param ... Ignored.
#' @export
print.combined_global_envelope2d <- function(x, ...) {
  GEprinthelper(attr(x, "level2_ge"), ...)
}


#' Plot method for the class 'global_envelope2d'
#'
#' Plot method for the class 'global_envelope2d'
#'
#'
#' Additional parameter \code{col} can be passed in \code{...} to \code{\link[spatstat]{plot.im}}.
#' If \code{col} not given, a \code{\link[spatstat]{colourmap}} of 255 grey values between the
#' minimum and maximum of the function values is used for each image separately.
#' If \code{col} is provided, the same specification will be used for each produced plot,
#' which may make it easier to compare the figures with each other.
#'
#' @inheritParams plot.global_envelope
#' @param x an 'global_envelope2d' object
#' @param plot_style Either "ggplot2" or "basic". (Similar to the argument in \code{\link{plot.global_envelope}}.)
#' @param fixedscales Logical. TRUE for the same scales for all images.
#' @param sign.col The color for the significant regions. Default to "red".
#' @param transparency A number between 0 and 1 (default 85/255, 33% transparency).
#' Similar to alpha of \code{\link[grDevices]{rgb}}. Used in plotting the significant regions.
#' @param contours Logical. Whether to draw contour lines to the observed function and the lower and upper envelope.
#' @param main The overall main. Default exists.
#' @param ... Additional parameters to be passed to \code{\link[spatstat]{plot.im}}.
#'
#' @export
#' @importFrom ggplot2 facet_wrap ggtitle
#' @importFrom gridExtra grid.arrange
plot.global_envelope2d <- function(x, plot_style = c("ggplot2", "basic"),
                                    fixedscales = TRUE, sign.col = "red", transparency = 85/255,
                                    contours = FALSE, main = NULL, digits = 3, ...) {
  plot_style <- match.arg(plot_style)
  switch(plot_style,
         basic = {
           sign.col <- as.numeric(col2rgb(sign.col, alpha = FALSE))
           transparency <- 255*transparency
           env2d_basic_plot(x, var='obs', sign.col=sign.col, transparency=transparency, contours=contours, ...)
           env2d_basic_plot(x, var='lo', sign.col=sign.col, transparency=transparency, contours=contours, ...)
           env2d_basic_plot(x, var='hi', sign.col=sign.col, transparency=transparency, contours=contours, ...)
           env2d_basic_plot(x, var='lo.sign', sign.col=sign.col, transparency=transparency, contours=contours, ...)
           env2d_basic_plot(x, var='hi.sign', sign.col=sign.col, transparency=transparency, contours=contours, ...)
         },
         ggplot2 = {
           if(is.null(main)) main <- env_main_default(x, digits=digits)
           w <- x$r$rx[2]-x$r$rx[1]
           h <- x$r$ry[2]-x$r$ry[1]
           dfs <- env2d_ggplot2_helper(x, fixedscales=fixedscales)
           if(fixedscales) {
             g <- env2d_ggplot2_helper_1(dfs, w, h, sign.col, transparency, contours)
             g <- g + facet_wrap(vars(label))
             g + ggtitle(main)
           } else {
             gs = env2d_ggplot2_helper_many_single_plots(dfs, w, h, sign.col, transparency, contours)
             p1 = grid.arrange(grobs=gs, nrow=ceiling(length(gs)/3), top=main)
           }
         })
}

#' Plot method for the class 'combined_global_envelope2d'
#'
#' Plot method for the class 'combined_global_envelope2d'
#'
#'
#' Additional parameter \code{col} can be passed in \code{...} to \code{\link[spatstat]{plot.im}}.
#' If \code{col} not given, a \code{\link[spatstat]{colourmap}} of 255 grey values between the
#' minimum and maximum of the function values is used for each image separately.
#' If \code{col} is provided, the same specification will be used for each produced plot,
#' which may make it easier to compare the figures with each other.
#'
#' If fixedscales is FALSE (or 0) all images will have separate scale.
#' If fixedscales is TRUE (or 1) each x[[i]] will have a common scale.
#' If fixedscales is 2 all images will have common scale.
#'
#' @param x an 'combined_global_envelope2d' object
#' @inheritParams plot.global_envelope2d
#' @param fixedscales 0, 1 or 2. See details.
#' @method plot combined_global_envelope2d
#' @export
#' @importFrom ggplot2 facet_grid vars
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices col2rgb
plot.combined_global_envelope2d <- function(x, plot_style = c("ggplot2", "basic"),
                                             fixedscales = 2, sign.col = "red", transparency = 85/255,
                                             contours = FALSE, main = NULL, digits = 3, ...) {
  plot_style <- match.arg(plot_style)
  switch(plot_style,
         basic = {
           sign.col <- as.numeric(col2rgb(sign.col, alpha = FALSE))
           transparency <- 255*transparency
           for(i in 1:length(x)) {
             main.plots <- paste(names(x)[i], ": ", c("Observed",
                                                      "Lower envelope",
                                                      "Upper envelope",
                                                      "Sign. below",
                                                      "Sign. above"), sep="")
             env2d_basic_plot(x[[i]], var='obs', sign.col=sign.col, transparency=transparency, contours=contours, main=main.plots[1], ...)
             env2d_basic_plot(x[[i]], var='lo', sign.col=sign.col, transparency=transparency, contours=contours, main=main.plots[2], ...)
             env2d_basic_plot(x[[i]], var='hi', sign.col=sign.col, transparency=transparency, contours=contours, main=main.plots[3], ...)
             env2d_basic_plot(x[[i]], var='lo.sign', sign.col=sign.col, transparency=transparency, contours=contours, main=main.plots[4], ...)
             env2d_basic_plot(x[[i]], var='hi.sign', sign.col=sign.col, transparency=transparency, contours=contours, main=main.plots[5], ...)
           }
         },
         ggplot2 = {
           if(is.null(names(x))) names(x) <- paste(1:length(x))
           if(is.null(main)) fullmain <- env_main_default(attr(x, "level2_ge"), digits=digits)
           else fullmain <- NULL
           w <- x[[1]]$r$rx[2]-x[[1]]$r$rx[1]
           h <- x[[1]]$r$ry[2]-x[[1]]$r$ry[1]
           dfs <- mapply(x, names(x), SIMPLIFY=FALSE, FUN=function(x, main) {
             env2d_ggplot2_helper(x, fixedscales=fixedscales, contours=contours, main=main, insertmain=!fixedscales)
           })
           if(fixedscales==2) {
             df <- do.call(rbind, dfs)
             g <- env2d_ggplot2_helper_1(df, w, h, sign.col, transparency, contours) + facet_grid(rows=vars(main), cols=vars(label))
             g + ggtitle(fullmain)
           } else if(fixedscales==1) {
             gs <- lapply(dfs, function(df) {
               env2d_ggplot2_helper_1(df, w, h, sign.col, transparency, contours) + facet_grid(rows=vars(main), cols=vars(label))
             })
             grid.arrange(grobs=gs, ncol=1, top=fullmain)
           } else if(fixedscales==0) {
             gs <- lapply(dfs, function(dfs2) {
               gs = env2d_ggplot2_helper_many_single_plots(dfs2, w, h, sign.col, transparency, contours)
               grid.arrange(grobs=gs, nrow=1)
             })
             grid.arrange(grobs=gs, ncol=1, top=fullmain)
           } else { stop("Invalid fixedscales argument.") }
         })
}
