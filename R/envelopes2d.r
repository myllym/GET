# Additional parameter \code{col} can be passed in \code{...} to \code{\link[spatstat]{plot.im}}.
# If \code{col} not given, a \code{\link[spatstat]{colourmap}} of 255 grey values between the
# minimum and maximum of the function values is used for each image separately.
# If \code{col} is provided, the same specification will be used for each produced plot,
# which may make it easier to compare the figures with each other.
#
# @inheritParams plot.global_envelope
# @param x an 'global_envelope2d' object
# @param plot_style Either "ggplot2" or "basic". (Similar to the argument in \code{\link{plot.global_envelope}}.)
# @param fixedscales Logical. TRUE for the same scales for all images.
# @param sign.col The color for the significant regions. Default to "red".
# @param transparency A number between 0 and 1 (default 85/255, 33% transparency).
# Similar to alpha of \code{\link[grDevices]{rgb}}. Used in plotting the significant regions.
# @param contours Logical. Whether to draw contour lines to the observed function and the lower and upper envelope.
# @param main The overall main. Default exists.
# @param ... Additional parameters to be passed to \code{\link[spatstat]{plot.im}}.
#' @importFrom ggplot2 facet_wrap ggtitle
#' @importFrom gridExtra grid.arrange
plot_global_envelope2d <- function(x, plot_style = c("ggplot2", "basic"),
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
           dfs <- env2d_ggplot2_helper(x, fixedscales=fixedscales)
           if(fixedscales) {
             g <- env2d_ggplot2_helper_1(dfs, sign.col, transparency, contours)
             g <- g + facet_wrap(vars(label))
             g + ggtitle(main)
           } else {
             gs = env2d_ggplot2_helper_many_single_plots(dfs, sign.col, transparency, contours)
             p1 = grid.arrange(grobs=gs, nrow=ceiling(length(gs)/3), top=main)
           }
         })
}

# Additional parameter \code{col} can be passed in \code{...} to \code{\link[spatstat]{plot.im}}.
# If \code{col} not given, a \code{\link[spatstat]{colourmap}} of 255 grey values between the
# minimum and maximum of the function values is used for each image separately.
# If \code{col} is provided, the same specification will be used for each produced plot,
# which may make it easier to compare the figures with each other.
#
# If fixedscales is FALSE (or 0) all images will have separate scale.
# If fixedscales is TRUE (or 1) each x[[i]] will have a common scale.
# If fixedscales is 2 all images will have common scale.
#
# @param fixedscales 0, 1 or 2. See details.
#' @importFrom ggplot2 facet_grid vars
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices col2rgb
plot_combined_global_envelope2d <- function(x, plot_style = c("ggplot2", "basic"),
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
           dfs <- mapply(x, names(x), SIMPLIFY=FALSE, FUN=function(x, main) {
             env2d_ggplot2_helper(x, fixedscales=fixedscales, contours=contours, main=main, insertmain=!fixedscales)
           })
           if(fixedscales==2) {
             df <- do.call(rbind, dfs)
             g <- env2d_ggplot2_helper_1(df, sign.col, transparency, contours) + facet_grid(rows=vars(main), cols=vars(label))
             g + ggtitle(fullmain)
           } else if(fixedscales==1) {
             gs <- lapply(dfs, function(df) {
               env2d_ggplot2_helper_1(df, sign.col, transparency, contours) + facet_grid(rows=vars(main), cols=vars(label))
             })
             grid.arrange(grobs=gs, ncol=1, top=fullmain)
           } else if(fixedscales==0) {
             gs <- lapply(dfs, function(dfs2) {
               gs = env2d_ggplot2_helper_many_single_plots(dfs2, sign.col, transparency, contours)
               grid.arrange(grobs=gs, nrow=1)
             })
             grid.arrange(grobs=gs, ncol=1, top=fullmain)
           } else { stop("Invalid fixedscales argument.") }
         })
}
