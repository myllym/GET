# 2d plots with ggplot2
#----------------------
globalVariables(c("main", "label"))

# Choose ggplot2 geom based on variables found in df
# @param varfill (Optional) Name of the variable used for 'fill' aesthetic.
# (fill is always specified, but for a fixed color it is given in ...)
#' @importFrom ggplot2 geom_tile geom_rect aes .data
choose_geom <- function(df, varfill, ...) {
  if(!is.null(df$width)) {
    if(all(df$width[1] == df$width) && all(df$height[1] == df$height)) {
      w <- df$width[1]
      h <- df$height[1]
      if(!missing(varfill))
        geom_raster_fixed(data=df, aes(x=.data$x, y=.data$y, fill=.data[[varfill]]), width=w, height=h, ...)
      else
        geom_raster_fixed(data=df, aes(x=.data$x, y=.data$y), width=w, height=h, ...)
    } else {
      if(!missing(varfill))
        geom_tile(data=df, aes(x=.data$x, y=.data$y, width=.data$width,
                               height=.data$height, fill=.data[[varfill]]), ...)
      else
        geom_tile(data=df, aes(x=.data$x, y=.data$y, width=.data$width,
                               height=.data$height), ...)
    }
  } else {
    if(!missing(varfill))
      geom_rect(data=df, aes(xmin=.data$xmin, ymin=.data$ymin, xmax=.data$xmax,
                             ymax=.data$ymax, fill=.data[[varfill]]), ...)
    else
      geom_rect(data=df, aes(xmin=.data$xmin, ymin=.data$ymin, xmax=.data$xmax,
                             ymax=.data$ymax), ...)
  }
}

# A helper function for env2d_ggplot2. Produces a ggplot with the significant region
#' @importFrom ggplot2 ggplot aes geom_contour coord_fixed .data labs
env2d_ggplot2_helper_1 <- function(df, sign.col, transparency, contours = TRUE) {
  g <- ggplot() + choose_geom(df, varfill='z')
  if(any(df$signif))
    g <- g + choose_geom(df[df$signif,], fill=sign.col, alpha=transparency)
  if(contours && !is.null(df$x)) g <- g + geom_contour(data=df[df$contour,], aes(x=.data$x, y=.data$y, z=.data$z))
  g <- g + coord_fixed(ratio=1)
  g <- g + labs(x="", y="", fill="")
  g
}

#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 facet_wrap ggtitle theme element_blank vars
env2d_ggplot2_helper <- function(x, fixedscales, contours = TRUE, main="", insertmain=TRUE) {
  namelist <- list(obs = "Observed",
                   lo = "Lower envelope" ,
                   hi = "Upper envelope" ,
                   lo.sign = "Sign.: below" ,
                   hi.sign = "Sign.: above" )
  if(!missing(main) && !is.null(main) && insertmain) {
    for (i in seq_along(namelist)) {
      namelist[[i]] <- paste(main, ": ", namelist[[i]])
    }
  }

  # If curve_set$r was created using a data.frame
  if(!is.null(x[['x']])) df <- x[, c("height", "width", "x", "y")]
  else if(!is.null(x[['xmin']])) df <- x[, c("xmax", "xmin", "ymax", "ymin")]
  else stop("Cannot detect curve_set r")

  adddf <- function(df, z, name, label=namelist[[name]], contour=FALSE, signif=FALSE) {
    df$z <- c(z)
    df$name <- name
    df$label <- factor(label)
    df$contour <- contour
    df$signif <- c(signif)
    df$main <- main
    df
  }
  alt <- get_alternative(x)
  dfs <- list()
  if(!is.null(x$obs)) {
    dfs <- c(dfs, list(adddf(df, x$obs, "obs", contour=contours)))
  }
  if(alt != "greater") {
    dfs <- c(dfs, list(adddf(df, x$lo, "lo", contour=contours)))
  }
  if(alt != "less") {
    dfs <- c(dfs, list(adddf(df, x$hi, "hi", contour=contours)))
  }
  if(!is.null(x$obs) && alt != "greater") {
    dfs <- c(dfs, list(adddf(df, x$obs, "lo.sign", signif=x$obs < x$lo)))
  }
  if(!is.null(x$obs) && alt != "less") {
    dfs <- c(dfs, list(adddf(df, x$obs, "hi.sign", signif=x$obs > x$hi)))
  }
  if(fixedscales)
    do.call(rbind, dfs)
  else
    dfs
}

# @param sign.col The color for the significant regions.
# @param transparency A number between 0 and 1.
# Similar to alpha of \code{\link[grDevices]{rgb}}. Used in plotting the significant regions.
#' @importFrom ggplot2 theme element_blank facet_wrap vars
env2d_ggplot2_helper_many_single_plots <- function(dfs, sign.col, transparency, contours = TRUE) {
  remove_axes_theme <- theme(axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.title.y=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank())
  lapply(dfs, function(df) {
    g <- env2d_ggplot2_helper_1(df, sign.col, transparency, contours)
    g <- g + facet_wrap(vars(label))
    g + remove_axes_theme
  })
}

# Plotting function for 2d global envelopes
#
# @inheritParams plot.global_envelope
# @param x A 'global_envelope' object for two-dimensional functions
# @param fixedscales Logical. TRUE for the same scales for all images.
# @param sign.col The color for the significant regions. Default to "red".
# @param transparency A number between 0 and 1 (default 85/255, 33% transparency).
# Similar to alpha of \code{\link[grDevices]{rgb}}. Used in plotting the significant regions.
# @param contours Logical. Whether to draw contour lines to the observed function and the lower and upper envelope.
# @param main The overall main. Default exists.
# @param ... Additional parameters to be passed to \code{\link[spatstat]{plot.im}}.
#' @importFrom ggplot2 facet_wrap ggtitle
#' @importFrom gridExtra grid.arrange
plot_global_envelope2d <- function(x, fixedscales = TRUE, sign.col = "red", transparency = 85/255,
                                   contours = FALSE, main = NULL, digits = 3, ...) {
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
}

# Plotting function for combined 2d global envelopes
#
# If fixedscales is FALSE (or 0) all images will have separate scale.
# If fixedscales is TRUE (or 1) each x[[i]] will have a common scale.
# If fixedscales is 2 all images will have common scale.
#
# @param fixedscales 0, 1 or 2. See details.
#' @importFrom ggplot2 facet_grid vars
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices col2rgb
plot_combined_global_envelope2d <- function(x, fixedscales = 2, sign.col = "red", transparency = 85/255,
                                            contours = FALSE, main = NULL, digits = 3, ...) {
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
}
