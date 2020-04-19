# 2d plots with ggplot2
#----------------------
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

# Add a part "what" of the envelope2d plot to the given ggplot object g.
# The ggplot object will have faceting variable "label".
addplot_global_envelope2d <- function(g, x, what, sign.col, transparency) {
  namelist <- list(obs = "Observed",
                   lo = "Lower envelope" ,
                   hi = "Upper envelope" ,
                   lo.sign = "Sign.: below" ,
                   hi.sign = "Sign.: above" )
  for(w in what) {
    df <- x
    df$label <- factor(namelist[w])
    if(w == "lo.sign") {
      w <- "obs"
      df$signif <- df$obs < df$lo
    } else if(w == "hi.sign") {
      w <- "obs"
      df$signif <- df$obs > df$hi
    }
    g <- g + choose_geom(df, varfill=w)
    if(any(df$signif))
      g <- g + choose_geom(df[df$signif,], fill=sign.col, alpha=transparency)
  }
  g
}

# Drop options that are not relevant to the alternative hypothesis.
checkarg_envelope2d_what <- function(x, what) {
  alt <- get_alternative(x)
  if(alt == "greater") what <- setdiff(what, c("lo", "lo.sign"))
  else if(alt == "less") what <- setdiff(what, c("hi", "hi.sign"))
  what
}

# Used also when creating a plot without fixedscales.
#' @importFrom ggplot2 ggplot facet_wrap labs coord_fixed
plot_global_envelope2d_fixedscales <- function(x, what=c("obs", "hi", "lo", "hi.sign", "lo.sign"),
                                   sign.col = "red", transparency = 85/255) {
  g <- ggplot()
  g <- addplot_global_envelope2d(g, x, what, sign.col, transparency)
  g + facet_wrap("label") + labs(x="", y="", fill="") + coord_fixed()
}

# Used also when creating a plot without fixedscales.
# @param sign.col The color for the significant regions. Default to "red".
# @param transparency The transparency of the significant regions.
# A number between 0 and 1 (default 85/255, 33% transparency).
#' @importFrom ggplot2 ggplot facet_grid labs coord_fixed
plot_combined_global_envelope2d_fixedscales <- function(x, what=c("obs", "hi", "lo", "hi.sign", "lo.sign"),
                                            sign.col = "red", transparency = 85/255) {
  g <- ggplot()
  for(i in seq_along(x)) {
    df <- x[[i]]
    df$main <- names(x)[i]
    g <- addplot_global_envelope2d(g, df, what, sign.col, transparency)
  }
  g + facet_grid(main ~ label) + labs(x="", y="", fill="") + coord_fixed()
}

# Plotting function for 2d global envelopes
#
# @inheritParams plot.global_envelope
# @param x A 'global_envelope' object for two-dimensional functions
# @param fixedscales Logical. TRUE for the same scales for all images.
# @param main The overall main.
# @param what Character vector specifying what information should be plotted.
#   A combination of:
#   Observed (\code{"obs"}), upper envelope (\code{"hi"}), lower envelope (\code{"lo"}),
#   observed with significantly higher values highlighted (\code{"hi.sign"}),
#   observed with significantly lower values highlighted (\code{"lo.sign"}).
# @param ... Additional parameters to be passed to \code{plot_global_envelope2d_fixedscales}.
#' @importFrom ggplot2 ggtitle
#' @importFrom gridExtra grid.arrange
plot_global_envelope2d <- function(x, fixedscales = TRUE, main,
                                   what=c("obs", "hi", "lo", "hi.sign", "lo.sign"), ...) {
  what <- match.arg(what, several.ok = TRUE)
  what <- checkarg_envelope2d_what(x, what)

  if(fixedscales) {
    g <- plot_global_envelope2d_fixedscales(x, what, ...)
    g + ggtitle(main)
  } else {
    gs <- lapply(what, function(w) plot_global_envelope2d_fixedscales(x, w, ...))
    grid.arrange(grobs=gs, nrow=ceiling(length(gs)/3), top=main)
  }
}

# Plotting function for combined 2d global envelopes
#
# If fixedscales is FALSE (or 0) all images will have separate scale.
# If fixedscales is TRUE (or 1) each x[[i]] will have a common scale.
# If fixedscales is 2 all images will have common scale.
#
# @inheritParams plot_global_envelope2d
# @param fixedscales 0, 1 or 2. See details.
#' @importFrom ggplot2 facet_grid facet_wrap ggtitle
#' @importFrom ggplot2 label_value theme element_blank
#' @importFrom gridExtra grid.arrange
plot_combined_global_envelope2d <- function(x, fixedscales = 2, main,
                                            what=c("obs", "hi", "lo", "hi.sign", "lo.sign"), ...) {
  what <- match.arg(what, several.ok = TRUE)
  what <- checkarg_envelope2d_what(x[[1]], what)

  if(fixedscales==2) {
    g <- plot_combined_global_envelope2d_fixedscales(x, what, ...)
    g + ggtitle(main)
  } else if(fixedscales==1) {
    gs <- lapply(seq_along(x), function(i) {
      env <- x[[i]]
      env$main <- names(x)[i]
      g <- plot_global_envelope2d_fixedscales(env, what, ...)
      g + facet_grid(main~label)
    })
    grid.arrange(grobs=gs, ncol=1, top=main)
  } else if(fixedscales==0) {
    f <- function(env, main, what) {
      env$main <- main
      g <- plot_global_envelope2d_fixedscales(env, what, ...)
      g + facet_wrap(c("main", "label"), labeller=function(l) label_value(l, multi_line = FALSE))
      g + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
    }
    gs <- outer(seq_along(x), c("obs", "hi"), FUN=Vectorize(function(i, w) {
      list(f(x[[i]], names(x)[i], w))
    }))
    grid.arrange(grobs=t(gs), top=main)
  } else { stop("Invalid fixedscales argument.") }
}
