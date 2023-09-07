# A helper function for reducing from multiple envelopes (lo/hi.xx, where xx are the levels)
# to a single, widest envelope. Currently employed for 2d plotting.
# Only lo/hi.xx are modified, no other components.
as_simple_envelope <- function(x) {
  simplify_lohi <- function(x, lohi="lo") {
    ind <- grep(lohi, names(x))
    if(length(ind) > 1) {
      tmp <- x[,ind[1]]
      x[,ind] <- NULL
      x[[lohi]] <- tmp
    }
    x
  }
  if(inherits(x, c("global_envelope", "global_envelope2d"))) {
    x <- simplify_lohi(x, "lo")
    x <- simplify_lohi(x, "hi")
  }
  if(inherits(x, c("combined_global_envelope", "combined_global_envelope2d"))) {
    for(i in seq_along(x)) {
      x[[i]] <- as_simple_envelope(x[[i]])
    }
  }
  x
}

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
#' @importFrom ggplot2 geom_point geom_contour scale_color_discrete guides scale_fill_gradient
#' @importFrom grDevices grey
addplot_global_envelope2d <- function(g, x, what, sign.col, transparency, sign.type) {
  namelist <- list(obs = "Observed",
                   lo = "Lower envelope" ,
                   hi = "Upper envelope" ,
                   lo.sign = "Sign.: below" ,
                   hi.sign = "Sign.: above" ,
                   obs.sign = "")
  for(w in what) {
    background <- w
    df <- x
    df$label <- factor(namelist[w])
    df$Direction <- NA
    df$Deviation <- NA
    df$signif.lo <- df$signif.hi <- FALSE
    #- Handle the one-sided cases
    if(all(is.finite(df$hi))) {
      if(all(is.finite(df$lo))) {
        df$scaler <- df$hi - df$lo
        whichenvelopes <- "both"
      }
      else {
        df$scaler <- df$hi - df$central
        whichenvelopes <- "hi only"
      }
    }
    else {
      if(all(is.finite(df$lo))) {
        df$scaler <- df$central - df$lo
        whichenvelopes <- "lo only"
      }
      else stop("Something seems to be wrong with lo and hi.")
    }
    if(w %in% c("lo.sign", "obs.sign")) {
      background <- "obs"
      df$signif.lo <- df$obs < df$lo
      df$Deviation.lo <- (df$lo - df$obs)/df$scaler
      if(sign.type == "circles") {
        df$Deviation[df$signif.lo] <- df$Deviation.lo[df$signif.lo]
      }
      df$Direction[df$signif.lo] <- "Below"
    }
    if(w %in% c("hi.sign", "obs.sign")) {
      background <- "obs"
      df$signif.hi <- df$obs > df$hi
      df$Deviation.hi <- (df$obs - df$hi)/df$scaler
      if(sign.type == "circles") {
        df$Deviation[df$signif.hi] <- df$Deviation.hi[df$signif.hi]
      }
      df$Direction[df$signif.hi] <- "Above"
    }
    g <- g + choose_geom(df, varfill=background)
    if(sign.type == "col") {
      if(any(df$signif.lo))
        g <- g + choose_geom(df[df$signif.lo,], fill=sign.col[1], alpha=transparency)
      if(any(df$signif.hi))
        g <- g + choose_geom(df[df$signif.hi,], fill=sign.col[2], alpha=transparency)
    } else if(sign.type == "circles") {
      if(any(df$signif.lo) || any(df$signif.hi)) {
        g <- g + geom_point(aes(x=.data$x, y=.data$y, size=.data$Deviation, colour=.data$Direction),
                            data=df[df$signif.lo | df$signif.hi,], shape=1)
      }
    } else if(sign.type == "contour") {
      if(any(df$signif.lo))
        g <- g + geom_contour(aes(x=.data$x, y=.data$y, z=.data$Deviation.lo), data=df, breaks=0, colour=sign.col[1])
      if(any(df$signif.hi))
        g <- g + geom_contour(aes(x=.data$x, y=.data$y, z=.data$Deviation.hi), data=df, breaks=0, colour=sign.col[2])
    }
  }
  g <- g + labs(size="Deviation")
  if(any(df$signif.lo) & any(df$signif.hi)) g <- g + scale_color_discrete(type = c(Below=sign.col[1], Above=sign.col[2]))
  else if(any(df$signif.lo)) g <- g + scale_color_discrete(type = c(Below=sign.col[1]))
  else if(any(df$signif.hi)) g <- g + scale_color_discrete(type = c(Above=sign.col[2]))
  if(w %in% c("lo.sign", "hi.sign")) g <- g + guides(colour = "none")
  g + scale_fill_gradient(low = grey(0.5), high = grey(1))
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
                                   sign.col, transparency = 85/255, sign.type) {
  g <- ggplot()
  g <- addplot_global_envelope2d(g, x, what, sign.col, transparency, sign.type)
  g + facet_wrap("label") + labs(x="", y="", fill="") + coord_fixed()
}

# Used also when creating a plot without fixedscales.
# @param sign.col Two colors for the significant regions (below and above).
# @param transparency The transparency of the significant regions.
# A number between 0 and 1 (default 85/255, 33% transparency).
#' @importFrom ggplot2 ggplot facet_grid labs coord_fixed
plot_combined_global_envelope2d_fixedscales <- function(x, what=c("obs", "hi", "lo", "hi.sign", "lo.sign"),
                                            sign.col, transparency = 85/255, sign.type) {
  g <- ggplot()
  for(i in seq_along(x)) {
    df <- x[[i]]
    df$main <- names(x)[i]
    g <- addplot_global_envelope2d(g, df, what, sign.col, transparency, sign.type)
  }
  g + facet_grid(main ~ label) + labs(x="", y="", fill="") + coord_fixed()
}

#' Plotting function for 2d global envelopes
#'
#' @description
#' If more than one envelope has been calculated (corresponding to several coverage/alpha),
#' only the largest one is plotted.
#' @param x A 'global_envelope' object for two-dimensional functions
#' @param fixedscales Logical. TRUE for the same scales for all images.
#' @param what Character vector specifying what information should be plotted for 2d functions.
#' A combination of:
#' Observed (\code{"obs"}), upper envelope (\code{"hi"}), lower envelope (\code{"lo"}),
#' observed with significantly higher values highlighted (\code{"hi.sign"}),
#' observed with significantly lower values highlighted (\code{"lo.sign"}),
#' observed with significantly (lower and higher) values highlighted (\code{"obs.sign"}).
#' Default to the last one. Combination c("obs", "lo", "hi", "lo.sign", "hi.sign") can
#' also be of interest (earlier default).
#' @param sign.type Either \code{"col"} for color showing the significant region, or
#' \code{"contour"} for colored contour showing the significant region, or
#' \code{"circles"} for plotting circles at locations where the observed function
#' exceeds the envelope: diameters proportional to (obs-hi)/(hi-lo) for values
#' above the envelope and (lo-obs)/(hi-lo) for values below the envelope.
#' In the one-sided (testing) case, the divisors are instead (hi-central) (case 'greater')
#' and (central-lo) (case 'less').
#' Default to \code{"circles"}.
#' @param sign.col A vector of length two giving the colors for significant parts
#' below the envelope (first value) and above the envelope (second value).
#' @param transparency A number between 0 and 1 (default 155/255, 60% transparency).
#' Similar to alpha of \code{\link[grDevices]{rgb}}. Used in plotting the significant regions for 2d
#' functions.
#' @param ... Ignored.
#' @inheritParams plot.global_envelope
#' @importFrom ggplot2 ggtitle
#' @importFrom gridExtra grid.arrange
#' @export
#' @seealso \code{\link{graph.flm}}
plot.global_envelope2d <- function(x, fixedscales = TRUE,
                                   what = c("obs.sign", "obs", "lo", "hi", "lo.sign", "hi.sign"),
                                   sign.type = c("circles", "contour", "col"),
                                   sign.col = c("blue", "red"), transparency = 155/255,
                                   digits = 3, ...) {
  x <- as_simple_envelope(x)
  if(missing(what)) what <- "obs.sign"
  else what <- match.arg(what, several.ok = TRUE)
  what <- checkarg_envelope2d_what(x, what)
  main <- env_main_default(x, digits=digits)
  sign.type <- match.arg(sign.type)
  if(length(sign.col)==1) sign.col <- rep(sign.col, times=2)

  if(fixedscales) {
    g <- plot_global_envelope2d_fixedscales(x, what, sign.col, transparency, sign.type)
    g + ggtitle(main)
  } else {
    gs <- lapply(what, function(w) { plot_global_envelope2d_fixedscales(x, w, sign.col, transparency, sign.type) })
    grid.arrange(grobs=gs, nrow=ceiling(length(gs)/3), top=main)
  }
}

#' Plotting function for combined 2d global envelopes
#'
#' @description
#' If fixedscales is FALSE (or 0) all images will have separate scale.
#' If fixedscales is TRUE (or 1) each x[[i]] will have a common scale.
#' If fixedscales is 2 all images will have common scale.
#'
#' If more than one envelope has been calculated (corresponding to several coverage/alpha),
#' only the largest one is plotted.
#' @inheritParams plot.global_envelope2d
#' @param fixedscales 0, 1 or 2. See details.
#' @param labels A character vector of suitable length giving the labels for the separate plots.
#' Default exists. This parameter allows replacing the default.
#' @importFrom ggplot2 facet_grid facet_wrap ggtitle
#' @importFrom ggplot2 label_value theme element_blank
#' @importFrom gridExtra grid.arrange
#' @export
#' @examples
#' data("abide_9002_23")
#' iset <- subset(abide_9002_23[['curve_set']], 1:50)
#' factors <- abide_9002_23[['factors']][1:50,]
#' \dontshow{
#' # Cut the data to reduce time
#' iset$r <- iset$r[1:7,]
#' iset$funcs <- iset$funcs[1:7, ]
#' }
#' res <- graph.flm(nsim = 19, # Increase nsim for serious analysis!
#'   formula.full = Y ~ Group + Sex + Age,
#'   formula.reduced = Y ~ Sex + Age,
#'   curve_sets = list(Y=iset), factors = factors,
#'   contrasts = FALSE, GET.args = list(type="area"))
#' plot(res)
#' plot(res, what=c("obs", "hi"))
#'
#' plot(res, what=c("hi", "lo"), fixedscales=1)
#'
#' plot(res, what=c("obs", "lo", "hi"), fixedscales=FALSE)
#'
#' if(requireNamespace("gridExtra", quietly=TRUE)) {
#'   # Edit style of "fixedscales = 2" plots
#'   plot(res, what=c("obs", "hi")) + ggplot2::theme_minimal()
#'   plot(res, what=c("obs", "hi")) + ggplot2::theme_bw()
#'
#'   # Edit style (e.g. theme) of "fixedscales = 1 or 0" plots
#'   gs <- lapply(res, function(x, what) { plot(x, what=what) +
#'      ggplot2::ggtitle("") }, what=c("obs", "hi"))
#'   gridExtra::grid.arrange(grobs=gs, ncol=1, top="My main")
#'
#'   gs <- outer(res, c("obs", "hi"), FUN=Vectorize(function(res, what)
#'     list(plot(res, what=what) + ggplot2::ggtitle("") +
#'       ggplot2::theme(axis.ticks=ggplot2::element_blank(),
#'       axis.text=ggplot2::element_blank(), axis.title=ggplot2::element_blank()))))
#'   gridExtra::grid.arrange(grobs=t(gs))
#' }
plot.combined_global_envelope2d <- function(x, fixedscales = 2, labels,
                                            what=c("obs.sign", "obs", "lo", "hi", "lo.sign", "hi.sign"),
                                            sign.type = c("circles", "contour", "col"),
                                            sign.col = c("blue", "red"), transparency = 155/255,
                                            digits = 3, ...) {
  x <- as_simple_envelope(x)
  if(missing(what)) what <- "obs.sign"
  else what <- match.arg(what, several.ok = TRUE)
  what <- checkarg_envelope2d_what(x[[1]], what)
  main <- env_main_default(x, digits=digits)
  if(missing('labels')) labels <- default_labels(x, labels)
  names(x) <- labels
  sign.type <- match.arg(sign.type)

  if(fixedscales==2) {
    g <- plot_combined_global_envelope2d_fixedscales(x, what, sign.col, transparency, sign.type)
    g + ggtitle(main)
  } else if(fixedscales==1) {
    gs <- lapply(seq_along(x), function(i) {
      env <- x[[i]]
      env$main <- names(x)[i]
      g <- plot_global_envelope2d_fixedscales(env, what, sign.col, transparency, sign.type)
      g + facet_grid(main~label)
    })
    grid.arrange(grobs=gs, ncol=1, top=main)
  } else if(fixedscales==0) {
    f <- function(env, main, what) {
      env$main <- main
      g <- plot_global_envelope2d_fixedscales(env, what, sign.col, transparency, sign.type)
      g + facet_wrap(c("main", "label"), labeller=function(l) label_value(l, multi_line = FALSE))
      g + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
    }
    gs <- outer(what, seq_along(x), FUN=Vectorize(function(w, i) {
      list(f(x[[i]], names(x)[i], w))
    }))
    grid.arrange(grobs=gs, top=main, nrow=ncol(gs))
  } else { stop("Invalid fixedscales argument.") }
}
