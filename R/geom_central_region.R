#' @importFrom stats ave
create_cr <- function(x, y, curveid, coverage=0.50, type="erl") {
  y <- y[order(curveid, x)]
  curve_len <- stats::ave(curveid, curveid, FUN = length)
  if (!all(curve_len == curve_len[1])) {
    stop("All curves must have the same argument values (after removing non-finite values).") # ggplot automatically removes infinities from aesthetics?
  }

  obs <- matrix(y, nrow = curve_len[1])

  if(IsNfuncTooSmall(Nfunc=ncol(obs), alpha=1-coverage)) {
    return(NULL)
    #return(data.frame(x=x, y=y, curveid=curveid))
  }

  cset <- create_curve_set(list(obs=obs))
  cr <- central_region(cset, coverage=coverage, type=type)
  cr <- data.frame(cr)
  #cr$x <- cr$r
  #cr$r <- NULL
  cr$x <- x[curveid==curveid[[1]]]
  # Set ymax and ymin as the outermost envelopes, because
  # apparently ylim is deduced from y, ymax, ymin
  level <- sort(1-coverage)
  lonames <- env_loname(level)
  hinames <- env_hiname(level)
  cr$ymax <- cr[[hinames[1]]]
  cr$ymin <- cr[[lonames[1]]]
  cr
}

#' @rdname GeomCentralRegion
#' @importFrom ggplot2 ggproto Stat
#' @export
StatCentralRegion <- ggproto("StatCentralRegion", Stat,
  compute_group = function(data, scales,
                           coverage = 0.50,
                           type = "erl") {
    cols_to_keep <- setdiff(names(data), c("x", "y", "curveid"))
    #rows_to_drop <- data$curveid %in% data$curveid[!is.finite(data$y)]
    #if(any(rows_to_drop)) warning("In StatCentralRegion: Dropping ", sum(rows_to_drop), " rows that are part of curves with non-finite values.")
    #data <- data[!rows_to_drop,]
    d <- create_cr(data$x, data$y, data$curveid, coverage, type)
    if (is.null(d)) {
      d <- data
      d$curveid <- as.numeric(as.factor(d$curveid))
      # position_dodge with geom_central_region needs original group given by ggplot
      # The lines version just needs to keep group integer and avoid collisions
      d$group <- -(100L*d$group + d$curveid)
      d$region <- FALSE
    } else {
      d <- cbind(d, unclass(data[1, cols_to_keep]))
      d$y <- d$central # y needs to be defined to avoid warning of "dropped aesthetic".
      # y needs to be a valid number of avoid "Discrete value"
      # y needs to be withing ymax and ymin to not affect their calculation
      d$region <- TRUE
    }
    d
  },
  required_aes = c("x", "y", "curveid"),
  dropped_aes = c("curveid")
)

#' Central region plot
#'
#' \code{ggproto} objects for central region plot. Not to be used directly.
#'
#' @seealso \code{\link{geom_central_region}}
#' @importFrom grid gList
#' @importFrom ggplot2 Geom GeomRibbon GeomLine
#' @export
GeomCentralRegion <- ggproto("GeomCentralRegion", Geom,
  draw_panel = function(data, panel_params, coord,
                        na.rm = FALSE,
                        coverage = 0.50,
                        type = "erl",
                        filled=TRUE,
                        drawcenterline=TRUE) {

    dataribbon <- data[data$region,]
    p1 <- list()
    if(nrow(dataribbon) > 0) {
      p1 <- NULL
      df <- dataribbon
      if(filled) {
        p2 <- GeomRibbon$draw_panel(data = df,
                              panel_params = panel_params,
                              coord = coord,
                              na.rm = na.rm)
        p1 <- c(p1, list(p2))
      } else {
        df$y <- df$ymax
        p2 <- GeomLine$draw_panel(data = df,
                              panel_params = panel_params,
                              coord = coord,
                              na.rm = na.rm)
        p1 <- c(p1, list(p2))
        df$y <- df$ymin
        p2 <- GeomLine$draw_panel(data = df,
                            panel_params = panel_params,
                            coord = coord,
                            na.rm = na.rm)
        p1 <- c(p1, list(p2))
      }
      if(drawcenterline) {
        df <- dataribbon
        df$y <- df$central
        df$linetype <- 2L
        p1 <- c(p1, list(GeomLine$draw_panel(data = df,
                                             panel_params = panel_params,
                                             coord = coord,
                                             na.rm = na.rm)))
      }
    }
    datalines <- data[!data$region,]
    p2 <- NULL
    # GeomLine will give a warning if the group variable is zero length
    if(nrow(datalines) > 0) {
      datalines$group <- paste(datalines$group, datalines$curveid, sep="_")
      p2 <- GeomLine$draw_panel(
        data = datalines,
        panel_params = panel_params,
        coord = coord,
        na.rm = na.rm
      )
    }
    do.call(grid::gList, c(p1, list(p2)))
  },
  required_aes = c("x"),
  ## TODO: y ymax ymin
  default_aes = aes(
    colour = "black",
    fill = "grey59",
    linewidth = 0.5,
    linetype = 1L,
    alpha = NA
  )
)

#' @importFrom grid gList
#' @importFrom ggplot2 Geom GeomRibbon GeomLine
#' @importFrom grDevices grey.colors
#' @rdname GeomCentralRegion
#' @export
GeomCentralRegionMulti <- ggproto("GeomCentralRegionMulti", Geom,
  draw_panel = function(data, panel_params, coord,
                        na.rm = FALSE,
                        coverage = c(0.95, 0.9, 0.7, 0.5),
                        type = "erl",
                        filled=TRUE,
                        drawcenterline=TRUE,
                        colours=grey.colors(length(coverage), start=0.9, end=0.5)) {

    level <- sort(1-coverage)
    dataribbon <- data[data$region,]
    p1 <- list()
    if(nrow(dataribbon) > 0) {
      lonames <- env_loname(level)
      hinames <- env_hiname(level)
      p1 <- NULL
      for(i in 1:length(level)) {
        df <- dataribbon
        df$ymax <- df[[hinames[i]]]
        df$ymin <- df[[lonames[i]]]
        if(filled) {
          df$fill <- colours[i]
          df$colour <- NA # Hide borders
          p2 <- GeomRibbon$draw_panel(data = df,
                                panel_params = panel_params,
                                coord = coord,
                                na.rm = na.rm)
          p1 <- c(p1, list(p2))
        } else {
          df$colour <- colours[i]
          df$y <- df$ymax
          p2 <- GeomLine$draw_panel(data = df,
                                panel_params = panel_params,
                                coord = coord,
                                na.rm = na.rm)
          p1 <- c(p1, list(p2))
          df$y <- df$ymin
          p2 <- GeomLine$draw_panel(data = df,
                              panel_params = panel_params,
                              coord = coord,
                              na.rm = na.rm)
          p1 <- c(p1, list(p2))
        }
      }
      if(drawcenterline) {
        df <- dataribbon
        df$y <- df$central
        df$colour <- "black"
        df$linetype <- 2L
        p1 <- c(p1, list(GeomLine$draw_panel(data = df,
                                             panel_params = panel_params,
                                             coord = coord,
                                             na.rm = na.rm)))
      }
    }
    datalines <- data[!data$region,]
    p2 <- NULL
    # GeomLine will give a warning if the group variable is zero length
    if(nrow(datalines) > 0) {
      datalines$group <- paste(datalines$group, datalines$curveid, sep="_")
      datalines$colour <- "black"
      p2 <- GeomLine$draw_panel(
        data = datalines,
        panel_params = panel_params,
        coord = coord,
        na.rm = na.rm
      )
    }
    do.call(grid::gList, c(p1, list(p2)))
  },
  required_aes = c("x"),
  ## TODO: y ymax ymin
  default_aes = aes(
    # colour = "black",
    # fill = "grey59",
    linewidth = 0.5,
    linetype = 1L,
    alpha = NA
  )
)


#' Central region plot
#'
#' \code{geom_central_region} and \code{stat_central_region} can be used to compute
#' and plot \code{central_region} from data arranged in a \code{data.frame}.
#'
#'
#' Plots of central regions (global envelopes) with the specified \code{coverage}
#' and \code{type} (see \code{\link{central_region}}).
#' When splitting the set of functions to groups by aesthetics or facets, see
#' examples, the central regions are constructed separately for each group,
#' each having the specified \code{coverage}.
#'
#' If Nfunc*(1-coverage) < 1, where Nfunc is the number of functions/curves,
#' the curves are plotted instead of any region.
#'
#' @section Aesthetics:
#' \code{geom_central_region} requires \code{x}, \code{y} and \code{curveid}.
#' Additionally \code{geom_central_region} uses the same aesthetics as
#' \code{\link{geom_ribbon}} if \code{filled==TRUE} and
#' \code{\link{geom_line}} otherwise.
#' For multiple coverage values additional aesthetics are not currently supported.
#'
#' @section Computed variables:
#' \code{stat_central_region} computes
#' \code{after_stat(ymax)} and \code{after_stat(ymin)} for the high and low value of the central region.
#'
#' For multiple coverages the variables use the same names as \code{\link{central_region}},
#' i.e. \code{hi.95} and \code{lo.95} for the region with 95% coverage.
#'
#' @seealso \code{\link{central_region}} for the basic computation and,
#' \code{\link{geom_ribbon}} for the default base geom.
#'
#' @inheritParams ggplot2::geom_ribbon
#' @inheritParams central_region
#' @param type The options and details for \code{type} are given in \code{\link{central_region}}.
#' @param filled Boolean. Should the ribbon be filled?
#' @param drawcenterline Boolean. Should the center line be drawn?
#' @param colours Colours for different coverage levels
#' @examples
#' require("ggplot2")
#' ## Generate some data
#' #------------------------------------------------------
#' # Simulate regression data according to the cubic model
#' # f(x) = 0.8x - 1.8x^2 + 1.05x^3 for x in [0,1]
#' par <- c(0,0.8,-1.8,1.05) # Parameters of the true polynomial model
#' res <- 100 # Resolution
#' x <- seq(0, 1, by=1/res); x2=x^2; x3=x^3;
#'
#' f <- par[1] + par[2]*x + par[3]*x^2 + par[4]*x^3 # The true function
#' d <- f + rnorm(length(x), 0, 0.04) # Data
#'
#' # Estimate polynomial regression model
#' reg <- lm(d ~ x + x2 + x3)
#' ftheta <- reg$fitted.values
#' resid0 <- reg$residuals
#'
#' # Bootstrap regression
#' B <- 200 # Number of bootstrap samples
#' df <- NULL
#' for(i in 1:B) {
#'   u <- sample(resid0, size=length(resid0), replace=TRUE)
#'   reg1 <- lm((ftheta+u) ~ x + x2 + x3)
#'   df <- rbind(df, data.frame(y=reg1$fitted.values, x=x, i=i,
#'     g=ifelse(i<14, "A", "B"), g2=ifelse(i<100, "A", "B")))
#' }
#'
#' ggplot(df) + geom_line(aes(x, y, group=i))
#' ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=0.50)
#' ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=0.50, filled=FALSE)
#' # Central regions for two groups as specified by 'g2'
#' ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i, col=g2), coverage=0.90, filled=FALSE)
#' ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=0.90) + facet_wrap(vars(g2))
#' \dontshow{
#' # If nr. of functions < 20, then the functions are drawn; otherwise the 100*coverage% central region
#' ggplot(df[df$i < 10,]) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=0.90)
#' # Central regions for two groups split by 'g'; <20 functions in the first group
#' ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i, col=g, fill=g), coverage=0.90)
#' ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i, col=g), coverage=0.90, filled=FALSE)
#' }
#'
#' # Central regions with multiple coverage levels
#' ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=c(0.2,0.4,0.6)) +
#'   theme_minimal()
#' ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i), coverage=seq(0.1, 0.9, length=20),
#'   colours=rainbow(20))
#'
#' # Colors for multiregions are not supported
#' ggplot(df) + geom_central_region(aes(x=x, y=y+0.1*(g2=="B"),
#'   curveid=i, col=as.factor(g2)), coverage=c(0.05, 0.2,0.4,0.6))
#' ggplot(df) + geom_central_region(aes(x=x, y=y, curveid=i),
#'   coverage=c(0.05, 0.2,0.4,0.6)) + facet_wrap(vars(g2))
#'
#' # Using stat_central_region with geom_linerange and geom_rect
#' ggplot(df) +
#'   geom_linerange(aes(curveid=i, x=x, y=y, ymax=after_stat(ymax), ymin=after_stat(ymin),
#'                group=g2, col=factor(g2)),
#'                stat="central_region", coverage = 0.90, position=position_dodge(0.01))
#' ggplot(within(df, {x = x+0.004*(g2=="B")})) +
#'   geom_rect(aes(curveid=i, x=x, y=y, xmax=after_stat(x), xmin=after_stat(x+0.004),
#'               ymax=after_stat(ymax), ymin=after_stat(ymin), group=g2, fill=factor(g2)),
#'               stat="central_region", coverage = 0.90)
#' # Non-finite values are not supported
#' ggplot(within(df, {y = ifelse(runif(length(y)) < 0.001, Inf, y)})) +
#'   geom_central_region(aes(x=x, y=y, curveid=i))
#' @export
geom_central_region <- function(mapping = NULL, data = NULL, stat = "CentralRegion",
                        position = "identity", ...,
                        coverage = 0.50,
                        type = "erl",
                        filled = TRUE,
                        drawcenterline = TRUE,
                        colours = grey.colors(length(coverage), start=0.9, end=0.5),
                        na.rm = FALSE,
                        show.legend = NA, inherit.aes = TRUE) {
  if(length(coverage)==1 && !missing(colours)) warning("'colors' only applies to multiple coverage levels.")

  dots <- list(...)
  # colours could be a function that makes a gradient by multiplying the colour given by aesthetic
  if(length(coverage) != 1) dots <- c(dots, list(colours=colours))
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = if(length(coverage)==1) {GeomCentralRegion} else {GeomCentralRegionMulti},
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = c(list(
      type = type,
      coverage = coverage,
      filled = filled,
      drawcenterline = drawcenterline,
      na.rm = na.rm
    ), dots)
  )
}

#' @rdname geom_central_region
#' @inheritParams ggplot2::stat_identity
#' @inheritParams geom_central_region
#' @export
stat_central_region <- function(mapping = NULL, data = NULL,
                                position = "identity", ...,
                                coverage = 0.50,
                                type = "erl",
                                na.rm = FALSE,
                                show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatCentralRegion,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      type = type,
      coverage = coverage,
      na.rm = na.rm,
      ...
    )
  )
}

