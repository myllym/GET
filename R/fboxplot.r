#' Functional boxplot
#'
#' Functional boxplot based on central region computed by a specified measure.
#' The options of the measures can be found in \code{\link{central_region}}.
#' @inheritParams central_region
#' @param factor The constant factor to inflate the central region to produce a functional boxplot and
#' determine fences for outliers. Default is 1.5 as in a classical boxplot.
#' @param ... Additional parameters to be passed to \code{\link{central_region}},
#' which is responsible for calculating the central region (global envelope) on which
#' the functional boxplot is based.
#' @export
#' @examples
#' if(requireNamespace("fda", quietly=TRUE)) {
#'   years <- paste(1:18)
#'   curves <- fda::growth[['hgtf']][years,]
#'   # Heights
#'   cset1 <- create_curve_set(list(r = as.numeric(years),
#'                                  obs = curves))
#'   plot(cset1, ylab="Height")
#'   bp <- fBoxplot(cset1, coverage=0.50, type="area", factor=1)
#'   plot(bp)
#'   # Use fbplot from fda:
#'   area_depth <- forder(cset1, measure = "area")
#'   fda::fbplot(fit = cset1$funcs, x=cset1$r, depth=area_depth, factor=1)
#'
#'   # Considering simultaneously heights and height differences
#'   cset2 <- create_curve_set(list(r = as.numeric(years[-1]),
#'              obs = curves[-1,] - curves[-nrow(curves),]))
#'   csets <- list(Height = cset1, Change = cset2)
#'   res <- fBoxplot(csets, type = 'area', factor = 1.5)
#'   plot(res, xlab = "Age (years)", ylab = "")
#' }
fBoxplot <- function(curve_sets, factor = 1.5, ...) {
  if(length(class(curve_sets)) == 1 && class(curve_sets) == "list") {
    if(!all(sapply(curve_sets, FUN=curve_set_is1d)))
      stop("r in curve_sets should be vectors.")
    if(any(sapply(curve_sets, FUN=curve_set_is1obs)))
      warning("One of the curves is observed, other simulated.\n Do you want to make a functional boxplot of all these?")
  }
  else { if(!curve_set_is1d(curve_sets)) stop("r in curve_sets should be a vector.") }
  res <- central_region(curve_sets, alternative = "two.sided", ...)
  if(inherits(res, "combined_global_envelope")) {
    dist <- factor * (attr(res, "level2_ge")$hi - attr(res, "level2_ge")$lo)
    attr(res, "level2_ge")$whisker.lo <- attr(res, "level2_ge")$lo - dist
    attr(res, "level2_ge")$whisker.hi <- attr(res, "level2_ge")$hi + dist
    attr(attr(res, "level2_ge"), "method") <- "Functional boxplot"
    class(attr(res, "level2_ge")) <- c("fboxplot", class(attr(res, "level2_ge")))
    outliers_id <- NULL
    for(i in 1:length(res)) {
      dist <- factor * (res[[i]]$hi - res[[i]]$lo)
      res[[i]]$whisker.lo <- res[[i]]$lo - dist
      res[[i]]$whisker.hi <- res[[i]]$hi + dist
      attr(res[[i]], "method") <- "Functional boxplot"
      class(res[[i]]) <- c("fboxplot", class(res[[i]]))
      # Outliers
      funcs <- curve_set_funcs(curve_sets[[i]])
      for(j in 1:ncol(funcs)) {
        if(any(funcs[,j] < res[[i]]$whisker.lo | funcs[,j] > res[[i]]$whisker.hi))
          outliers_id <- c(outliers_id, j)
      }
    }
    outliers_id <- unique(outliers_id)
    if(!is.null(outliers_id)) {
      outliers <- lapply(curve_sets, FUN = function(x) { curve_set_funcs(x)[,outliers_id] })
      if(is.vector(outliers[[1]])) outliers <- lapply(outliers, FUN = function(x) matrix(x, ncol=1))
      if(is.null(colnames(outliers[[1]])))
        for(i in 1:length(outliers)) {
          n <- colnames(curve_set_funcs(curve_sets[[i]]))
          if(!is.null(n))
            colnames(outliers[[i]]) <- n[outliers_id]
          else
            colnames(outliers[[i]]) <- outliers_id
        }
    }
    else
      outliers <- NULL
    attr(res, "outliers") <- outliers
    attr(res, "factor") <- factor
    attr(res, "method") <- "Combined functional boxplot"
    attr(res, "call") <- match.call()
    class(res) <- c("combined_fboxplot", class(res))
  }
  else {
    dist <- factor * (res$hi - res$lo)
    res$whisker.lo <-  res$lo - dist
    res$whisker.hi <- res$hi + dist
    # Outliers
    outliers_id <- NULL
    funcs <- curve_set_funcs(curve_sets)
    for(j in 1:ncol(funcs)) {
      if(any(funcs[,j] < res$whisker.lo | funcs[,j] > res$whisker.hi))
        outliers_id <- c(outliers_id, j)
    }
    if(!is.null(outliers_id)) {
      outliers <- funcs[,outliers_id]
      if(is.null(colnames(outliers)))
        colnames(outliers) <- outliers_id
    }
    else
      outliers <- NULL
    attr(res, "outliers") <- outliers
    attr(res, "outpoint") <- outliers_id
    attr(res, "method") <- "Functional boxplot"
    attr(res, "factor") <- factor
    attr(res, "call") <- match.call()
    class(res) <- c("fboxplot", class(res))
  }
  res
}

#' Plot method for the class 'fboxplot'
#'
#' @param x an 'fboxplot' object
#' @inheritParams plot.global_envelope
#' @param plot_outliers Logical. If TRUE, then the functions outside the functional boxplot are drawn,
#'  with legend if \code{legend} is TRUE.
#' @param legend Logical. See \code{plot_outliers.}
#' @param ... Additional arguments to be passed to \code{\link{plot.global_envelope}}.
#' @export
plot.fboxplot <- function(x, main, xlab, ylab, digits = 3,
                          legend = TRUE, plot_outliers = TRUE, ...) {
  if(missing(main)) main <- env_main_default(x, digits=digits)
  d <- plotdefaultlabs(x, xlab, ylab)
  fboxplot_ggplot(x, main=main, xlab=d$xlab, ylab=d$ylab, plot_outliers=plot_outliers, legend=legend)
}

#' Plot method for the class 'combined_fboxplot'
#'
#' @param x an 'combined_fboxplot' object
#' @inheritParams plot.combined_global_envelope
#' @inheritParams plot.fboxplot
#' @param ... Additional arguments to be passed to \code{\link{plot.combined_global_envelope}}.
#'
#' @export
plot.combined_fboxplot <- function(x, main, xlab, ylab, labels, scales = "free",
                                   digits = 3, ncol = 2 + 1*(length(x)==3),
                                   plot_outliers = TRUE, legend = TRUE, ...) {
  alt <- get_alternative(x[[1]])
  if(missing(main)) main <- env_main_default(attr(x, "level2_ge"), digits=digits, alternative=alt)
  if(missing(labels)) labels <- default_labels(x, labels)
  d <- plotdefaultlabs(attr(x, "level2_ge"), xlab, ylab)
  # Combined plot, level 1 plots
  fboxplot_combined_ggplot(x, main=main, xlab=d$xlab, ylab=d$ylab,
                           labels=labels, scales=scales,
                           max_ncols_of_plots=ncol, legend=legend)
}
