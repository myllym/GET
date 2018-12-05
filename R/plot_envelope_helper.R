# A helper function for avoiding repetition of code (in envelope tests).
pick_attributes <- function(curve_set, alternative="two.sided") {
    # saving for attributes / plotting purposes
    lo.name <- "lower critical boundary for %s"
    hi.name <- "upper critical boundary for %s"
    switch(alternative,
            two.sided = {},
            less = {
                hi.name <- "infinite upper boundary"
            },
            greater = {
                lo.name <- "infinite lower boundary"
            })
    if(inherits(curve_set, 'envelope')) {
        fname <- attr(curve_set, "fname")
        labl <- attr(curve_set, "labl")
        desc <- attr(curve_set, "desc")
        desc[4] <- lo.name
        desc[5] <- hi.name
        ylab <- attr(curve_set, "ylab")
        einfo <- attr(curve_set, "einfo")
    }
    else {
        fname <- "T"
        if(is.vector(curve_set[['obs']])) {
          labl <- c("r", "T[obs](r)", "T[0](r)", "T[lo](r)", "T[hi](r)")
          desc <- c("distance argument r",
                    "observed value of %s for data pattern",
                    "central curve under the null hypothesis",
                    lo.name, hi.name)
        }
        else {
          labl <- c("r", "T[0](r)", "T[lo](r)", "T[hi](r)")
          desc <- c("distance argument r",
                    "central curve under the null hypothesis",
                    lo.name, hi.name)
        }
        ylab <- "T(r)"
        einfo <- list(Yname="curve_set[['obs']]", valname="", csr=FALSE, crs.theo=FALSE,
                      use.theory=with(curve_set, exists("theo")), pois=FALSE,
                      simtype="other",
                      global=TRUE, ginterval=range(curve_set[['r']]),
                      VARIANCE=FALSE, alternative=alternative)
    }
    list(fname=fname, labl=labl, desc=desc, ylab=ylab, einfo=einfo)
}

# A helper function to check whether the xaxis needs to be reticked with new values due to
# combined tests. Called also from plot.global_envelope for checking if fv plot style is available.
retick_xaxis <- function(x) {
  if(class(x)[1] != "list") x <- list(x)
  if(any(unlist(lapply(x, FUN=function(x) { !("global_envelope" %in% class(x) | "fboxplot" %in% class(x) | "curve_set" %in% class(x)) }))))
    stop("x should consist of global_envelope objects.\n")
  r_values_ls <- lapply(x, FUN=function(x) x$r)
  r_values <- do.call(c, r_values_ls, quote=FALSE)
  nr <- length(r_values)
  list(retick_xaxis = !(length(x) == 1 & all(r_values[-1] - r_values[-nr] > 0)),
       r_values_ls = r_values_ls, r_values = r_values)
}

# Define breaking r values and labels on x-axis for plotting several
# global envelopes or curve sets jointly.
# @param x A global_envelope object of a list of global_envelope objects.
# Also curve_set objects allowed with a restricted use.
# @param nticks Number of ticks per a sub test.
combined_global_envelope_rhelper <- function(x, nticks = 5) {
  if(class(x)[1] != "list") x <- list(x)
  retick <- retick_xaxis(x)
  r_values_ls <- retick$r_values_ls
  r_values <- retick$r_values
  nr <- length(r_values)
  if(!retick$retick_xaxis) {
      new_r_values <- NULL
      r_break_values <- NULL
      loc_break_values <- NULL
      r_values_newstart_id <- NULL
  }
  else {
    if(length(x) > 1 & any(unlist(lapply(r_values_ls, function(x) { !all(x[-1] - x[-length(x)] > 0) })))) {
      warning(paste("Something strange. The r values are not increasing in a", class(x[[1]])[1], "object.\n", sep=""))
    }
    new_r_values <- 1:nr # to be used in plotting
    # Define where the functions start when they are put all together
    if(length(x) == 1) { # Find from the r-values
      r_values_newstart_id <- which(!(r_values[1:(nr-1)] < r_values[2:nr])) + 1
    }
    else { # Define directly from the r_values_ls
      r_values_newstart_id <- NULL
      r_values_newstart_id[1] <- length(r_values_ls[[1]]) + 1
      if(length(r_values_ls) > 2) {
        for(i in 2:(length(r_values_ls)-1))
          r_values_newstart_id <- c(r_values_newstart_id, r_values_newstart_id[i-1] + length(r_values_ls[[1]]))
      }
    }
    # r-values for labeling ticks
    r_starts <- r_values[c(1, r_values_newstart_id)]
    r_ends <- r_values[c(r_values_newstart_id - 1, nr)]
    r_break_values <- NULL
    # indeces for ticks in the running numbering from 1 to nr
    loc_starts <- (1:nr)[c(1, r_values_newstart_id)]
    loc_ends <- (1:nr)[c(r_values_newstart_id - 1, nr)]
    loc_break_values <- NULL
    nslots <- length(r_starts) # number of combined tests/slots
    for(i in 1:(nslots-1)) {
      r_break_values <- c(r_break_values, seq(r_starts[i], r_ends[i], length=nticks)[1:(nticks-1)])
      loc_break_values <- c(loc_break_values, seq(loc_starts[i], loc_ends[i], length=nticks)[1:(nticks-1)])
    }
    r_break_values <- c(r_break_values, seq(r_starts[nslots], r_ends[nslots], length=nticks))
    loc_break_values <- c(loc_break_values, seq(loc_starts[nslots], loc_ends[nslots], length=nticks))
  }
  if(class(x[[1]])[1] %in% c("global_envelope", "fboxplot")) {
    if(!is.null(x[[1]]$obs))
      x_vec <- data.frame(r = r_values,
                          obs = do.call(c, lapply(x, FUN = function(x) x$obs), quote=FALSE),
                          central = do.call(c, lapply(x, FUN = function(x) x$central), quote=FALSE),
                          lo = do.call(c, lapply(x, FUN = function(x) x$lo), quote=FALSE),
                          hi = do.call(c, lapply(x, FUN = function(x) x$hi), quote=FALSE))
    else
      x_vec <- data.frame(r = r_values,
                          central = do.call(c, lapply(x, FUN = function(x) x$central), quote=FALSE),
                          lo = do.call(c, lapply(x, FUN = function(x) x$lo), quote=FALSE),
                          hi = do.call(c, lapply(x, FUN = function(x) x$hi), quote=FALSE))
  }
  else x_vec <- NULL

  list(x_vec = x_vec,
       retick_xaxis = retick$retick_xaxis,
       new_r_values = new_r_values,
       r_break_values = r_break_values, loc_break_values = loc_break_values,
       r_values_newstart_id = r_values_newstart_id)
}

# An internal GET function for setting the default main for a global envelope plot.
# @param x An 'global_envelope' object.
env_main_default <- function(x, digits=3) {
  if(!is.null(attr(x, "p_interval"))) {
    if(attr(x, "alternative") == "two.sided")
      main <- paste(attr(x, "method"), ": p-interval = (",
                    round(attr(x, "p_interval")[1], digits=digits),", ",
                    round(attr(x, "p_interval")[2], digits=digits), ")", sep="")
    else
      main <- paste(attr(x, "method"), ": p-interval = (",
                    round(attr(x, "p_interval")[1], digits=digits),", ",
                    round(attr(x, "p_interval")[2], digits=digits), ") \n",
                    "Alternative = \"", attr(x, "alternative"), "\"\n", sep="")
  }
  else {
    if(!is.null(attr(x, "p"))) {
      if(attr(x, "alternative") == "two.sided")
        main <- paste(attr(x, "method"), ": p = ", round(attr(x, "p"), digits=digits), sep="")
      else
        main <- paste(attr(x, "method"), ": p = ", round(attr(x, "p"), digits=digits), "\n",
                      "Alternative = \"", attr(x, "alternative"), "\"\n", sep="")
    }
    else {
      switch(class(x)[1],
             global_envelope = {
               main <- paste(100*(1-attr(x, "alpha")), "% central region (", attr(x, "type"), ")", sep="")
             },
             fboxplot = {
               main <- paste("Functional boxplot based on ", 100*(1-attr(x, "alpha")), "% central region (", attr(x, "type"), ")", sep="")
             })
    }
  }
  main
}

# An internal GET function for setting the default ylim for a global envelope plot.
# @param x An 'global_envelope' object or a list of them.
# @param use_ggplot2 TRUE/FALSE, If TRUE, then default ylim are for \code{\link{env_ggplot}}.
# Otherwise the ylim are for \code{\link{env_basic_plot}}.
env_ylim_default <- function(x, use_ggplot2) {
  if(class(x)[1] != "list") x <- list(x)
  if(!use_ggplot2)
    switch(attr(x[[1]], "alternative"),
            two.sided = {
              ylim <- c(do.call(min, lapply(x, function(x) min(x[['obs']],x[['lo']],x[['hi']],x[['central']])), quote=FALSE),
                        do.call(max, lapply(x, function(x) max(x[['obs']],x[['lo']],x[['hi']],x[['central']])), quote=FALSE))
            },
            less = {
              ylim <- c(do.call(min, lapply(x, function(x) min(x[['obs']],x[['lo']],x[['central']])), quote=FALSE),
                        do.call(max, lapply(x, function(x) max(x[['obs']],x[['lo']],x[['central']])), quote=FALSE))
            },
            greater = {
              ylim <- c(do.call(min, lapply(x, function(x) min(x[['obs']],x[['hi']],x[['central']])), quote=FALSE),
                        do.call(max, lapply(x, function(x) max(x[['obs']],x[['hi']],x[['central']])), quote=FALSE))
            })
  else ylim <- NULL
  ylim
}

# An internal GET function for making a dotplot style "global envelope plot".
#
# An internal GET function for making a dotplot style "global envelope plot".
#
# @param x An 'envelope_test' object.
# @param main See \code{\link{plot.default}}.
# @param ylim See \code{\link{plot.default}}.
# @param xlab See \code{\link{plot.default}}.
# @param ylab See \code{\link{plot.default}}.
# @param color_outside Logical. Whether to color the places where the data function goes
# outside the envelope. Currently red color is used.
# @param labels Labels for the tests at x-axis.
# @param add Whether to add the plot to an existing plot (TRUE) or to draw a new plot (FALSE).
# @param arrows.col Color for the doplot arrows. If not given, 1 (black) is used.
# @param curve_sets If provided, then curves going outside the envelope are plotted.
# @param ... Additional parameters to be passed to the function \code{\link{plot}}.
#' @importFrom graphics plot
#' @importFrom graphics arrows
#' @importFrom graphics points
#' @importFrom graphics axis
env_dotplot <- function(x, main, ylim, xlab, ylab, color_outside = TRUE,
                        labels = NULL, add = FALSE, arrows.col, curve_sets = NULL, ...) {
    nr <- length(x[['r']])
    if(is.null(labels)) labels <- paste(round(x[['r']], digits=2))
    if(missing(arrows.col)) arrows.col <- 1
    if(nr > 10) warning("Dotplot style meant for low dimensional test vectors.\n")
    if(!add) graphics::plot(1:nr, x[['central']], main=main, ylim=ylim, xlab=xlab, ylab=ylab, cex=0.5, pch=16, xaxt="n", ...)
    else graphics::points(1:nr, x[['central']], main=main, ylim=ylim, xlab=xlab, ylab=ylab, cex=0.5, pch=16, xaxt="n", ...)
    if(attr(x, "alternative")!="greater")
        graphics::arrows(1:nr, x[['lo']], 1:nr, x[['central']], code = 1, angle = 75, length = .1, col=arrows.col)
    else
        graphics::arrows(1:nr, x[['lo']], 1:nr, x[['central']], code = 1, angle = 75, length = .1, col=grey(0.8))
    if(attr(x, "alternative")!="less")
        graphics::arrows(1:nr, x[['hi']], 1:nr, x[['central']], code = 1, angle = 75, length = .1, col=arrows.col)
    else
        graphics::arrows(1:nr, x[['hi']], 1:nr, x[['central']], code = 1, angle = 75, length = .1, col=grey(0.8))
    graphics::axis(1, 1:nr, labels=labels)
    if(!is.null(x[['obs']])) {
      graphics::points(1:nr, x[['obs']], pch='x')
      if(color_outside) {
        outside <- x[['obs']] < x[['lo']] | x[['obs']] > x[['hi']]
        graphics::points((1:nr)[outside], x[['obs']][outside], pch='x', col="red")
      }
    }
    if(!is.null(curve_sets)) {
      for(i in 1:ncol(curve_sets$obs)) {
        if(any(curve_sets$obs[,i] < x[['lo']] | curve_sets$obs[,i] > x[['hi']])) {
          graphics::points(1:nr, curve_sets$obs[,i], pch='x', col=grey(0.7), type="b")
        }
      }
    }
}


# An internal GET function for making a basic "global envelope plot".
#
# An internal GET function for making a basic "global envelope plot".
#
# @inheritParams env_dotplot
# @param separate_yaxes Logical (default FALSE). By default also the combined envelope plots have
# a common y-axis. If TRUE, then separate y-axes are used for different parts of a combined test.
# @param max_ncols_of_plots If separate_yaxes is TRUE, then max_ncols_of_plots gives the maximum
# number of columns for figures. Default 2.
# @param env.col The color for the envelope lines. Default 1 (black).
# @param nticks The number of ticks on the xaxis, if the xaxis is re-ticked for combined tests.
# @param ... Additional parameters to be passed to the function \code{\link{plot}}.
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics lines
#' @importFrom graphics axis
#' @importFrom graphics abline
env_basic_plot <- function(x, main, ylim, xlab, ylab, color_outside=TRUE,
                           separate_yaxes = FALSE, max_ncols_of_plots = 2, add = FALSE, env.col = 1,
                           nticks = 5, curve_sets = NULL, ...) {
    if(class(x)[1] != "list") x <- list(x)
    # Handle combined tests; correct labels on x-axis
    # a) if x is a list of global_envelope objects
    # b) if x[['r']] contains repeated values (when length(x) == 1)
    rdata <- combined_global_envelope_rhelper(x, nticks=nticks)
    alt <- attr(x[[1]], "alternative")
    x <- rdata$x_vec
    # Plot
    if(!separate_yaxes) {
        if(!rdata$retick_xaxis) {
            if(!add) graphics::plot(x[['r']], x[['central']], main=main, ylim=ylim, xlab=xlab, ylab=ylab,
                                    type="l", lty=3, lwd=2, ...)
            else graphics::lines(x[['r']], x[['central']], lty=3, lwd=2, ...)
        }
        else {
          x[['r']] <- 1:length(x[['r']])
          if(!add) graphics::plot(x[['r']], x[['central']], main=main, ylim=ylim, xlab=xlab, ylab=ylab,
                                  type="l", lty=3, lwd=2, xaxt="n", ...)
          else graphics::lines(x[['r']], x[['central']], lty=3, lwd=2, xaxt="n", ...)
          if(!is.null(curve_sets) && class(curve_sets)[1] == "list") curve_sets <- combine_curve_sets(curve_sets)
        }
        if(alt != "greater") lines(x[['r']], x[['lo']], lty=2, col=env.col) else lines(x[['r']], x[['lo']], lty=2, col=grey(0.8))
        if(alt != "less") lines(x[['r']], x[['hi']], lty=2, col=env.col) else lines(x[['r']], x[['hi']], lty=2, col=grey(0.8))
        if(!is.null(x[['obs']])) {
          lines(x[['r']], x[['obs']], lty=1)
          if(color_outside) {
              outside <- x[['obs']] < x[['lo']] | x[['obs']] > x[['hi']]
              graphics::points(x[['r']][outside], x[['obs']][outside], col="red")
          }
        }
        # curves/outliers
        if(!is.null(curve_sets)) {
          for(i in 1:ncol(curve_sets$obs)) {
            if(any(curve_sets$obs[,i] < x[['lo']] | curve_sets$obs[,i] > x[['hi']]))
              lines(x[['r']], curve_sets$obs[,i], col=grey(0.7))
          }
        }
        if(rdata$retick_xaxis) {
            graphics::axis(1, rdata$loc_break_values, labels=paste(round(rdata$r_break_values, digits=2)))
            graphics::abline(v = rdata$new_r_values[rdata$r_values_newstart_id], lty=3)
        }
    }
    else {
        n_of_plots <- as.integer(1 + length(rdata$r_values_newstart_id))
        ncols_of_plots <- min(n_of_plots, max_ncols_of_plots)
        nrows_of_plots <- ceiling(n_of_plots / ncols_of_plots)
        if(length(xlab)!=n_of_plots) {
            if(length(xlab)==1) xlab <- rep(xlab, times=n_of_plots)
            else warning("The length of the vector xlab is unreasonable.\n")
        }
        if(length(ylab)!=n_of_plots) {
            if(length(ylab)==1) ylab <- rep(ylab, times=n_of_plots)
            else warning("The length of the vector ylab is unreasonable.\n")
        }
        graphics::par(mfrow=c(nrows_of_plots, ncols_of_plots))
        cat("Note: \"main\" and \"ylim\" ignored as separate plots are produced.\n")
        tmp_indeces <- c(1, rdata$r_values_newstart_id, length(rdata$new_r_values)+1)
        if(!is.null(curve_sets)) {
          curve_sets <- combine_curve_sets(curve_sets)
        }
        for(i in 1:n_of_plots) {
            ylim <- c(min(x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                          x[['lo']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                          x[['hi']][tmp_indeces[i]:(tmp_indeces[i+1]-1)]),
                      max(x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                          x[['lo']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                          x[['hi']][tmp_indeces[i]:(tmp_indeces[i+1]-1)]))
            if(!add)
              graphics::plot(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             x[['central']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             main="", xlab=xlab[i], ylab=ylab[i],
                             type="l", lty=3, lwd=2, ylim=ylim, ...)
            else
              graphics::lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             x[['central']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             lty=3, lwd=2, ...)
            if(alt != "greater")
                graphics::lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], x[['lo']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2)
            else
                graphics::lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], x[['lo']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2, col=grey(0.8))
            if(alt != "less")
                graphics::lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], x[['hi']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2)
            else
                graphics::lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], x[['hi']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2, col=grey(0.8))
            if(!is.null(x[['obs']])) {
              graphics::lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=1)
              if(color_outside) {
                  outside <- x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)] < x[['lo']][tmp_indeces[i]:(tmp_indeces[i+1]-1)] | x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)] > x[['hi']][tmp_indeces[i]:(tmp_indeces[i+1]-1)]
                  graphics::points(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)][outside], x[['obs']][tmp_indeces[i]:(tmp_indeces[i+1]-1)][outside], col="red")
              }
            }
            # curves/outliers
            if(!is.null(curve_sets)) {
              for(j in 1:ncol(curve_sets$obs)) {
                if(any(curve_sets$obs[,j] < x[['lo']] | curve_sets$obs[,j] > x[['hi']]))
                  lines(x[['r']][tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                        curve_sets$obs[tmp_indeces[i]:(tmp_indeces[i+1]-1),j], col=grey(0.7))
              }
            }
        }
    }
}


# An internal GET function for making a ggplot2 style "global envelope plot".
#
# An internal GET function for making a ggplot2 style "global envelope plot".
#
# @param x An 'global_envelope' object.
# @param base_size Base font size, to be passed to theme style when \code{plot_type="ggplot2"}.
# @param main See \code{\link{plot.default}}.
# @param ylim See \code{\link{plot.default}}.
# @param xlab See \code{\link{plot.default}}.
# @param ylab See \code{\link{plot.default}}.
# @param separate_yaxes Logical (default FALSE). By default also the combined envelope plots have
# a common y-axis. If TRUE, then separate y-axes are used for different parts of a combined test.
# @param labels Labels for the separate plots. Ignored if separate_yaxes is FALSE.
# @param max_ncols_of_plots If separate_yaxes is TRUE, then max_ncols_of_plots gives the maximum
# number of columns for figures. Default 2.
# @param labels Labels for components of the combined tests.
# @param nticks The number of ticks on the xaxis, if the xaxis is re-ticked for combined tests.
# @param curve_sets If provided, then curves going outside the envelope are plotted.
# @param x2 Another 'global_envelope' object, which is plotted within x, i.e. x2 is assumed to be narrower
# of the two envelopes.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 aes_
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_linetype_manual
#' @importFrom ggplot2 scale_size_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 labs
env_ggplot <- function(x, base_size, main, ylim, xlab, ylab,
                       separate_yaxes = TRUE, max_ncols_of_plots = 2,
                       labels = NULL, nticks = 5, curve_sets = NULL, x2 = NULL) {
    if(class(x)[1] != "list") x <- list(x)
    if(!is.null(x2)) {
      if(class(x2)[1] != "list") x2 <- list(x2)
      if(length(x) != length(x2)) {
        warning("Unsuitable x2. Setting it to NULL.\n")
        x2 <- NULL
      }
      else {
        for(i in 1:length(x)) {
          if(!all(x[[i]][['r']] == x2[[i]][['r']])) stop("The two envelopes are for different r-values.\n")
          if(!all(x[[i]][['central']] == x2[[i]][['central']])) warning("The two envelopes have different central functions!\n")
        }
      }
      rdata <- combined_global_envelope_rhelper(x2, nticks=nticks)
      x2 <- rdata$x_vec
    }
    # Handle combined tests; correct labels on x-axis
    # a) if x is a list of global_envelope objects
    # b) if x[['r']] contains repeated values (when length(x) == 1)
    rdata <- combined_global_envelope_rhelper(x, nticks=nticks)
    alt <- attr(x[[1]], "alternative")
    x <- rdata$x_vec

    linetype.values <- c('dashed', 'solid')
    size.values <- c(0.2, 0.2)

    if(!is.null(curve_sets)) {
      if(class(curve_sets)[1] == "list") curve_sets <- combine_curve_sets(curve_sets)
      counter <- 0
      outliers <- NULL
      for(j in 1:ncol(curve_sets$obs)) {
        if(any(curve_sets$obs[,j] < x[['lo']] | curve_sets$obs[,j] > x[['hi']])) {
          outliers <- c(outliers, curve_sets$obs[,j])
          counter <- counter + 1
        }
      }
    }

    if(!separate_yaxes | is.null(rdata$r_values_newstart_id)) {
        if(rdata$retick_xaxis) x[['r']] <- 1:length(x[['r']])
        if(is.null(x[['obs']])) {
          df <- data.frame(r = x[['r']],
                           curves = x[['central']],
                           type = factor("Central function", levels = "Central function"),
                           lower = x[['lo']],
                           upper = x[['hi']],
                           main = main)
          if(!is.null(x2)) {
            df$lower2 <- x2[['lo']]
            df$upper2 <- x2[['hi']]
          }
        }
        else {
          df <- data.frame(r = rep(x[['r']], times=2),
                           curves = c(x[['central']], x[['obs']]),
                           type = factor(rep(c("Central function", "Data function"),
                                             each=length(x[['r']])),
                                         levels=c("Data function", "Central function")),
                           lower = rep(x[['lo']], times=2),
                           upper = rep(x[['hi']], times=2),
                           main = main)
          if(!is.null(x2)) {
            df$lower2 <- rep(x2[['lo']], times=2)
            df$upper2 <- rep(x2[['hi']], times=2)
          }
        }
        if(is.null(x2)) {
          p <- ( ggplot2::ggplot()
                 + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower, ymax = ~upper),
                                        fill = 'grey59', alpha = 1)
          )
        }
        else {
          p <- ( ggplot2::ggplot()
                 + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower, ymax = ~upper),
                                        fill = 'grey80', alpha = 1)
                 + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower2, ymax = ~upper2),
                                        fill = 'grey59', alpha = 1)
          )
        }
        p <- ( p + ggplot2::geom_line(data = df, ggplot2::aes_(x = ~r, y = ~curves, group = ~type,
                                    linetype = ~type, size = ~type))
                    + ggplot2::facet_grid('~ main', scales = 'free')
                    + ggplot2::scale_y_continuous(name = ylab, limits = ylim)
                    + ggplot2::scale_linetype_manual(values = linetype.values, name = '')
                    + ggplot2::scale_size_manual(values = size.values, name = '')
                    + ThemePlain(base_size=base_size)
                    )
        if(!is.null(outliers)) {
          outliers.df <- data.frame(r = rep(x[['r']], times=counter),
                                    curves = outliers,
                                    id = rep(1:counter, each=length(x[['r']])))
          p <- p + ggplot2::geom_line(data = outliers.df, ggplot2::aes_(x = ~r, y = ~curves, group = ~id))
        }
        if(rdata$retick_xaxis) {
            p <- p + ggplot2::scale_x_continuous(name = xlab,
                                                 breaks = rdata$loc_break_values,
                                                 labels = paste(round(rdata$r_break_values, digits=2)),
                                                 limits = range(rdata$new_r_values))
            p <- p + ggplot2::geom_vline(xintercept = rdata$new_r_values[rdata$r_values_newstart_id],
                                         linetype = "dotted")
        }
        else p <- p + ggplot2::scale_x_continuous(name = xlab)
    }
    else {
        if(!is.null(ylim)) cat("Note: \"ylim\" ignored as separate yaxes are used.\n")
        n_of_plots <- as.integer(1 + length(rdata$r_values_newstart_id))
        ncols_of_plots <- min(n_of_plots, max_ncols_of_plots)
        nrows_of_plots <- ceiling(n_of_plots / ncols_of_plots)
        if(is.null(labels)) labels <- paste(1:n_of_plots)
        if(length(labels)!=n_of_plots) {
            if(length(labels)==1) {
                labels <- paste(labels, " - ", 1:n_of_plots, sep="")
                warning(paste("Consider giving labels as a vector of length ", n_of_plots,
                              " containing the label for each test function/vector used.\n", sep=""))
            }
            else {
                warning("The length of the vector labels is unreasonable. Setting labels to empty.\n")
                labels <- rep("", times=n_of_plots)
            }
        }

        tmp_indeces <- c(1, rdata$r_values_newstart_id, length(rdata$new_r_values)+1)
        func_labels <- NULL
        for(i in 1:(length(tmp_indeces)-1)) {
            func_labels <- c(func_labels, rep(labels[i], times=tmp_indeces[i+1]-tmp_indeces[i]))
        }

        if(is.null(x[['obs']])) {
          df <- data.frame(r = x[['r']],
                           curves =x[['central']],
                           type = factor("Central function", levels="Central function"),
                           lower = x[['lo']],
                           upper = x[['hi']],
                           main = main,
                           test_function = factor(func_labels, levels=labels))
          if(!is.null(x2)) {
            df$lower2 <- x2[['lo']]
            df$upper2 <- x2[['hi']]
          }
        }
        else {
          df <- data.frame(r = rep(x[['r']], times=2),
                           curves = c(x[['central']], x[['obs']]),
                           type = factor(rep(c("Central function", "Data function"),
                                             each=length(x[['r']])),
                                         levels=c("Data function", "Central function")),
                           lower = rep(x[['lo']], times=2),
                           upper = rep(x[['hi']], times=2),
                           main = main,
                           test_function = factor(func_labels, levels=labels))
          if(!is.null(x2)) {
            df$lower2 <- rep(x2[['lo']], times=2)
            df$upper2 <- rep(x2[['hi']], times=2)
          }
        }
        if(is.null(x2)) {
          p <- ( ggplot2::ggplot()
                 + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower, ymax = ~upper),
                                        fill = 'grey59', alpha = 1)
          )
        }
        else {
          p <- ( ggplot2::ggplot()
                 + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower, ymax = ~upper),
                                        fill = 'grey80', alpha = 1)
                 + ggplot2::geom_ribbon(data = df, ggplot2::aes_(x = ~r, ymin = ~lower2, ymax = ~upper2),
                                        fill = 'grey59', alpha = 1)
          )
        }
        p <- (p + ggplot2::geom_line(data = df, ggplot2::aes_(x = ~r, y = ~curves, group = ~type,
                                    linetype = ~type, size = ~type))
                + ggplot2::facet_wrap(~ test_function, scales="free",
                                      nrow=nrows_of_plots, ncol=ncols_of_plots)
                + ggplot2::scale_y_continuous(name = ylab)
                + ggplot2::scale_linetype_manual(values = linetype.values, name = '')
                + ggplot2::scale_size_manual(values = size.values, name = '')
                + ThemePlain(base_size=base_size)
                + ggplot2::labs(title=main)
        )
        if(!is.null(outliers)) {
          outliers.df <- data.frame(r = rep(x[['r']], times=counter),
                                    curves = outliers,
                                    id = rep(1:counter, each=length(x[['r']])),
                                    test_function = factor(func_labels, levels=labels))
          p <- p + ggplot2::geom_line(data = outliers.df, ggplot2::aes_(x = ~r, y = ~curves, group = ~id))
        }
    }
    # Return
    p
}
