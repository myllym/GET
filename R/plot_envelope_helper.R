#' Check r values of a curve_set object for plotting
#'
#' Check r values of a curve_set object to find out if there is
#' one or several test functions. Find out breaking r values for plotting.
#' @param x A curve_set or envelope_test object.
#' @seealso \code{\link{create_curve_set}}
curve_set_check_r <- function(x) {
    if(!with(x, exists('r'))) stop("The argument \'x\' should contain the element 'r'.\n")
    # Handle combined tests; correct labels on x-axis if x[['r']] contains repeated values
    r_values <- x[['r']]
    nr <- length(r_values)
    if( !all(r_values[-1] - r_values[-nr] > 0) ) {
        retick_xaxis <- TRUE
        new_r_values <- 1:nr # to be used in plotting
        r_values_newstart_id <- which(!(r_values[1:(nr-1)] < r_values[2:nr])) + 1
        # number of ticks per a sub test
        nticks <- 5
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
    else {
        retick_xaxis <- FALSE
        new_r_values <- NULL
        r_break_values <- NULL
        loc_break_values <- NULL
        r_values_newstart_id <- NULL
    }
    
    list(r_values = r_values,
         retick_xaxis = retick_xaxis,
         new_r_values = new_r_values,
         r_break_values = r_break_values, loc_break_values = loc_break_values,
         r_values_newstart_id = r_values_newstart_id)
}


#' An internal GET function for setting the default main for a global envelope plot.
#' @param x An 'envelope_test' object.
env_main_default <- function(x) {
    if(with(x, exists('p_interval')))
        if(x$alternative == "two.sided")
            main <- with(x, paste(method, ": p-interval = (",
                            round(p_interval[1],3),", ", round(p_interval[2],3), ")", sep=""))
        else
            main <- with(x, paste(method, ": p-interval = (",
                            round(p_interval[1],3),", ", round(p_interval[2],3), ") \n",
                            "Alternative = \"", alternative, "\"\n", sep=""))
    else {
        if(x$alternative == "two.sided")
            main <- with(x, paste(method, ": p = ", round(p,3), sep=""))
        else
            main <- with(x, paste(method, ": p = ", round(p,3), "\n",
                            "Alternative = \"", alternative, "\"\n", sep=""))
    }
    main
}

#' An internal GET function for setting the default ylim for a global envelope plot.
#' @param x An 'envelope_test' object.
#' @param use_ggplot2 TRUE/FALSE, If TRUE, then default ylim are for \code{\link{env_ggplot}}.
#' Otherwise the ylim are for \code{\link{env_basic_plot}}.
env_ylim_default <- function(x, use_ggplot2) {
    if(!use_ggplot2 || x$alternative != "two.sided")
        ylim <- with(x, c(min(data_curve,lower,upper,central_curve),
                        max(data_curve,lower,upper,central_curve)))
    else ylim <- NULL
    ylim
}

#' An internal GET function for making a dotplot style "global envelope plot".
#'
#' An internal GET function for making a dotplot style "global envelope plot".
#'
#' @param x An 'envelope_test' object.
#' @param main See \code{\link{plot.default}}.
#' @param ylim See \code{\link{plot.default}}.
#' @param xlab See \code{\link{plot.default}}.
#' @param ylab See \code{\link{plot.default}}.
#' @param color_outside Logical. Whether to color the places where the data function goes
#' outside the envelope. Currently red color is used.
#' @param labels Labels for the tests at x-axis.
#' @param ... Additional parameters to be passed to the function \code{\link{plot}}.
env_dotplot <- function(x, main, ylim, xlab, ylab, color_outside=TRUE, labels, ...) {
    nr <- length(x[['r']])
    if(missing(labels)) labels <- paste(round(x[['r']], digits=2))
    if(nr > 10) warning("Dotplot style meant for low dimensional test vectors.\n")
    with(x, {
                plot(1:nr, central_curve, main=main, ylim=ylim, xlab=xlab, ylab=ylab, cex=0.5, pch=16, xaxt="n", ...)
                if(alternative!="greater")
                    arrows(1:nr, lower, 1:nr, central_curve, code = 1, angle = 75, length = .1)
                else
                    arrows(1:nr, lower, 1:nr, central_curve, code = 1, angle = 75, length = .1, col=grey(0.8))
                if(alternative!="less")
                    arrows(1:nr, upper, 1:nr, central_curve, code = 1, angle = 75, length = .1)
                else
                    arrows(1:nr, upper, 1:nr, central_curve, code = 1, angle = 75, length = .1, col=grey(0.8))
                axis(1, 1:nr, labels=labels)
                points(1:nr, data_curve, pch='x')
                if(color_outside) {
                    outside <- data_curve < lower | data_curve > upper
                    points((1:nr)[outside], data_curve[outside], pch='x', col="red")
                }
            }
    )
}


#' An internal GET function for making a basic "global envelope plot".
#'
#' An internal GET function for making a basic "global envelope plot".
#'
#' @param x An 'envelope_test' object.
#' @param main See \code{\link{plot.default}}.
#' @param ylim See \code{\link{plot.default}}.
#' @param xlab See \code{\link{plot.default}}.
#' @param ylab See \code{\link{plot.default}}.
#' @param color_outside Logical. Whether to color the places where the data function goes
#' outside the envelope. Currently red color is used.
#' @param separate_yaxes Logical (default FALSE). By default also the combined envelope plots have
#' a common y-axis. If TRUE, then separate y-axes are used for different parts of a combined test.
#' @param max_ncols_of_plots If separate_yaxes is TRUE, then max_ncols_of_plots gives the maximum
#' number of columns for figures. Default 2.
#' @param ... Additional parameters to be passed to the function \code{\link{plot}}.
#' @importFrom graphics par
env_basic_plot <- function(x, main, ylim, xlab, ylab, color_outside=TRUE,
                           separate_yaxes=FALSE, max_ncols_of_plots=2, ...) {
    # Handle combined tests; correct labels on x-axis if x[['r']] contains repeated values
    nr <- length(x[['r']])
    rdata <- curve_set_check_r(x)
    # Plot
    if(!separate_yaxes) {
        if(rdata$retick_xaxis) x[['r']] <- 1:nr
        with(x, {
                    if(!rdata$retick_xaxis)
                        plot(r, data_curve, main=main, ylim=ylim, xlab=xlab, ylab=ylab,
                                type="l", lty=1, lwd=2, ...)
                    else
                        plot(r, data_curve, ylim=ylim, main=main, xlab=xlab, ylab=ylab,
                                type="l", lty=1, lwd=2, xaxt="n", ...)
                    if(alternative!="greater") lines(r, lower, lty=2) else lines(r, lower, lty=2, col=grey(0.8))
                    if(alternative!="less") lines(r, upper, lty=2) else lines(r, upper, lty=2, col=grey(0.8))
                    lines(r, central_curve, lty=3)
                    if(color_outside) {
                        outside <- data_curve < lower | data_curve > upper
                        points(r[outside], data_curve[outside], col="red")
                    }
                    if(rdata$retick_xaxis) {
                        axis(1, rdata$loc_break_values, labels=paste(round(rdata$r_break_values, digits=2)))
                        abline(v = rdata$new_r_values[rdata$r_values_newstart_id], lty=3)
                    }
                }
        )
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
        with(x, {
                    cat("Note: \"main\" and \"ylim\" ignored as separate plots are produced.\n")
                    tmp_indeces <- c(1, rdata$r_values_newstart_id, length(rdata$new_r_values)+1)
                    for(i in 1:n_of_plots) {
                        ylim <- c(min(data_curve[tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                                      lower[tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                                      upper[tmp_indeces[i]:(tmp_indeces[i+1]-1)]),
                                  max(data_curve[tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                                      lower[tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                                      upper[tmp_indeces[i]:(tmp_indeces[i+1]-1)]))
                        plot(r[tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             data_curve[tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             main="", xlab=xlab[i], ylab=ylab[i],
                             type="l", lty=1, lwd=2, ylim=ylim, ...)
                        if(alternative!="greater")
                            lines(r[tmp_indeces[i]:(tmp_indeces[i+1]-1)], lower[tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2)
                        else lines(r[tmp_indeces[i]:(tmp_indeces[i+1]-1)], lower[tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2, col=grey(0.8))
                        if(alternative!="less")
                            lines(r[tmp_indeces[i]:(tmp_indeces[i+1]-1)], upper[tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2)
                        else lines(r[tmp_indeces[i]:(tmp_indeces[i+1]-1)], upper[tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=2, col=grey(0.8))
                        lines(r[tmp_indeces[i]:(tmp_indeces[i+1]-1)], central_curve[tmp_indeces[i]:(tmp_indeces[i+1]-1)], lty=3)
                        if(color_outside) {
                            outside <- data_curve[tmp_indeces[i]:(tmp_indeces[i+1]-1)] < lower[tmp_indeces[i]:(tmp_indeces[i+1]-1)] | data_curve[tmp_indeces[i]:(tmp_indeces[i+1]-1)] > upper[tmp_indeces[i]:(tmp_indeces[i+1]-1)]
                            points(r[tmp_indeces[i]:(tmp_indeces[i+1]-1)][outside], data_curve[tmp_indeces[i]:(tmp_indeces[i+1]-1)][outside], col="red")
                        }
                    }
                }
        )
    }
}


#' An internal GET function for making a ggplot2 style "global envelope plot".
#'
#' An internal GET function for making a ggplot2 style "global envelope plot".
#'
#' @param x An 'envelope_test' object.
#' @param base_size Base font size, to be passed to theme style when \code{use_ggplot2 = TRUE}.
#' @param main See \code{\link{plot.default}}.
#' @param ylim See \code{\link{plot.default}}.
#' @param xlab See \code{\link{plot.default}}.
#' @param ylab See \code{\link{plot.default}}.
#' @param separate_yaxes Logical (default FALSE). By default also the combined envelope plots have
#' a common y-axis. If TRUE, then separate y-axes are used for different parts of a combined test.
#' @param labels Labels for the separate plots. Ignored if separate_yaxes is FALSE.
#' @param max_ncols_of_plots If separate_yaxes is TRUE, then max_ncols_of_plots gives the maximum
#' number of columns for figures. Default 2.
#' @import ggplot2
env_ggplot <- function(x, base_size, main, ylim, xlab, ylab, separate_yaxes=FALSE, max_ncols_of_plots=2,
                       labels=NULL) {
    # Handle combined tests; correct labels on x-axis if x[['r']] contains repeated values
    nr <- length(x[['r']])
    rdata <- curve_set_check_r(x)

    linetype.values <- c('solid', 'dashed')
    size.values <- c(0.2, 0.2)

    if(!separate_yaxes | is.null(rdata$r_values_newstart_id)) {
        if(rdata$retick_xaxis) x[['r']] <- 1:nr
        with(x, {
                    df <- data.frame(r = rep(r, times=2),
                            curves = c(data_curve, central_curve),
                            type = factor(rep(c("Data function", "Central function"), each=length(r)), levels=c("Data function", "Central function")),
                            lower = rep(lower, times=2),
                            upper = rep(upper, times=2),
                            main = factor(rep(main, times=length(r)))
                    )
                    p <- (ggplot2::ggplot()
                                + ggplot2::geom_ribbon(data = df, ggplot2::aes(x = r, ymin = lower, ymax = upper),
                                        fill = 'grey59', alpha = 1)
                                + ggplot2::geom_line(data = df, ggplot2::aes(x = r, y = curves, group = type,
                                                linetype = type, size = type))
                                + ggplot2::facet_grid('~ main', scales = 'free')
                                + ggplot2::scale_y_continuous(name = ylab, limits = ylim)
                                + ggplot2::scale_linetype_manual(values = linetype.values, name = '')
                                + ggplot2::scale_size_manual(values = size.values, name = '')
                                + ThemePlain(base_size=base_size)
                                )
                    if(rdata$retick_xaxis) {
                        p <- p + ggplot2::scale_x_continuous(name = xlab,
                                breaks = rdata$loc_break_values,
                                labels = paste(round(rdata$r_break_values, digits=2),
                                        limits = range(rdata$new_r_values)))
                        p <- p + ggplot2::geom_vline(xintercept = rdata$new_r_values[rdata$r_values_newstart_id], linetype = "dotted")
                    }
                    else p <- p + ggplot2::scale_x_continuous(name = xlab)
                    print(p)
                    return(invisible(p))
                }
        )
    }
    else {
        cat("Note: \"ylim\" ignored as separate plots are produced.\n")
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

        with(x, {
                    df <- data.frame(r = rep(r, times=2),
                            curves = c(data_curve, central_curve),
                            type = factor(rep(c("Data function", "Central function"), each=length(r)), levels=c("Data function", "Central function")),
                            lower = rep(lower, times=2),
                            upper = rep(upper, times=2),
                            main = factor(rep(main, times=length(r))),
                            test_function = factor(rep(func_labels, times=2), levels=labels)
                    )
                    p <- (ggplot2::ggplot()
                                + ggplot2::geom_ribbon(data = df, ggplot2::aes(x = r, ymin = lower, ymax = upper),
                                        fill = 'grey59', alpha = 1)
                                + ggplot2::geom_line(data = df, ggplot2::aes(x = r, y = curves, group = type,
                                                linetype = type, size = type))
                                + ggplot2::facet_wrap(~ test_function, scales="free",
                                                      nrow=nrows_of_plots, ncol=ncols_of_plots)
                                + ggplot2::scale_y_continuous(name = ylab)
                                + ggplot2::scale_linetype_manual(values = linetype.values, name = '')
                                + ggplot2::scale_size_manual(values = size.values, name = '')
                                + ThemePlain(base_size=base_size)
                                + labs(title=main)
                                )
                    print(p)
                    return(invisible(p))
                }
        )
    }
}

#' A helper function for plotting two global envelopes into a same plot.
#'
#' @param env1 An 'envelope_test' or 'adjusted_envelope_test' object. In essence this object
#' must contain arguments '$r', '$data_curve', '$lower', '$upper' and '$central_curve$.
#' @param env2 An 'envelope_test' or 'adjusted_envelope_test' object. In essence this object
#' must contain arguments '$r', '$data_curve', '$lower', '$upper' and '$central_curve$.
#' @param base_size Base font size, to be passed to theme style.
#' @param main See \code{\link{plot.default}}.
#' @param ylim See \code{\link{plot.default}}.
#' @param xlab See \code{\link{plot.default}}.
#' @param ylab See \code{\link{plot.default}}.
#' @export
#' @import ggplot2
two_envelopes_ggplot <- function(env1, env2, base_size=15, main, ylim, xlab, ylab) {
    if(!(class(env1) %in% c("envelope_test", "adjusted_envelope_test")) |
       !(class(env2) %in% c("envelope_test", "adjusted_envelope_test"))) stop("env1 and/or env2 is not desired object type.\n")
    if(!all(env1$r == env2$r)) stop("The two envelopes are for different r-values.\n")
    linetype.values <- c('solid', 'dashed')
    size.values <- c(0.2, 0.2)
    if(missing(xlab)) xlab <- expression(italic(r))
    if(missing(ylab)) ylab <- expression(italic(T(r)))
    if(missing(main)) main <- "Rank envelope test"
    if(missing(ylim)) ylim <- c(min(env1$data_curve,env1$lower,env1$upper,env1$central_curve,
                                    env2$data_curve,env2$lower,env2$upper,env2$central_curve),
                                max(env1$data_curve,env1$lower,env1$upper,env1$central_curve,
                                    env2$data_curve,env2$lower,env2$upper,env2$central_curve))
    df <- data.frame(r = rep(env1$r, times=2),
            curves = c(env1$data_curve, env1$central_curve),
            type = factor(rep(c("Data function", "Central function"), each=length(env1$r)), levels=c("Data function", "Central function")),
            lower = rep(env1$lower, times=2),
            upper = rep(env1$upper, times=2),
            lower2 = rep(env2$lower, times=2),
            upper2 = rep(env2$upper, times=2),
            main = factor(rep(main, times=length(env1$r)))
    )
    p <- (ggplot2::ggplot()
                + ggplot2::geom_ribbon(data = df, ggplot2::aes_string(x = 'r', ymin = 'lower2', ymax = 'upper2'),
                        fill = 'grey80', alpha = 1)
                + ggplot2::geom_ribbon(data = df, ggplot2::aes_string(x = 'r', ymin = 'lower', ymax = 'upper'),
                        fill = 'grey59', alpha = 1)
                + ggplot2::geom_line(data = df, ggplot2::aes_string(x = 'r', y = 'curves', group = 'type',
                                linetype = 'type', size = 'type'))
                + ggplot2::facet_grid('~ main', scales = 'free')
                + ggplot2::scale_x_continuous(name = xlab)
                + ggplot2::scale_y_continuous(name = ylab, limits = ylim)
                + ggplot2::scale_linetype_manual(values = linetype.values, name = '')
                + ggplot2::scale_size_manual(values = size.values, name = '')
                + ThemePlain(base_size=base_size)
                )
    print(p)
    invisible(p)
}
