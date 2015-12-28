#' Check r values of a \code{\link{curve_set}} object for plotting
#'
#' Check r values of a \code{\link{curve_set}} object to find out if there is
#' one or several test functions. Find out breaking r values for plotting.
#' @param x A curve_set or envelope_test object.
curve_set_check_r <- function(x) {
    if(!with(x, exists('r'))) stop("The argument \'x\' should contain the element 'r'.\n")
    # Handle combined tests; correct labels on x-axis if x[['r']] contains repeated values
    r_values <- x[['r']]
    nr <- length(r_values)
    if( length(unique(x[['r']])) < nr ) {
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


#' An internal spptest function for setting the default main for a global envelope plot.
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

#' An internal spptest function for setting the default ylim for a global envelope plot.
#' @param x An 'envelope_test' object.
env_ylim_default <- function(x, use_ggplot2) {
    if(!use_ggplot2 || x$alternative != "two.sided")
        ylim <- with(x, c(min(data_curve,lower,upper,central_curve),
                        max(data_curve,lower,upper,central_curve)))
    else ylim <- NULL
    ylim
}

#' An internal spptest function for making a dotplot style "global envelope plot".
#'
#' An internal spptest function for making a dotplot style "global envelope plot".
#'
#' @param x An 'envelope_test' object.
#' @param main See \code{\link{plot.default}}.
#' @param ylim See \code{\link{plot.default}}.
#' @param xlab See \code{\link{plot.default}}.
#' @param ylab See \code{\link{plot.default}}.
#' @param color_outside Logical. Whether to color the places where the data function goes
#' outside the envelope. Currently red color is used.
#' @param ... Additional parameters to be passed to the function \code{\link{plot}}.
env_dotplot <- function(x, main, ylim, xlab, ylab, color_outside, ...) {
    nr <- length(x[['r']])
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
                axis(1, 1:nr, label=paste(round(r, digits=2)))
                points(1:nr, data_curve, pch='x')
                if(color_outside) {
                    outside <- data_curve < lower | data_curve > upper
                    points((1:nr)[outside], data_curve[outside], pch='x', col="red")
                }
            }
    )
}


#' An internal spptest function for making a dotplot style "global envelope plot".
#'
#' An internal spptest function for making a dotplot style "global envelope plot".
#'
#' @param x An 'envelope_test' object.
#' @param main See \code{\link{plot.default}}.
#' @param ylim See \code{\link{plot.default}}.
#' @param xlab See \code{\link{plot.default}}.
#' @param ylab See \code{\link{plot.default}}.
#' @param color_outside Logical. Whether to color the places where the data function goes
#' outside the envelope. Currently red color is used.
#' @param separate_yaxis Logical (default FALSE). By default also the combined envelope plots have
#' a common y-axis. If TRUE, then separate y-axes are used for different parts of a combined test.
#' @param max_ncols_of_plots If separate_yaxis is TRUE, then max_ncols_of_plots gives the maximum
#' number of columns for figures. Default 2.
#' @param ... Additional parameters to be passed to the function \code{\link{plot}}.
env_basic_plot <- function(x, main, ylim, xlab, ylab, color_outside,
                           separate_yaxis=FALSE, max_ncols_of_plots=2, ...) {
    # Handle combined tests; correct labels on x-axis if x[['r']] contains repeated values
    nr <- length(x[['r']])
    rdata <- curve_set_check_r(x)
    if(rdata$retick_xaxis) x[['r']] <- 1:nr
    # Plot
    if(!separate_yaxis) {
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
                        axis(1, rdata$loc_break_values, label=paste(round(rdata$r_break_values, digits=2)))
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
        par(mfrow=c(nrows_of_plots, ncols_of_plots))
        with(x, {
                    cat("Note: \"main\" and \"ylim\" ignored as separate plots are produced.\n")
                    tmp_indeces <- c(1, rdata$r_values_newstart_id, length(rdata$new_r_values)+1)
                    for(i in 1:n_of_plots) {
                        plot(r[tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             data_curve[tmp_indeces[i]:(tmp_indeces[i+1]-1)],
                             main="", xlab=xlab[i], ylab=ylab[i],
                             type="l", lty=1, lwd=2, ...)
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
