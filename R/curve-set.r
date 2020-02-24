#' Turn an \code{\link[spatstat]{envelope}} object into a curve_set object.
#'
#' @param env An \code{\link[spatstat]{envelope}} object. The envelope()
#'   functions must have been called with savefuns = TRUE.
#' @return A corresponding curve_set object.
#' @param ... Do not use. (For internal use only.)
#' @export
envelope_to_curve_set <- function(env, ...) {
  if(!inherits(env, 'envelope')) {
    stop('env is not an envelope object.')
  }

  r <- env[['r']]
  n_r <- length(r)
  if(n_r < 1L) {
    stop('env[["r"]] must exist.')
  }

  obs <- env[['obs']]
  if(length(obs) != n_r) {
    stop('env[["obs"]] and env[["r"]] must have the same length.')
  }

  simfuns_fv <- attr(env, 'simfuns', exact = TRUE)
  if(length(simfuns_fv) < 1L) {
    stop('env did not include the simulated curve set. envelope must ',
         'be run with savefuns = TRUE.')
  }
  simulation_df <- as.data.frame(simfuns_fv)
  simulation_df_names <- names(simulation_df)

  if(!('r' %in% simulation_df_names)) {
    stop('The env attribute simfuns did not include r.')
  }
  if(!identical(r, simulation_df[['r']])) {
    stop('env[["r"]] must be equal to ',
         'attributes(env)[["simfuns"]][["r"]].')
  }

  # Dimensions: r_idx, sim_idx.
  sim_m <- as.matrix(simulation_df[, !(simulation_df_names %in% 'r')])
  if(nrow(sim_m) != n_r) {
    stop('simulated curves must have the same length as r.')
  }

  theo <- env[['theo']]
  n_theo <- length(theo)
  # If theo exists it must be correct.
  if(n_theo > 0L && n_theo != n_r) {
    stop('env[["theo"]] and env[["r"]] must have the same length.')
  }

  res <- list(r = r, obs = obs, sim_m = sim_m)
  if(n_theo > 0L) {
    res[['theo']] <- theo
  }
  res[['is_residual']] <- FALSE

  res <- create_curve_set(res, ...)
  res
}

# Turn an \code{\link[fda.usc]{fdata}} object into a curve_set object.
#
# @param fdata An \code{\link[fda.usc]{fdata}} object.
# @return A corresponding curve_set object.
# @param ... Ignored.
fdata_to_curve_set <- function(fdata, ...) {
  if(!inherits(fdata, 'fdata')) {
    stop('fdata is not a fdata object.')
  }
  if(!is.null(fdata[['fdata2d']]) && fdata[['fdata2d']])
    stop("No support for converting this type of functional data, only functional data observed in a single grid supported.\n")

  r <- fdata[['argvals']]
  n_r <- length(r)
  if(n_r < 1L) {
    stop('fdata[["argvals"]] must exist.')
  }

  obs <- t(fdata[['data']])
  if(nrow(obs) != n_r) {
    stop('Number of functions in fdata[["data"]] and the length of fdata[["argvals"]] must be equal.')
  }

  res <- list(r = r, obs = obs)

  res <- create_curve_set(res, ...)
  res
}

# Check the content validity of a potential curve_set object.
#
# @param curve_set An object to be checked.
# @param allow_Inf_values Logical. Can be used to allow infinite or nonnumeric
# values in an \code{\link[spatstat]{envelope}} object at the first place, if those are cropped
# away (in \code{\link{crop_curves}}).
check_curve_set_content <- function(curve_set, allow_Inf_values = FALSE) {
  possible_names <- c('r', 'obs', 'sim_m', 'theo', 'is_residual')

  n <- length(curve_set)
  if(n < 1L) {
    stop('curve_set must have some elements.')
  }
  if(!is.list(curve_set)) {
    stop('curve_set must be a list.')
  }

  name_vec <- names(curve_set)
  if(length(name_vec) != n) {
    stop('Every element of curve_set must be named.')
  }
  if(!all(name_vec %in% possible_names)) {
    stop('curve_set should contain only elements with some of these ',
         ' names: ', paste(possible_names, collapse = ', '))
  }

  r <- curve_set[['r']]
  n_r <- length(r)
  if(n_r < 1L) {
    stop('curve_set[["r"]] must have at least one element.')
  }
  if(!is.vector(r)) {
    stop('curve_set[["r"]] must be a vector.')
  }
  if(!all(is.numeric(r)) || !all(is.finite(r))) {
    stop('curve_set[["r"]] must have only finite numeric values.')
  }

  obs <- curve_set[['obs']]
  if(!is.vector(obs) & !is.matrix(obs)) {
    stop('curve_set[["obs"]] must be a vector or a matrix.')
  }
  if((is.vector(obs) && length(obs) != n_r) | (is.matrix(obs) && nrow(obs) != n_r)) {
    stop('curve_set[["obs"]] must have as many values or rows as ',
         'curve_set[["r"]].')
  }
  if(!(allow_Inf_values | ( all(is.numeric(obs)) && all(is.finite(obs)) ))) {
    stop('curve_set[["obs"]] must have only finite numeric values.')
  }

  sim_m <- curve_set[['sim_m']]
  if(length(sim_m) > 0L) {
    dim_sim_m <- dim(sim_m)
    if(!is.matrix(sim_m)) {
      stop('curve_set[["sim_m"]] must be a matrix.')
    }
    if(dim_sim_m[1] != n_r) {
      stop('curve_set[["sim_m"]] must have as many rows as there are ',
           'elements in curve_set[["r"]].')
    }
    if(dim_sim_m[2] < 1L) {
      stop('curve_set[["sim_m"]] must have at least one column.')
    }
    if(!(allow_Inf_values | ( all(is.numeric(sim_m)) && all(is.finite(sim_m)) ))) {
      stop('curve_set[["sim_m"]] must have only finite numeric values.')
    }
  }

  theo <- curve_set[['theo']]
  n_theo <- length(theo)
  if(n_theo > 0L) {
    if(n_theo != n_r) {
      stop('curve_set[["theo"]] must have as many values as ',
           'curve_set[["r"]].')
    }
    if(!is.vector(theo)) {
      stop('curve_set[["theo"]] must be a vector.')
    }
    if(!(allow_Inf_values | ( all(is.numeric(theo)) && all(is.finite(theo)) ))) {
      stop('curve_set[["theo"]] must have only finite numeric ',
           'values.')
    }
  }

  is_residual <- curve_set[['is_residual']]
  n_is_residual <- length(is_residual)
  if(n_is_residual > 0L && (n_is_residual != 1L ||
                             !is.logical(is_residual) ||
                             !is.finite(is_residual))) {
    stop('curve_set[["is_residual"]] must be either TRUE or FALSE.')
  }

  if(n_is_residual > 0L && is_residual && n_theo > 0L) {
    stop('A residual curve set must not contain a theoretical curve.')
  }

  curve_set
}

# Convert an envelope object to a curve_set object.
#
# If given an envelope object, convert it into a curve_set object. If given
# a curve_set object, check its correctness and give it back.
#
# @inheritParams crop_curves
# @param ... Allows to pass arguments to \code{\link{check_curve_set_content}}
# and \code{\link{envelope_to_curve_set}} (to be passed further through
# \code{\link{create_curve_set}} to \code{\link{check_curve_set_content}}).
# @return If an \code{\link[spatstat]{envelope}} object was given, return a
#   corresponding curve_set object. If a curve_set object was given, return
#   it unharmed.
#' @importFrom methods is
convert_envelope <- function(curve_set, ...) {
  if(inherits(curve_set, 'envelope')) {
    curve_set <- envelope_to_curve_set(curve_set, ...)
  }
  else if(!methods::is(curve_set, 'curve_set')) {
    stop('curve_set must either have class "envelope" (from spatstat) ',
         'or class "curve_set".')
  }
  check_curve_set_content(curve_set, ...)
  curve_set
}

# Convert a fdata object to a curve_set object.
#
# If given a fdata object, convert it into a curve_set object. If given
# a curve_set object, check its correctness and give it back.
#
# @inheritParams crop_curves
# @param ... Allows to pass arguments to \code{\link{check_curve_set_content}}
# and \code{\link{envelope_to_curve_set}} (to be passed further through
# \code{\link{create_curve_set}} to \code{\link{check_curve_set_content}}).
# @return If an \code{\link[fda.usc]{fdata}} object was given, return a
#   corresponding curve_set object. If a curve_set object was given, return
#   it unharmed.
#' @importFrom methods is
convert_fdata <- function(curve_set, ...) {
  if(inherits(curve_set, 'fdata')) {
    curve_set <- fdata_to_curve_set(curve_set, ...)
  } else if(!methods::is(curve_set, 'curve_set')) {
    stop('curve_set must either have class "fdata" (from fda.usc) ',
         'or class "curve_set".')
  }
  check_curve_set_content(curve_set, ...)
  curve_set
}

# Check that the curve_set object portrays residual curves.
# @param curve_set A 'curve_set' object
check_residualness <- function(curve_set) {
  is_residual <- curve_set[['is_residual']]
  if(length(is_residual) < 1L || !is_residual) {
    stop('curve_set must consist of residual curves. Run function ',
         'residual first.')
  }
  invisible(curve_set)
}

#' Create a curve set
#'
#' Create a curve set out of a list in the right form.
#'
#'
#' \code{obs} must be either
#' \itemize{
#' \item a vector containing the data function, or
#' \item a matrix containing the s data functions, in which case it is assumed that
#' each column corresponds to a data function.
#' }
#'
#' If given, \code{r} describes the 1- or 2-dimensional argument values where the curves have been observed (or
#' simulated). It must be either
#' \itemize{
#' \item a vector,
#' \item a data.frame with columns "x", "y", "width" and "height",
#' where the width and height give the width and height of the pixels placed at x and y, or
#' \item a data.frame with columns "xmin", "xmax", "ymin" and "ymax" giving the corner coordinates of the pixels
#' where the data have been observed.
#' }
#'
#' If \code{obs} is a vector, \code{sim_m} must be a matrix containing the simulated functions.
#' Each column is assumed to correspond to a function, and the number of rows must match the
#' length of \code{obs}. If \code{obs} is a matrix, \code{sim_m} is ignored.
#'
#' If \code{obs} is a vector, \code{theo} can be given and it should then correspond
#' to a theoretical function (e.g., under the null hypothesis). If present, its length must match the length of
#' \code{obs}.
#' @param curve_set A list containing the element obs, and optionally
#'   the elements r, sim_m and theo. See details.
#' @param ... Do not use. (For internal use only.)
#' @return The given list adorned with the proper class name.
#' @export
create_curve_set <- function(curve_set, ...) {
  check_curve_set_content(curve_set, ...)
  class(curve_set) <- 'curve_set'
  curve_set
}

#' Check class.
#'
#' @param x An object to be checked.
#' @export
is.curve_set <- function(x) inherits(x, 'curve_set')


#' Print method for the class 'curve_set'
#'
#' @param x an 'curve_set' object
#' @param ... Passed to \code{\link[utils]{str}}.
#'
#' @export
#' @importFrom utils str
print.curve_set <- function(x, ...) {
  str(x, ...)
}

#' Plot method for the class 'curve_set'
#'
#' @param x An \code{curve_set} object
#' @param plot_style Either "ggplot2" or "basic".
#' @param ylim The y limits of the plot with the default being the minimum and maximum over all curves.
#' @param xlab The label for the x-axis. Default "r".
#' @param ylab The label for the y-axis. Default "obs".
#' @param col_obs Color for 'obs' in the argument \code{x}.
#' @param col_sim Color for 'sim_m' in the argument \code{x}.
#' @param ... Additional parameters to be passed to plot and lines.
#' @inheritParams plot.global_envelope
#'
#' @export
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom grDevices grey
#' @importFrom graphics axis
#' @importFrom graphics abline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 labs
plot.curve_set <- function(x, plot_style = c("ggplot2", "basic"),
                           ylim, xlab = "r", ylab = "obs", main = NULL,
                           col_obs = 1, col_sim = grDevices::grey(0.7),
                           base_size = 11, ...) {
  plot_style <- match.arg(plot_style)
  if(with(x, is.matrix(obs))) {
    funcs <- x[['obs']]
    col_sim <- col_obs
  }
  else funcs <- cbind(x[['obs']], x[['sim_m']])
  rdata <- combined_global_envelope_rhelper(x)
  if(rdata$retick_xaxis) {
    rvalues <- rdata$new_r_values
  }
  else rvalues <- x$r
  nr <- length(rvalues)
  if(missing('ylim')) {
    if(plot_style == "basic") ylim <- with(x, c(min(funcs), max(funcs)))
    else ylim <- NULL
  }

  switch(plot_style,
         ggplot2 = {
             df <- data.frame(r = rvalues, f = c(funcs), id = rep(1:ncol(funcs), each=nrow(funcs)))
             p <- ( ggplot(data=df) + geom_line(aes_(x = ~r, y = ~f, group = ~id))
                   + scale_y_continuous(name = ylab, limits = ylim)
                   + labs(title=main)
                   + ThemePlain(base_size=base_size))
             if(rdata$retick_xaxis) {
               p <- p + scale_x_continuous(name = xlab,
                                           breaks = rdata$loc_break_values,
                                           labels = paste(round(rdata$r_break_values, digits=2)),
                                           limits = range(rdata$new_r_values))
               p <- p + geom_vline(xintercept = rdata$new_r_values[rdata$r_values_newstart_id],
                                  linetype = "dotted")
             }
             else p <- p + scale_x_continuous(name = xlab)
             p
         },
         basic = {
             # Plot
             if(!rdata$retick_xaxis)
                 graphics::plot(rvalues, funcs[,1], type="l", ylim=ylim, col=col_obs, xlab=xlab, ylab=ylab, main=main, ...)
             else
                 graphics::plot(rvalues, funcs[,1], type="l", ylim=ylim, xaxt="n", col=col_obs, xlab=xlab, ylab=ylab, main=main, ...)
             for(i in 2:ncol(funcs)) graphics::lines(rvalues, funcs[,i], col=col_sim)
             graphics::lines(rvalues, funcs[,1], type="l", col=col_obs, ...)
             if(rdata$retick_xaxis) {
                 graphics::axis(1, rdata$loc_break_values, labels=paste(round(rdata$r_break_values, digits=2)))
                 graphics::abline(v = rdata$new_r_values[rdata$r_values_newstart_id], lty=3)
             }
         })
}

# Combine curve sets.
#
# Combine curve sets to a one curve set, e.g. for testing by means of several test functions.
# @param x A list of curve sets or \code{\link[spatstat]{envelope}} objects.
# @param equalr Whether to demand equal lengths of r vectors of the different curve sets
# @return A curve set that is a combination of the curve sets given in 'x'.
combine_curve_sets <- function(x, equalr = TRUE) {
  cset <- NULL
  x <- check_curve_set_dimensions(x, equalr=equalr)
  name_vec <- lapply(x, FUN=names)
  if('r' %in% name_vec[[1]]) {
    cset$r <- c(sapply(x, FUN=function(curve_set) { curve_set['r'] }), recursive=TRUE)
  }
  if('obs' %in% name_vec[[1]]) {
    if(is.matrix(x[[1]]$obs)) {
      cset$obs <- matrix(nrow=sum(sapply(x, FUN=function(curve_set) {dim(curve_set$obs)[1]})), ncol=dim(x[[1]]$obs)[2])
      for(i in 1:dim(x[[1]]$obs)[2]) {
        cset$obs[,i] <- c(sapply(x, FUN=function(curve_set) { curve_set$obs[,i] }), recursive=TRUE)
      }
    }
    else {
      cset$obs <- c(sapply(x, FUN=function(curve_set) { curve_set['obs'] }), recursive=TRUE)
    }
  }
  if('sim_m' %in% name_vec[[1]]) {
    # check_curve_set_dimensions was used above to check that dimensions match
    # Combine
    cset$sim_m <- matrix(nrow=sum(sapply(x, FUN=function(curve_set) {dim(curve_set$sim_m)[1]})), ncol=dim(x[[1]]$sim_m)[2])
    for(i in 1:dim(x[[1]]$sim_m)[2]) {
      cset$sim_m[,i] <- c(sapply(x, FUN=function(curve_set) { curve_set$sim_m[,i] }), recursive=TRUE)
    }
  }
  if('theo' %in% name_vec[[1]])
    cset$theo <- c(sapply(x, FUN=function(curve_set) { curve_set['theo'] }), recursive=TRUE)
  # check_curve_set_dimensions has checked that 'is_residual' is the same for all curve sets.
  if('is_residual' %in% name_vec[[1]])
    cset$is_residual <- x[[1]]$is_residual
  # Return the combined set of curves
  create_curve_set(cset)
}

# Check curve set dimensions
#
# Check that x contains list of curve sets or \code{\link[spatstat]{envelope}} objects.
# If the latter, then convert the objects to curve sets.
# Check that the curve sets have same elements and dimensions of them (numbers of r-values can differ for equalr=FALSE).
# @inheritParams combine_curve_sets
check_curve_set_dimensions <- function(x, equalr=FALSE) {
  x <- lapply(x, FUN=convert_envelope)
  name_vec <- lapply(x, FUN=function(x) { n <- names(x); n[n!="theo"] })
  # Possible_names in curve_sets are 'r', 'obs', 'sim_m', 'theo' and 'is_residual'.
  # Check that all curve sets contain the same elements
  if(!all(sapply(name_vec, FUN=identical, y=name_vec[[1]])))
    stop("The curve sets in \'x\' contain different elements.\n")
  if(equalr) {
    # Check equality of length of r
    if(!all(sapply(x, FUN=function(curve_set) { length(curve_set$r) == length(x[[1]]$r) })))
      stop("The lengths of 'r' in curve sets differ.\n")
  }
  # Check that 'is_residual' is the same for all curve sets.
  if('is_residual' %in% name_vec[[1]]) {
    if(!all(sapply(x, FUN=function(curve_set) { curve_set$is_residual == x[[1]]$is_residual })))
      stop("The element \'is_residual\' should be the same for each curve set.\n")
  }
  # Check that the number of functions in 'obs' in curve sets are equal.
  if('obs' %in% name_vec[[1]]) {
    if(all(sapply(x, FUN=is.matrix))) {
      if(!all(sapply(x, FUN=function(curve_set) { dim(curve_set$obs)[2] == dim(x[[1]]$obs)[2] })))
        stop("The numbers of functions in 'obs' in curve sets differ.\n")
    }
  }
  # Check that the number of functions in 'sim_m' in curve sets are equal.
  if('sim_m' %in% name_vec[[1]]) {
    if(!all(sapply(x, FUN=function(curve_set) { dim(curve_set$sim_m)[2] == dim(x[[1]]$sim_m)[2] })))
      stop("The numbers of functions in 'sim_m' in curve sets differ.\n")
  }
  x
}

# A helper function to return all the functions from a curve set in a matrix.
# If obs is a vector, then the returned matrix will contain the obs vector on its first row.
# If obs is a matrix, then the returned matrix will be a transpose of obs.
data_and_sim_curves <- function(curve_set) {
  if(with(curve_set, is.matrix(obs))) funcs <- t(curve_set[['obs']]) # all rows data (sim_m ignored)
  else {
    funcs <- rbind(curve_set[['obs']], t(curve_set[['sim_m']])) # first row data, rest simulations
    rownames(funcs) <- c('obs', paste("sim", 1:(nrow(funcs)-1), sep=""))
  }
  funcs
}
# Return all the functions from a curve set in a matrix.
# Columns: functions, Rows: values of the functions at an argument value. (same as in curve_set$obs)
curve_set_funcs <- function(curve_set) {
  if(with(curve_set, is.matrix(obs))) funcs <- curve_set[['obs']]
  else {
    funcs <- cbind(curve_set[['obs']], curve_set[['sim_m']])
    colnames(funcs) <- c('obs', paste("sim", 1:(ncol(funcs)-1), sep=""))
  }
  funcs
}

# A helper function to obtain the mean of functions in curve_set.
#
# Calculate the mean of all the functions in the curve_set (obs and sim_m).
curve_set_mean <- function(curve_set) {
  funcs <- data_and_sim_curves(curve_set)
  colMeans(funcs)
}

# A helper function to obtain the median of functions in curve_set.
#
# Calculate the median of all the functions in the curve_set (obs and sim_m).
#' @importFrom stats median
curve_set_median <- function(curve_set) {
  funcs <- data_and_sim_curves(curve_set)
  sapply(1:ncol(funcs), function(i) { median(funcs[,i]) })
}

# A helper function to obtain the sd of functions in curve_set.
#
# No matter whether obs is a matrix or a vector, the quantiles are calculated
# from all the functions in 'obs' and 'sim_m'.
#' @importFrom stats sd
curve_set_sd <- function(curve_set) {
  funcs <- data_and_sim_curves(curve_set)
  sapply(1:ncol(funcs), FUN = function(i) { sd(funcs[,i]) })
}

# A helper function to obtain the quantiles of functions in curve_set.
#
# No matter whether obs is a matrix or a vector, the quantiles are calculated
# from all the functions in 'obs' and 'sim_m'.
# @param ... Additional parameters passed to \code{\link[stats]{quantile}}.
#' @importFrom stats quantile
curve_set_quant <- function(curve_set, probs, ...) {
  funcs <- data_and_sim_curves(curve_set)
  # Dimensions: 2, r_idx.
  sapply(1:ncol(funcs), FUN = function(i) { quantile(funcs[,i], probs=probs, ...) })
}

# Number of arguments where the curves have been evaluated
curve_set_narg <- function(curve_set) {
  r <- curve_set[['r']]
  if(is.vector(r)) length(r)
  else nrow(r)
}

# Return curve_set$r as a data.frame
curve_set_rdf <- function(curve_set) {
  r <- curve_set[['r']]
  if(is.vector(r)) data.frame(r=r)
  else r
}
