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
  simfuns_fv <- attr(env, 'simfuns', exact = TRUE)
  if(is.null(simfuns_fv)) {
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

  theo <- env[['theo']]

  res <- list(r = r, obs = env[['obs']], sim_m = sim_m)
  if(!is.null(theo)) {
    res[['theo']] <- theo
  }

  create_curve_set(res, ...)
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
  r <- fdata[['argvals']]
  if(is.list(r))
    stop("No support for converting this type of functional data, only functional data observed in a single grid supported.")

  create_curve_set(list(r = r, obs = t(fdata[['data']])), ...)
}

# Check the content validity of a potential curve_set object.
#
# @param curve_set An object to be checked.
# @param allow_Inf_values Logical. Can be used to allow infinite or nonnumeric
# values in an \code{\link[spatstat]{envelope}} object at the first place, if those are cropped
# away (in \code{\link{crop_curves}}).
check_curve_set_content <- function(curve_set, allow_Inf_values = FALSE) {
  if(inherits(curve_set, "curve_set")) {
    if(!allow_Inf_values) {
      if(!all(is.finite(curve_set[['funcs']]))) {
        stop('All curves must have only finite values.')
      }
      if(!all(is.finite(curve_set[['theo']]))) {
        stop('The theoretical values should be finite.')
      }
    }
    return(curve_set)
  }

  possible_names <- c('r', 'obs', 'sim_m', 'theo')

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

  obs <- curve_set[['obs']]
  if(is.vector(obs)) narg <- length(obs)
  else if(is.matrix(obs)) narg <- nrow(obs)
  else {
    stop('curve_set[["obs"]] must be a vector or a matrix.')
  }
  if(!(allow_Inf_values | ( all(is.numeric(obs)) && all(is.finite(obs)) ))) {
    stop('curve_set[["obs"]] must have only finite numeric values.')
  }

  r <- curve_set[['r']]
  if(is.null(r)) n_r <- narg
  else if(is.data.frame(r)) {
    if(!(identical(sort(names(r)), c("height", "width", "x", "y"))
       || identical(sort(names(r)), c("xmax", "xmin", "ymax", "ymin"))))
      stop('The names of curve_set[[\'r\']] must be either c("height", "width", "x", "y") or c("xmax", "xmin", "ymax", "ymin").\n')
    n_r <- nrow(r)
    if(!all(sapply(r, is.numeric)) || !all(sapply(r, is.finite))) {
      stop('curve_set[["r"]] must have only finite numeric values.')
    }
  }
  else if(is.vector(r)) {
    if(!all(is.numeric(r)) || !all(is.finite(r))) {
      stop('curve_set[["r"]] must have only finite numeric values.')
    }
    n_r <- length(r)
  }
  else {
    stop("If r given, it must be either a vector or a data.frame.")
  }
  if(n_r != narg) {
    stop('curve_set[["r"]] must be compatible with curve_set[["obs"]].')
  }

  sim_m <- curve_set[['sim_m']]
  if(!is.null(sim_m)) {
    dim_sim_m <- dim(sim_m)
    if(!is.matrix(sim_m)) {
      stop('curve_set[["sim_m"]] must be a matrix.')
    }
    if(dim_sim_m[1] != narg) {
      stop('curve_set[["sim_m"]] must have as many rows as there are ',
           'elements in curve_set[["obs"]].')
    }
    if(dim_sim_m[2] < 1L) {
      stop('curve_set[["sim_m"]] must have at least one column.')
    }
    if(!(allow_Inf_values | ( is.numeric(sim_m) && all(is.finite(sim_m)) ))) {
      stop('curve_set[["sim_m"]] must have only finite numeric values.')
    }
  }

  theo <- curve_set[['theo']]
  if(!is.null(theo)) {
    n_theo <- length(theo)
    if(is.null(sim_m)) stop("theo can be provided only with sim_m.")
    if(n_theo != narg) {
      stop('curve_set[["theo"]] must have as many values as ',
           'curve_set[["obs"]].')
    }
    if(!is.vector(theo)) {
      stop('curve_set[["theo"]] must be a vector.')
    }
    if(!(allow_Inf_values | ( is.numeric(theo) && all(is.finite(theo)) ))) {
      stop('curve_set[["theo"]] must have only finite numeric ',
           'values.')
    }
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
  if(!curve_set_isresidual(curve_set)) {
    stop('curve_set must consist of residual curves. Run function ',
         'residual first.')
  }
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
#' @param ... For expert use only.
#' @return An object of class \code{curve_set}.
#' @export
create_curve_set <- function(curve_set, ...) {
  check_curve_set_content(curve_set, ...)
  cset <- list()
  if(!is.null(curve_set[['r']])) cset$r <- curve_set[['r']]
  is1obs <- is.vector(curve_set[['obs']])
  if(is1obs) {
    cset$funcs <- cbind(curve_set[['obs']], curve_set[['sim_m']])
    colnames(cset$funcs) <- c('obs', paste("sim", 1:(ncol(cset$funcs)-1), sep=""))
  }
  else cset$funcs <- curve_set[['obs']]
  cset$is1obs <- is1obs
  if(!is.null(curve_set[['theo']])) cset$theo <- curve_set[['theo']]
  class(cset) <- 'curve_set'
  cset
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
  funcs <- curve_set_funcs(x)
  if(!curve_set_is1obs(x)) col_sim <- col_obs
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
  x <- check_curve_set_dimensions(x, equalr=equalr)
  cset <- x[[1]]
  # Combine
  cset$funcs <- do.call(rbind, lapply(x, FUN=function(curve_set) { curve_set[['funcs']] }))
  if(!is.null(x[[1]][['r']])) {
    if(is.vector(x[[1]]$r)) cset$r <- do.call(c, lapply(x, FUN=function(curve_set) { curve_set[['r']] }))
    else cset$r <- do.call(rbind, lapply(x, FUN=function(curve_set) { curve_set[['r']] }))
  }
  if(!is.null(x[[1]][['theo']]))
    cset$theo <- do.call(c, lapply(x, FUN=function(curve_set) { curve_set[['theo']] }))
  # Return the combined set of curves
  cset
}

# Check curve set dimensions
#
# Check that x contains list of curve sets or \code{\link[spatstat]{envelope}} objects.
# If the latter, then convert the objects to curve sets.
# Check that the curve sets have same elements and dimensions of them (numbers of r-values can differ for equalr=FALSE).
# @inheritParams combine_curve_sets
check_curve_set_dimensions <- function(x, equalr=FALSE) {
  x <- lapply(x, FUN=convert_envelope)
  checkequal <- function(f) {
    all(sapply(x, FUN=function(curve_set) { f(curve_set) == f(x[[1]]) }))
  }
  if(equalr) {
    if(!checkequal(curve_set_narg))
      stop("The lengths of 'r' in curve sets differ.")
  }
  if(!checkequal(curve_set_isresidual))
    stop("The element \'is_residual\' should be the same for each curve set.")
  if(!checkequal(curve_set_nfunc))
    stop("The numbers of functions in curve sets differ.")
  if(!checkequal(curve_set_is1obs))
    stop("The numbers of observed functions in curve sets differ.")
  x
}

# A helper function to return all the functions from a curve set in a matrix.
# Each row corresponds to a function.
data_and_sim_curves <- function(curve_set) {
  t(curve_set[['funcs']])
}
# A helper function to return all the functions from a curve set in a matrix.
# Each column corresponds to a function.
curve_set_funcs <- function(curve_set) {
  curve_set[['funcs']]
}

# A helper function to obtain the mean of functions in curve_set.
#
# Calculate the mean of all the functions in the curve_set.
curve_set_mean <- function(curve_set) {
  rowMeans(curve_set[['funcs']])
}


# A helper function to apply the function f to the curve_set functions at all argument values
applyr <- function(curve_set, f) {
  funcs <- curve_set[['funcs']]
  sapply(1:nrow(funcs), function(i) { f(funcs[i,]) })
}

# A helper function to obtain the median of all the functions in curve_set.
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

# Number of curves
curve_set_nfunc <- function(curve_set) {
  if(is.matrix(curve_set$obs)) ncol(curve_set[['obs']])
  else ncol(curve_set[['sim_m']]) + 1
}

# Return curve_set$r as a data.frame
curve_set_rdf <- function(curve_set) {
  r <- curve_set[['r']]
  if(is.vector(r)) data.frame(r=r)
  else r
}

curve_set_isresidual <- function(curve_set) {
  x <- curve_set[['is_residual']]
  if(is.null(x)) FALSE
  else x
}

# Check whether there is also one observed function (and many simulated ones), i.e.
# the case of Monte Carlo and permutation tests.
curve_set_is1obs <- function(curve_set) {
  is.vector(curve_set[['obs']])
}

curve_set_1obs <- function(curve_set) {
  if(curve_set_is1obs(curve_set)) curve_set[['obs']]
  else stop("Internal error.")
}
