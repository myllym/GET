# Turn an \code{envelope} object of \pkg{spatstat} into a curve_set object.
#
# @param env An \code{envelope} object of \pkg{spatstat}. The envelope()
#   functions must have been called with savefuns = TRUE.
# @return A corresponding curve_set object.
# @param ... Do not use. (For internal use only.)
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
# values in an \code{envelope} object of \pkg{spatstat} at the first place, if those are cropped
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
      stop('The names of curve_set[[\'r\']] must be either c("height", "width", "x", "y") or c("xmax", "xmin", "ymax", "ymin").')
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
# @return If an \code{envelope} object of \pkg{spatstat} was given, return a
#   corresponding curve_set object. If a curve_set object was given, return
#   it unharmed.
convert_envelope <- function(curve_set, ...) {
  if(inherits(curve_set, 'envelope')) {
    curve_set <- envelope_to_curve_set(curve_set, ...)
  }
  else if(!inherits(curve_set, 'curve_set')) {
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
convert_fdata <- function(curve_set, ...) {
  if(inherits(curve_set, 'fdata')) {
    curve_set <- fdata_to_curve_set(curve_set, ...)
  } else if(!inherits(curve_set, 'curve_set')) {
    stop('curve_set must either have class "fdata" (from fda.usc) ',
         'or class "curve_set".')
  }
  check_curve_set_content(curve_set, ...)
  curve_set
}

# Convert an envelope or fdata object to a curve_set object.
convert_to_curveset <- function(curve_set, ...) {
  if(inherits(curve_set, 'envelope')) {
    curve_set <- envelope_to_curve_set(curve_set, ...)
  } else if(inherits(curve_set, 'fdata')) {
    curve_set <- fdata_to_curve_set(curve_set, ...)
  } else if(!inherits(curve_set, 'curve_set')) {
    stop('curve_set must either have class "envelope" (from spatstat) ',
         'or "fdata" (from fda.usc) or "curve_set".')
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

#' Create a curve_set object
#'
#' Create a curve_set object out of a list in the right form.
#'
#'
#' The function is used to clump together the functional data in the form
#' that can be handled by the other \pkg{GET} functions (\code{\link{forder}},
#' \code{\link{central_region}}, \code{\link{global_envelope_test}} etc.).
#' The function \code{create_curve_set} takes care of checking the content of
#' the data, and saves relevant information of the curves for global envelope
#' methods to be used in particular for plotting the results with graphical
#' interpretation.
#'
#' \code{obs} must be either
#' \itemize{
#' \item a vector containing the data function/vector, or
#' \item a matrix containing the s data functions/vectors, in which case it is assumed that
#' each column corresponds to a data function/vector.
#' }
#'
#' If given, \code{r} describes the 1- or 2-dimensional argument values where the functions/vectors
#' have been observed (or simulated). It must be either
#' \itemize{
#' \item a vector,
#' \item a data.frame with columns "x", "y", "width" and "height",
#' where the width and height give the width and height of the pixels placed at x and y, or
#' \item a data.frame with columns "xmin", "xmax", "ymin" and "ymax" giving the corner
#' coordinates of the pixels where the data have been observed.
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
#' @return An object of class \code{curve_set} containing the data.
#' If the argument values are two-dimensional, then the \code{curve_set} is additionally
#' a \code{curve_set2d} object.
#' @export
#' @seealso \code{\link{plot.curve_set}}, \code{\link{plot.curve_set2d}}
#' @examples
#' # 1d
#' cset <- create_curve_set(list(r=1:10, obs=matrix(runif(10*5), ncol=5)))
#' plot(cset)
#' # 2d
#' cset <- create_curve_set(list(r=data.frame(x=c(rep(1:3, 3), 4), y=c(rep(1:3, each=3), 1),
#'                                            width=1, height=1),
#'                               obs=matrix(runif(10*5), ncol=5)))
#' plot(cset)
create_curve_set <- function(curve_set, ...) {
  check_curve_set_content(curve_set, ...)
  is1obs <- is.vector(curve_set[['obs']])
  if(is1obs) {
    funcs <- cbind(curve_set[['obs']], curve_set[['sim_m']])
    colnames(funcs) <- c('obs', paste("sim", 1:(ncol(funcs)-1), sep=""))
  }
  else funcs <- curve_set[['obs']]
  if(!is.null(curve_set[['r']])) r <- curve_set[['r']]
  else r <- 1:nrow(funcs)
  cset <- list(r = r, funcs = funcs, is1obs = is1obs)
  if(!is.null(curve_set[['theo']])) cset$theo <- curve_set[['theo']]
  class(cset) <- 'curve_set'
  if(!is.vector(r)) class(cset) <- c('curve_set2d', class(cset))
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
  dtext <- paste0("(", ifelse(curve_set_is1d(x), 1, 2), "d)")
  cat("A curve_set", dtext, " object with ", curve_set_nfunc(x), " curves observed at ",
      curve_set_narg(x), " argument values", sep="")
  if(curve_set_is1obs(x)) cat("\n(1 observed,", curve_set_nfunc(x)-1, "simulated)")
  cat(".\n")
  cat("Contains: \n")
  cat("$ r     : ")
  str(x$r)
  cat("$ funcs : ")
  str(x$funcs)
}

#' Plot method for the class 'curve_set'
#'
#' @param x An \code{curve_set} object.
#' @param idx Indices of functions to highlight with color \code{col_idx}.
#' Default to the observed function, if there is just one.
#' The legend of curves' colours is shown if indices are given or \code{x} contains one observed function.
#' See examples to remove the legend if desired.
#' @param col_idx A color for the curves to highlight, or a vector of the same length as \code{idx}
#' containing the colors for the highlighted functions. Default exists.
#' @param idx_name A variable name to be printed with the highlighted functions' idx. Default to empty.
#' @param col The basic color for the curves (which are not highlighted).
#' @param ... Ignored.
#' @seealso \code{\link{create_curve_set}}
#'
#' @export
#' @importFrom ggplot2 ggplot geom_line aes_ scale_color_manual labs
#' @importFrom viridisLite viridis
#' @examples
#' cset <- create_curve_set(list(r=1:10, obs=matrix(runif(10*5), ncol=5)))
#' plot(cset)
#' # Highlight some functions
#' plot(cset, idx=c(1,3))
#' plot(cset, idx=c(1,3), col_idx=c("black", "red"))
#' # Change legend
#' plot(cset, idx=c(1,3), col_idx=c("black", "red"), idx_name="Special functions")
#' plot(cset, idx=c(1,3)) + ggplot2::theme(legend.position = "bottom")
#' # Add labels
#' plot(cset, idx=c(1,3)) + ggplot2::labs(x="x", y="Value")
#' # and title
#' plot(cset) + ggplot2::labs(title="Example curves", x="x", y="Value")
#' # A curve_set with one observed function (other simulated)
#' if(requireNamespace("mvtnorm", quietly=TRUE)) {
#'   cset <- create_curve_set(list(obs=c(-1.6, 1.6),
#'             sim_m=t(mvtnorm::rmvnorm(200, c(0,0), matrix(c(1,0.5,0.5,1), 2, 2)))))
#'   plot(cset)
#'   # Remove legend
#'   plot(cset) + ggplot2::theme(legend.position = "none")
#' }
plot.curve_set <- function(x, idx, col_idx, idx_name = "", col = 'grey70', ...) {
  if(!all(x$r[-1] - x$r[-curve_set_narg(x)] > 0))
    warning("r values non-increasing. Plot not valid.")

  if(missing(idx)) {
    if(curve_set_is1obs(x))
      idx <- 1
    else
      idx <- NULL
  }
  else {
    if(!is.numeric(idx)) stop("idx should be numeric.")
  }
  if(missing(col_idx)) {
    if(length(idx) == 1)
      col_idx <- 1
    else if(length(idx) > 1)
      col_idx <- viridis(length(idx))
    else
      col_idx <- NULL
  }
  funcs <- curve_set_funcs(x)
  nfunc <- curve_set_nfunc(x)
  if(!is.null(idx)) {
    if(length(col_idx) == 1) col_idx <- rep(col_idx, times=length(idx))
    else if(length(col_idx) != length(idx))
      stop("If given, the length of col_idx should be 1 or equal to the length of idx.")
    # Labels for the highlighted functions
    if(!is.null(colnames(funcs))) idx_labs <- colnames(funcs)[idx]
    else idx_labs <- idx
    # Arrange the function such that the highlighted functions are the last ones
    # (ggplot plots them to the top)
    id_v <- idx_v <- 1:ncol(funcs)
    idx_v[idx] <- idx_labs
    if(length(idx) <= nfunc) {
      if(curve_set_is1obs(x) & length(idx) == 1) othername <- "sim"
      else othername <- "Other"
      idx_labs <- c(idx_labs, othername)
      idx_v[-idx] <- othername
    }
    id_v_levels <- c(id_v[-idx], idx)
  }
  else {
    id_v <- 1:ncol(funcs)
    id_v_levels <- id_v
  }

  df <- data.frame(r = x$r, funcs = c(funcs),
                   id =  factor(rep(id_v, each=nrow(funcs)), levels = id_v_levels))
  if(!is.null(idx)) df$idx = factor(rep(idx_v, each=nrow(funcs)), levels = idx_labs)

  p <- ( ggplot() )
  if(!is.null(idx)) {
    col_values <- c(col_idx, col)
    names(col_values) <- idx_labs
    p <- ( p + geom_line(data=df, aes_(x = ~r, y = ~funcs, group = ~id, col = ~idx))
           + scale_color_manual(values = col_values)
           + labs(col = idx_name) )
  }
  else {
    p <- ( p + geom_line(data=df, aes_(x = ~r, y = ~funcs, group = ~id), col = col) )
  }
  p
}

#' Plot method for the class 'curve_set2d'
#'
#' Plot method for the class 'curve_set2d', i.e. two-dimensional functions
#'
#' @param x An \code{curve_set2d} object
#' @param idx Indices of 2d functions to plot.
#' @param ... Ignored.
#' @inheritParams plot.combined_global_envelope
#'
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 labs
#' @examples
#' data("abide_9002_23")
#' plot(abide_9002_23$curve_set, idx=c(1, 27))
plot.curve_set2d <- function(x, idx=1, ncol = 2 + 1*(length(idx)==3), ...) {
  funcs <- curve_set_funcs(x)
  rdf <- curve_set_rdf(x)
  df <- do.call(rbind, lapply(idx, function(i) data.frame(rdf, idx=i, f=funcs[,i])))
  return(ggplot() + choose_geom(df, varfill='f') +
           facet_wrap("idx", ncol=ncol) + labs(x="x", y="y", fill=""))
}

# Combine curve sets.
#
# Combine curve sets to a one curve set, e.g. for testing by means of several test functions.
# @param x A list of curve sets or \code{envelope} (from \pkg{spatstat}) objects.
# @param equalr Whether to demand equal lengths of r vectors of the different curve sets
# @return A curve set that is a combination of the curve sets given in 'x'.
combine_curve_sets <- function(x, equalr = TRUE) {
  x <- check_curve_set_dimensions(x, equalr=equalr)
  cset <- x[[1]]
  # Combine
  cset$funcs <- do.call(rbind, lapply(x, FUN=function(curve_set) { curve_set[['funcs']] }))
  if(!is.null(x[[1]][['r']])) {
    if(is.vector(x[[1]]$r))
      cset$r <- do.call(c, lapply(x, FUN=function(curve_set) { curve_set[['r']] }))
    else
      cset$r <- do.call(rbind, lapply(x, FUN=function(curve_set) { curve_set[['r']] }))
  }
  if(!is.null(x[[1]][['theo']]))
    cset$theo <- do.call(c, lapply(x, FUN=function(curve_set) { curve_set[['theo']] }))
  # Return the combined set of curves
  cset
}

# Check curve set dimensions
#
# Check that x contains list of curve sets or \code{envelope} (from \pkg{spatstat}) objects.
# If the latter, then convert the objects to curve sets.
# Check that the curve sets have same elements and dimensions of them (numbers of r-values can differ for equalr=FALSE).
# @inheritParams combine_curve_sets
check_curve_set_dimensions <- function(x, equalr=FALSE) {
  x <- lapply(x, FUN=convert_to_curveset)
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
# A helper function to give the names of the functions
curve_set_funcnames <- function(curve_set) {
  colnames(curve_set[['funcs']])
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
  applyr(curve_set, median)
}

# A helper function to obtain the sd of all the functions in curve_set.
#' @importFrom stats sd
curve_set_sd <- function(curve_set) {
  applyr(curve_set, sd)
}

# A helper function to obtain the quantiles of all the functions in curve_set.
# @param ... Additional parameters passed to \code{\link[stats]{quantile}}.
#' @importFrom stats quantile
curve_set_quant <- function(curve_set, probs, ...) {
  applyr(curve_set, function(x) { quantile(x, probs=probs, ...) }) # Dimensions: 2, r_idx.
}

# Number of arguments where the curves have been evaluated
curve_set_narg <- function(curve_set) {
  nrow(curve_set[['funcs']])
}

# Number of curves
curve_set_nfunc <- function(curve_set) {
  ncol(curve_set[['funcs']])
}

# Return curve_set$r as a data.frame
curve_set_rdf <- function(curve_set) {
  r <- curve_set[['r']]
  if(is.vector(r)) data.frame(r=r)
  else r
}

curve_set_is1d <- function(curve_set) {
  !inherits(curve_set, 'curve_set2d')
}

curve_set_is2d <- function(curve_set) {
  inherits(curve_set, 'curve_set2d')
}

curve_set_isresidual <- function(curve_set) {
  x <- curve_set[['is_residual']]
  if(is.null(x)) FALSE
  else x
}

# Check whether there is also one observed function (and many simulated ones), i.e.
# the case of Monte Carlo and permutation tests.
curve_set_is1obs <- function(curve_set) {
  curve_set[['is1obs']]
}

curve_set_1obs <- function(curve_set) {
  if(curve_set_is1obs(curve_set)) curve_set[['funcs']][,1]
  else stop("Internal error.")
}

#' Subsetting curve sets
#'
#' Return subsets of curve sets which meet conditions.
#' @param x A \code{curve_set} object.
#' @param subset A logical expression indicating curves to keep.
#' @param ... Ignored.
#' @export
#' @examples
#' if(require("fda.usc", quietly=TRUE)) {
#'   # Prepare data
#'   data("poblenou")
#'   Free <- poblenou$df$day.festive == 1 |
#'     as.integer(poblenou$df$day.week) >= 6
#'   MonThu <- poblenou$df$day.festive == 0 & poblenou$df$day.week %in% 1:4
#'   Friday <- poblenou$df$day.festive == 0 & poblenou$df$day.week == 5
#'
#'   # Data as a curve_set
#'   cset <- create_curve_set(list(r=0:23,
#'              obs=t(poblenou[['nox']][['data']])))
#'   plot(subset(cset, MonThu))
#'   plot(subset(cset, Friday))
#'   plot(subset(cset, Free))
#' }
subset.curve_set <- function(x, subset, ...) {
  x[['funcs']] <- x[['funcs']][, subset]
  x
}

#' @export
`[.curve_set` <- function(x, i, j) {
  x$r <- if(is.data.frame(x$r)) x$r[i,]
    else x$r[i]
  x$funcs <- x$funcs[i,j,drop=FALSE]
  if(!missing(j) && !(1 %in% j))
    x$is1obs <- FALSE
  x
}
