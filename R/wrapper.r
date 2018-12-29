#--------------------------------------------------------#
# Variogram and residual variogram with global envelopes #
#--------------------------------------------------------#

# A helper function for permuting the data variables
# @inheritParams GET.variogram
#' @importFrom sp coordinates
#' @importFrom gstat variogram
#' @importFrom stats residuals
#' @importFrom stats lm
#' @importFrom stats na.exclude
permvariogram <- function(object, data, vars, perm=TRUE, ...) {
  if(length(vars)>1) stop("Only one variable allowed. No test for correlation between variables implemented.\n")
  args <- list(...)
  # Treat coordinates/locations in order to do permutations of the data
  if(class(data)[1] == "SpatialPointsDataFrame") locations <- sp::coordinates(data)
  else {
    if(!("locations" %in% args)) stop("Either data must be provided with coordinates or locations must be given separately. See ?variogram.\n")
  }
  if(!is.data.frame(data)) data <- as.data.frame(data)

  permdata <- data
  if(perm) {
    newids <- sample(1:nrow(data), size=nrow(data), replace=FALSE)
    permdata[, vars] <- data[newids, vars]
  }
  sp::coordinates(permdata) <- locations
  if(inherits(object, "gstat")) {
    for(i in seq(along.with = object$data)) object$data[[i]]$data <- permdata
    v <- gstat::variogram(object, ...)
  }
  else {
    v <- gstat::variogram(object, data=permdata, ...)
  }
  v
}

#' Variogram and residual variogram with global envelopes
#'
#' The function accompanies the function \code{\link[gstat]{variogram}} with global envelopes
#' that are based on permutations of the variable(s) or residuals for which the variogram is calculated.
#' Therefore, one can inspect the hypothesis of "no spatial autocorrelation" of the variable or the residuals
#' of the fitted model.
#'
#' @param object An object of class \code{gstat} or a \code{variogram.formula}.
#' In the first case, direct (residual) variograms are calculated for the variable
#' defined in object. Only one variable allowed.
#' In the second case, a formula defining the response vector and (possible) regressors,
#' in case of absence of regressors, use e.g. z~1. See \code{\link[gstat]{variogram}}.
#' @param nsim The number of permutations.
#' @param data A data frame where the names in formula are to be found. If NULL,
#' the data are assumed to be found in the \code{object}.
#' @param ... Additional parameters to be passed to \code{\link[gstat]{variogram}}.
#' @param GET.param A list of parameters to be passed to \code{global_envelope_test}. If NULL,
#' the default global envelope test is used.
#' @inheritParams graph.fanova
#' @importFrom plyr is.formula
#' @importFrom stats formula
#' @export
#' @examples
#' if(require(sp, quietly=TRUE) & require(gstat, quietly=TRUE)) {
#'   # Examples from gstat complemented with global envelopes
#'   #-------------------------------------------------------
#'   data(meuse)
#'   sp::coordinates(meuse) <- ~x+y
#'   # topsoil zinc concentration, mg kg-1 soil ("ppm")
#'   bubble(meuse, "zinc",
#'          col=c("#00ff0088", "#00ff0088"), main = "zinc concentrations (ppm)")
#'   # Variogram can be calculated as follows by the function variogram of the gstat library.
#'   # The function variogram takes a formula as its first argument:
#'   # log(zinc)~1 means that we assume a constant trend for the variable log(zinc).
#'   lzn.vgm <- variogram(object=log(zinc)~1, data=meuse)
#'   plot(lzn.vgm)
#'   # Variogram with global envelopes is as easy:
#'   lzn.vgm.GET <- GET.variogram(object=log(zinc)~1, data=meuse)
#'   plot(lzn.vgm.GET)
#'
#'   # Instead of the constant mean, denoted by ~1, a mean function can
#'   # be specified, e.g. using ~sqrt(dist) as a predictor variable:
#'   lznr.vgm <- variogram(log(zinc)~sqrt(dist), meuse)
#'   # In this case, the variogram of residuals with respect
#'   # to a fitted mean function are shown.
#'   plot(lznr.vgm)
#'   # The variogram with global envelopes (obtained by permuting the residuals):
#'   lznr.vgm.GET <- GET.variogram(object=log(zinc)~sqrt(dist), data=meuse)
#'   plot(lznr.vgm.GET)
#'
#'   # Directional variograms
#'   lzn.dir = variogram(object=log(zinc)~1, data=meuse, alpha=c(0, 45, 90, 135))
#'   plot(lzn.dir)
#'   # with global envelopes
#'   lzn.dir.GET <- GET.variogram(object=log(zinc)~1, data=meuse, alpha=c(0, 45, 90, 135))
#'   plot(lzn.dir.GET, plot_style="ggplot2", base_size=10)
#'
#'   # Use instead gstat objects
#'   g <- gstat(id = "ln.zinc", formula = log(zinc)~1, data = meuse)
#'   # or: g <- gstat(id = "ln.zinc", formula = log(zinc)~sqrt(dist), data = meuse)
#'   # The variogram
#'   plot(variogram(g))
#'   # The variogram with global envelopes:
#'   g.GET <- GET.variogram(object = g)
#'   plot(g.GET)
#' }
GET.variogram <- function(object, nsim = 999, data = NULL, ..., GET.param = NULL, savefuns = TRUE) {
  if(!inherits(object, "formula") & !inherits(object, "gstat")) stop("object does not have the formula or gstat class.\n")
  # Check data w.r.t. formula
  if(inherits(object, "formula")) {
    if(is.null(data)) stop("The argument \'data\' must be provided as \'object\' is a formula.\n")
    vars <- all.vars(object)
  }
  if(inherits(object, "gstat")) {
    if(length(object$data) > 1) stop("Only one dependent variable allowed.\n")
    if(is.null(data)) data <- object$data[[1]]$data
    vars <- all.vars(object$data[[1]]$formula)
  }
  if(!all(vars %in% names(data))) stop("The variables found in the given gstat object do not exist in the first \"data\" component.\n")
  # The case of regression model(s):
  if(inherits(object, "formula") & object[[3]] != 1) { # there is a regression model:
    data$resid <- stats::residuals(stats::lm(formula=object, data=data, na.action = stats::na.exclude))
    object <- resid ~ 1
    vars <- all.vars(object) # Update formula variables
  }
  if(inherits(object, "gstat")) {
    if(object$data[[1]]$formula[[3]] != 1) { # there is a regression model:
      data$resid <- stats::residuals(stats::lm(formula=object$data[[1]]$formula,
                                               data=data,
                                               na.action = stats::na.exclude))
      object$data[[1]]$formula <- stats::formula(paste("resid ~ 1", sep=""))
    }
    # Update formula variables
    vars <- all.vars(object$data[[1]]$formula)
  }
  # Calculate variograms for data and permutations
  obs <- permvariogram(object=object, data=data, vars=vars, perm=FALSE, ...)
  fun <- function(i, ...) {
    permvariogram(object=object, data=data, vars=vars, perm=TRUE, ...)$gamma
  }
  sim <- sapply(1:nsim, FUN = fun, ..., simplify = TRUE)
  cset <- create_curve_set(list(r = obs$dist,
                                obs = obs$gamma,
                                sim_m = sim))
  if(!is.null(GET.param)) {
    tmpstring <- paste("res <- global_envelope_test(cset", sep="")
    for(i in 1:length(GET.param))
      tmpstring <- paste(tmpstring, paste(", ", names(GET.param)[i], "=GET.param[[", i, "]]", sep=""), sep="")
    tmpstring <- paste(tmpstring, ")", sep="")
    eval(parse(text = tmpstring))
  }
  else res <- global_envelope_test(cset)
  attr(res, "xlab") <- attr(res, "xexp") <- "distance"
  attr(res, "ylab") <- attr(res, "yexp") <- attr(obs, "what")
  if(length(levels(obs$id)) == 1) labels <- ""
  else labels <- levels(obs$id)
  if(length(unique(obs$dir.hor)) > 1) {
    labels <- paste(labels, unique(obs$dir.hor), sep="")
  }
  attr(res, "labels") <- labels
  attr(res, "variogram") <- obs
  if(savefuns) attr(res, "curve_set") <- cset
  attr(res, "call") <- match.call()
  res
}