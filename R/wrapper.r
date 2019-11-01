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
#' @param GET.args A named list of additional arguments to be passed to \code{\link{global_envelope_test}}.
#' @inheritParams graph.fanova
#' @importFrom plyr is.formula
#' @importFrom stats formula
#' @export
#' @examples
#' if(require("sp", quietly=TRUE) & require("gstat", quietly=TRUE)) {
#'   # Examples from gstat complemented with global envelopes
#'   #-------------------------------------------------------
#'   data(meuse)
#'   coordinates(meuse) <- ~x+y
#'   # topsoil zinc concentration, mg kg-1 soil ("ppm")
#'   bubble(meuse, "zinc",
#'          col=c("#00ff0088", "#00ff0088"), main="zinc concentrations (ppm)")
#'   # Variogram can be calculated as follows by the function variogram of the gstat library.
#'   # The function variogram takes a formula as its first argument:
#'   # log(zinc)~1 means that we assume a constant trend for the variable log(zinc).
#'   lzn.vgm <- variogram(object=log(zinc)~1, data=meuse)
#'   plot(lzn.vgm)
#'   # Variogram with global envelopes is as easy:
#'   \donttest{lzn.vgm.GET <- GET.variogram(object=log(zinc)~1, data=meuse)}
#'   \dontshow{lzn.vgm.GET <- GET.variogram(object=log(zinc)~1, data=meuse, nsim=4, GET.args=list(alpha=0.2))}
#'   plot(lzn.vgm.GET)
#'
#'   # Instead of the constant mean, denoted by ~1, a mean function can
#'   # be specified, e.g. using ~sqrt(dist) as a predictor variable:
#'   lznr.vgm <- variogram(log(zinc)~sqrt(dist), meuse)
#'   # In this case, the variogram of residuals with respect
#'   # to a fitted mean function are shown.
#'   plot(lznr.vgm)
#'   # The variogram with global envelopes (obtained by permuting the residuals):
#'   \donttest{lznr.vgm.GET <- GET.variogram(object=log(zinc)~sqrt(dist), data=meuse)}
#'   \dontshow{lznr.vgm.GET <- GET.variogram(object=log(zinc)~sqrt(dist), data=meuse, nsim=4, GET.args=list(alpha=0.2))}
#'   plot(lznr.vgm.GET)
#'
#'   # Directional variograms
#'   lzn.dir <- variogram(object=log(zinc)~1, data=meuse, alpha=c(0, 45, 90, 135))
#'   plot(lzn.dir)
#'   # with global envelopes
#'   \donttest{lzn.dir.GET <- GET.variogram(object=log(zinc)~1, data=meuse, alpha=c(0, 45, 90, 135))}
#'   \dontshow{lzn.dir.GET <- GET.variogram(object=log(zinc)~1, data=meuse, nsim=4, alpha=c(0, 45, 90, 135), GET.args=list(alpha=0.2))}
#'   plot(lzn.dir.GET, base_size=10)
#'
#'   # Use instead gstat objects
#'   g <- gstat(id="ln.zinc", formula=log(zinc)~1, data=meuse)
#'   # or: g <- gstat(id="ln.zinc", formula=log(zinc)~sqrt(dist), data=meuse)
#'   # The variogram
#'   plot(variogram(g))
#'   # The variogram with global envelopes:
#'   \donttest{g.GET <- GET.variogram(object=g)}
#'   \dontshow{g.GET <- GET.variogram(object=g, nsim=4, GET.args=list(alpha=0.2))}
#'   plot(g.GET)
#' }
GET.variogram <- function(object, nsim = 999, data = NULL, ..., GET.args = NULL, savefuns = TRUE) {
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

  obs.s <- split(obs, f=list(id=obs$id, dir.hor=obs$dir.hor))
  sim.s <- lapply(by(sim, INDICES=list(id=obs$id, dir.hor=obs$dir.hor), FUN=identity),
                  as.matrix)
  csets <- NULL
  for(i in 1:length(obs.s)) {
    csets[[names(obs.s)[[i]]]] <- create_curve_set(list(r = obs.s[[i]]$dist,
                                                        obs = obs.s[[i]]$gamma,
                                                        sim_m = sim.s[[i]]))
  }
  res <- do.call(global_envelope_test, c(list(curve_sets=csets, nstep=1), GET.args))

  attr(res, "xlab") <- attr(res, "xexp") <- "distance"
  attr(res, "ylab") <- attr(res, "yexp") <- attr(obs, "what")
  if(length(levels(obs$id)) == 1) labels <- ""
  else labels <- levels(obs$id)
  if(length(unique(obs$dir.hor)) > 1) {
    labels <- paste(labels, unique(obs$dir.hor), sep="")
  }
  attr(res, "labels") <- labels
  attr(res, "variogram") <- obs
  if(savefuns) attr(res, "curve_set") <- csets
  attr(res, "call") <- match.call()
  res
}


#-----------------------------------------------------------#
# n sample test of correspondence of distribution functions #
#-----------------------------------------------------------#

# Ecdf means
# y is a vector for which groups gives grouping
#' @importFrom stats ecdf
ecdfmeans.m <- function(x, groups, r) {
  ecdf.ls <- by(x, INDICES=groups, FUN=stats::ecdf, simplify=FALSE)
  sapply(ecdf.ls, FUN = function(x) { x(r) }, simplify=TRUE)
}

# Ecdf contrasts
# y is a vector for which groups gives grouping
#' @importFrom stats ecdf
ecdfcontrasts.m <- function(x, groups, r) {
  k <- nlevels(groups)
  gnames <- levels(groups)
  ecdf.ls <- by(x, INDICES=groups, FUN=stats::ecdf, simplify=FALSE)
  cont <- matrix(0, nrow=length(r), ncol=k*(k-1)/2)
  cont.names <- vector(length=k*(k-1)/2)
  counter <- 1
  for(i in 1:(k-1)) for(j in (i+1):k) {
    cont[, counter] <- ecdf.ls[[i]](r) - ecdf.ls[[j]](r)
    cont.names[counter] <- paste(gnames[i], gnames[j], sep="-")
    counter <- counter+1
  }
  colnames(cont) <- cont.names
  cont
}

#' Graphical n sample test of correspondence of distribution functions
#'
#' Compare the distributions of two (or more) groups.
#'
#'
#' A global envelope test can be performed to investigate whether the n distribution functions
#' differ from each other significally and how do they differ. This test is a generalization of
#' the two-sample Kolmogorov-Smirnov test with a graphical interpretation.
#' We assume that the curves in the different groups are an i.i.d. samples from the distribution
#' \deqn{F_i(r), i=1, \dots, n}{F_i(r), i=1, ..., n},
#' and we want to test the hypothesis
#' \deqn{F_1(r)= \dots = F_n(r)}{F_1(r)= ... = F_n(r)}.
#' If \code{contrasts = FALSE} (default), then the test statistic is taken to be
#' \deqn{\mathbf{T} = (\hat{F}_1(r), \dots, \hat{F}_n(r))}{T = (\hat{F}_1(r), \dots, \hat{F}_n(r))}
#' where \eqn{\hat{F}_i(r) = (\hat{F}_i(r_1), \dots, \hat{F}_i(r_k))}{\hat{F}_i(r) = (\hat{F}_i(r_1), ..., \hat{F}_i(r_k))}
#' is the ecdf of the $i$th sample evaluated at argument values
#' \eqn{r = (r_1,\dots,r_k)}{r = (r_1, ...,r_k)}.
#' This is our recommended test function for the test.
#' Another possibility is given by \code{contrasts = TRUE}, and then the test statistic is
#' \deqn{\mathbf{T} = (\hat{F}_1(r)-\hat{F}_2(r), \dots, \hat{F}_{n-1}(r)-\hat{F}_n(r))}{T = (\hat{F}_1(r)-\hat{F}_2(r), ..., \hat{F}_{n-1}(r)-\hat{F}_n(r))}
#'
#' The simulations under the null hypothesis that the distributions are the same can be obtained
#' by permuting the individuals of the groups. The default number of permutation, if nsim is not specified,
#' is n*1000 - 1 for the case \code{contrasts = FALSE} and
#' (n*(n-1)/2)*1000 - 1 for the case \code{contrasts = TRUE},
#' where n is the length of x.
#'
#' @param x A list (of length n) of values in the n groups.
#' @param r The sequence of argument values at which the distribution functions are compared.
#' The default is 100 equally spaced values between the minimum and maximum over all groups.
#' @inheritParams graph.fanova
#' @export
#' @examples
#' if(require(fda, quietly=TRUE)) {
#'   # Heights of boys and girls at age 10
#'   f.a <- growth$hgtf["10",] # girls at age 10
#'   m.a <- growth$hgtm["10",] # boys at age 10
#'   # Empirical cumulative distribution functions
#'   plot(ecdf(f.a))
#'   plot(ecdf(m.a), col=grey(0.7), add=TRUE)
#'   # Create a list of the data
#'   fm.list <- list(Girls=f.a, Boys=m.a)
#'   \donttest{
#'   res_m <- GET.necdf(fm.list)
#'   plot(res_m)
#'   res_c <- GET.necdf(fm.list, contrasts = TRUE)
#'   plot(res_c)
#'   }
#'   \dontshow{
#'   # The test with lower number of simulations
#'   res_m <- GET.necdf(fm.list, nsim=4, alpha=0.2)
#'   plot(res_m)
#'   res_c <- GET.necdf(fm.list, contrasts = TRUE, nsim=4, alpha=0.2)
#'   plot(res_c)
#'   }
#'
#'   # Heights of boys and girls at age 14
#'   f.a <- growth$hgtf["14",] # girls at age 14
#'   m.a <- growth$hgtm["14",] # boys at age 14
#'   # Empirical cumulative distribution functions
#'   plot(ecdf(f.a))
#'   plot(ecdf(m.a), col=grey(0.7), add=TRUE)
#'   # Create a list of the data
#'   fm.list <- list(Girls=f.a, Boys=m.a)
#'   \donttest{
#'   res_m <- GET.necdf(fm.list)
#'   plot(res_m)
#'   res_c <- GET.necdf(fm.list, contrasts = TRUE)
#'   plot(res_c)
#'   }
#'   \dontshow{
#'   # The test with lower number of simulations
#'   res_m <- GET.necdf(fm.list, nsim=4, alpha=0.2)
#'   plot(res_m)
#'   res_c <- GET.necdf(fm.list, contrasts = TRUE, nsim=4, alpha=0.2)
#'   plot(res_c)
#'   }
#' }
GET.necdf <- function(x, r = seq(min(unlist((lapply(x, min)))), max(unlist((lapply(x, max)))), length=100),
                      contrasts = FALSE, nsim, ...) {
  if(!is.list(x) && length(x)<2) stop("At least two groups should be provided.\n")
  x.lengths <- as.numeric(lapply(x, FUN = length))
  if(!is.null(names(x))) groups <- rep(names(x), times=x.lengths)
  else groups <- rep(1:length(x), times=x.lengths)
  groups <- factor(groups, levels=unique(groups))
  gnames <- levels(groups)
  if(missing(nsim)) {
    if(!contrasts) {
      nsim <- length(x)*1000 - 1
    }
    else {
      J <- length(x)
      nsim <- (J*(J-1)/2)*1000 - 1
    }
    message("Creating ", nsim, " permutations.\n", sep="")
  }
  # setting the 'fun', "means" or "contrasts"
  if(!contrasts) fun <- ecdfmeans.m
  else fun <- ecdfcontrasts.m
  x <- unlist(x)
  # Observed difference between the ecdfs
  obs <- fun(x, groups, r)
  # Simulations by permuting to which groups each value belongs to
  sim <- replicate(nsim, fun(x, sample(groups, size=length(groups), replace=FALSE), r), simplify = "array")
  complabels <- colnames(obs)

  csets <- list()
  for(i in 1:ncol(obs)) {
    csets[[i]] <- create_curve_set(list(r = r,
                                        obs = obs[,i],
                                        sim_m = sim[,i,]))
  }
  names(csets) <- complabels
  # GET
  res <- global_envelope_test(csets, alternative="two.sided", ..., nstep=1)
  attr(res, "xlab") <- "x"
  attr(res, "xexp") <- quote(x)
  if(!contrasts) {
    attr(res, "ylab") <- "F(x)"
    attr(res, "yexp") <- quote(hat(F)(x))
  }
  else {
    attr(res, "ylab") <- "diff. F(x)"
    attr(res, "yexp") <- quote(hat(F)[i](x)-hat(F)[j](x))
  }
  attr(res, "contrasts") <- contrasts
  attr(res, "labels") <- complabels
  attr(res, "call") <- match.call()
  res
}
