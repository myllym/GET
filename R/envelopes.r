get_alternative <- function(global_envelope) {
  attr(global_envelope, "alternative")
}

# It should be:
# small_significant=TRUE for 'rank', 'erl', 'cont' and 'area' -> ordering decreasing
# small_significant=FALSE for 'qdir', 'st', 'unscaled' -> ordering increasing (decreasing=FALSE)
critical <- function(distance, alpha, Nfunc, small_significant) {
  distancesorted <- sort(distance, decreasing=small_significant)
  distancesorted[floor((1-alpha)*Nfunc)]
}

make_envelope_object <- function(type, curve_set, LB, UB, T_0,
                                 picked_attr, isenvelope,
                                 Malpha, alpha, distance) {
  Nfunc <- curve_set_nfunc(curve_set)
  if(curve_set_is1obs(curve_set)) {
    df <- data.frame(curve_set_rdf(curve_set), obs=curve_set_1obs(curve_set),
                     central=T_0, lo=LB, hi=UB)
  }
  else {
    df <- data.frame(curve_set_rdf(curve_set), central=T_0, lo=LB, hi=UB)
  }
  if(isenvelope) {
    res <- spatstat::fv(x=df, argu = picked_attr[['argu']],
                        ylab = picked_attr[['ylab']], valu = "central", fmla = ". ~ r",
                        alim = c(min(curve_set[['r']]), max(curve_set[['r']])),
                        labl = picked_attr[['labl']], desc = picked_attr[['desc']],
                        unitname = NULL, fname = picked_attr[['fname']], yexp = picked_attr[['yexp']])
    attr(res, "shade") <- c("lo", "hi")
    attr(res, "alternative") <- picked_attr[['alternative']]
  }
  else {
    res <- df
    attrnames <- names(picked_attr)
    for(n in attrnames) attr(res, n) <- picked_attr[[n]]
  }
  class(res) <- c("global_envelope", class(res))
  if(curve_set_is2d(curve_set)) class(res) <- c("global_envelope2d", class(res))
  attr(res, "method") <- "Global envelope"
  attr(res, "type") <- type
  attr(res, "alpha") <- alpha
  attr(res, "M") <- distance
  attr(res, "M_alpha") <- Malpha
  res
}

# Functionality for central regions based on a curve set
# @param ... Ignored.
individual_central_region <- function(curve_set, type = "erl", coverage = 0.50,
                                      alternative = c("two.sided", "less", "greater"),
                                      probs = c((1-coverage)/2, 1-(1-coverage)/2),
                                      quantile.type = 7,
                                      central = "median") {
  isenvelope <- inherits(curve_set, "envelope")
  if(!is.numeric(coverage) || (coverage < 0 | coverage > 1)) stop("Unreasonable value of coverage.")
  alpha <- 1 - coverage
  if(!(type %in% c("rank", "erl", "cont", "area", "qdir", "st", "unscaled")))
    stop("No such type for global envelope.")
  alternative <- match.arg(alternative)
  small_significant <- TRUE
  if(type %in% c("qdir", "st", "unscaled")) {
    small_significant <- FALSE
    if(alternative != "two.sided") {
      warning("For qdir, st and unscaled envelopes only the two.sided alternative is valid.")
      alternative <- "two.sided"
    }
  }
  check_probs(probs)
  if(!(central %in% c("mean", "median"))) {
    central <- "median"
    warning("Invalid option fiven for central. Using central = median.")
  }
  picked_attr <- pick_attributes(curve_set, alternative=alternative) # saving for attributes / plotting purposes
  curve_set <- convert_to_curveset(curve_set)

  # Measures for functional ordering
  measure <- type
  scaling <- ""
  switch(type,
         qdir = {
           measure <- "max"
           scaling <- "qdir"
         },
         st = {
           measure <- "max"
           scaling <- "st"
         },
         unscaled = {
           measure <- "max"
           scaling <- "none"
         })
  distance <- forder(curve_set, measure=measure, scaling=scaling,
                     alternative=alternative, probs=probs, quantile.type=quantile.type)

  all_curves <- data_and_sim_curves(curve_set) # all the functions
  Nfunc <- length(distance) # Number of functions
  nr <- curve_set_narg(curve_set)
  # Define the central curve T_0
  T_0 <- get_T_0(curve_set)

  # Check reasonability of Nfunc vs alpha
  if(Nfunc*alpha < 1-.Machine$double.eps^0.5)
    stop("Number of functions s is only ", Nfunc, ", but alpha is ", alpha,
         ". So, s*alpha is ", Nfunc*alpha, ".", sep="")

  # The critical value
  Malpha <- critical(distance, alpha, Nfunc, small_significant)

  #-- 100(1-alpha)% global envelope
  switch(type,
         rank = {
           LB <- array(0, nr)
           UB <- array(0, nr)
           for(i in 1:nr){
             Hod <- sort(all_curves[,i])
             LB[i]<- Hod[Malpha]
             UB[i]<- Hod[Nfunc-Malpha+1]
           }
         },
         erl =,
         cont =,
         area = {
           j <- distance >= Malpha
           LB <- array(0, nr)
           UB <- array(0, nr)
           for(i in 1:nr){
             lu <- range(all_curves[j,i])
             LB[i]<- lu[1]
             UB[i]<- lu[2]
           }
         },
         qdir = {
           curve_set_res <- residual(curve_set, use_theo=TRUE)
           quant_m <- curve_set_quant(curve_set_res, probs=probs, type=quantile.type)
           LB <- T_0 - Malpha*abs(quant_m[1,])
           UB <- T_0 + Malpha*abs(quant_m[2,])
         },
         st = {
           sdX <- curve_set_sd(curve_set)
           LB <- T_0 - Malpha*sdX
           UB <- T_0 + Malpha*sdX
         },
         unscaled = {
           LB <- T_0 - Malpha
           UB <- T_0 + Malpha
         })

  switch(alternative,
         "two.sided" = {},
         "less" = { UB <- Inf },
         "greater" = { LB <- -Inf })

  res <- make_envelope_object(type, curve_set, LB, UB, T_0,
                              picked_attr, isenvelope,
                              Malpha, alpha, distance)
  attr(res, "call") <- match.call()
  res
}

# Functionality for global envelope tests based on a curve set (individual central region + p-values)
individual_global_envelope_test <- function(curve_set, type = "erl", alpha = 0.05,
                                            alternative = c("two.sided", "less", "greater"),
                                            ties = "erl",
                                            probs = c(0.025, 0.975), quantile.type = 7,
                                            central = "mean") {
  alternative <- match.arg(alternative)
  tmp <- convert_to_curveset(curve_set)
  if(!curve_set_is1obs(tmp))
    stop("The curve_set does not contain one observed function. Testing does not make sense.\n Did you want to construct a central region of your data? See the function central_region.")
  if(!is.numeric(alpha) || (alpha < 0 | alpha > 1)) stop("Unreasonable value of alpha.")
  res <- individual_central_region(curve_set, type=type, coverage=1-alpha,
                                   alternative=alternative,
                                   probs=probs, quantile.type=quantile.type,
                                   central=central)
  # The type of the p-value
  possible_ties <- c('midrank', 'random', 'conservative', 'liberal', 'erl')
  if(!(ties %in% possible_ties)) stop("Unreasonable ties argument!")

  # Measures for functional ordering
  distance <- attr(res, "M")

  #-- Calculate the p-values
  switch(type,
         rank = {
           u <- -distance
           #-- p-interval
           p_low <- estimate_p_value(x=u[1], sim_vec=u[-1], ties='liberal')
           p_upp <- estimate_p_value(x=u[1], sim_vec=u[-1], ties='conservative')
           #-- unique p-value
           if(ties == "erl") {
             distance_lexo <- forder(curve_set, measure="erl", alternative=alternative)
             u_lexo <- -distance_lexo
             p <- estimate_p_value(x=u_lexo[1], sim_vec=u_lexo[-1], ties="conservative")
           }
           else p <- estimate_p_value(x=u[1], sim_vec=u[-1], ties=ties)
         },
         erl = {
           u_lexo <- -distance
           p <- estimate_p_value(x=u_lexo[1], sim_vec=u_lexo[-1], ties="conservative")
         },
         cont = {
           u_cont <- -distance
           p <- estimate_p_value(x=u_cont[1], sim_vec=u_cont[-1], ties="conservative")
         },
         area = {
           u_area <- -distance
           p <- estimate_p_value(x=u_area[1], sim_vec=u_area[-1], ties="conservative")
         },
         qdir = {
           p <- estimate_p_value(x=distance[1], sim_vec=distance[-1])
         },
         st = {
           p <- estimate_p_value(x=distance[1], sim_vec=distance[-1])
         },
         unscaled = {
           p <- estimate_p_value(x=distance[1], sim_vec=distance[-1])
         })

  # Change the "method" attribute
  attr(res, "method") <- paste(attr(res, "method"), " test", sep="")
  # Add attributes related to p-values
  attr(res, "p") <- p
  if(type == "rank") {
    attr(res, "p_interval") <- c(p_low, p_upp)
    attr(res, "ties") <- ties
  }
  # Update "call" attribute
  attr(res, "call") <- match.call()
  res
}

# Functionality for combined_central_region and combined_global_envelope_test (two-step procedure)
combined_CR_or_GET <- function(curve_sets, CR_or_GET = c("CR", "GET"), coverage, ...) {
  ntests <- length(curve_sets)
  if(ntests < 1) stop("Only one curve_set, no combining to be done.")
  check_curve_set_dimensions(curve_sets) # Do not catch the curve_set here
  CR_or_GET <- match.arg(CR_or_GET)

  # 1) First stage: Calculate the functional orderings individually for each curve_set
  res_ls <- lapply(curve_sets, FUN = function(x) { individual_central_region(x, ...) })
  type <- attr(res_ls[[1]], "type")

  # 2) Second stage: ERL central region/test
  # Create a curve_set for the ERL test
  k_ls <- lapply(res_ls, FUN = function(x) attr(x, "M"))
  k_mat <- do.call(rbind, k_ls, quote=FALSE)
  Nfunc <- ncol(k_mat)
  # Construct the one-sided ERL central region
  if(type %in% c("qdir", "st", "unscaled")) alt2 <- "greater"
  else alt2 <- "less"
  switch(CR_or_GET,
         CR = {
           curve_set_u <- create_curve_set(list(r=1:ntests, obs=k_mat))
           res_erl <- individual_central_region(curve_set_u, type="erl", coverage=coverage, alternative=alt2)
         },
         GET = {
           curve_set_u <- create_curve_set(list(r=1:ntests, obs=k_mat[,1], sim_m=k_mat[,-1]))
           res_erl <- individual_global_envelope_test(curve_set_u, type="erl", alpha=1-coverage, alternative=alt2)
         }
  )
  res_erl <- envelope_set_labs(res_erl, xlab="Function", ylab="ERL measure")
  attr(res_erl, "labels") <- names(curve_sets)

  # 3) The 100(1-alpha)% global combined ERL envelope
  distance_lexo_sorted <- sort(attr(res_erl, "M"), decreasing=TRUE)
  Malpha <- distance_lexo_sorted[floor(coverage*Nfunc)]
  # Indices of the curves from which to calculate the convex hull
  curves_for_envelope_ind <- which(attr(res_erl, "M") >= Malpha)
  # Curves
  curve_sets <- lapply(curve_sets, FUN=convert_to_curveset)
  all_curves_l <- lapply(curve_sets, function(x) { data_and_sim_curves(x) })
  # Curves from which to calculate the convex hull
  curves_for_envelope_l <- lapply(all_curves_l, function(x) { x[curves_for_envelope_ind,] })
  # Bounding curves
  LB <- lapply(curves_for_envelope_l, FUN = function(x) { apply(x, MARGIN=2, FUN=min) })
  UB <- lapply(curves_for_envelope_l, FUN = function(x) { apply(x, MARGIN=2, FUN=max) })
  # Update the bounding curves (lo, hi) and Malpha to the first level central regions
  for(i in 1:ntests) {
    if(get_alternative(res_ls[[i]]) != "greater") res_ls[[i]]$lo <- LB[[i]]
    if(get_alternative(res_ls[[i]]) != "less") res_ls[[i]]$hi <- UB[[i]]
    attr(res_ls[[i]], "alpha") <- attr(res_ls[[i]], "M_alpha") <- NULL
    attr(res_ls[[i]], "method") <- paste0("1/", ntests, "th of a combined global envelope test")
  }
  if(!is.null(names(curve_sets))) names(res_ls) <- names(curve_sets)

  # Return
  attr(res_ls, "level2_ge") <- res_erl
  attr(res_ls, "level2_curve_set") <- curve_set_u
  switch(CR_or_GET,
         CR = {
           attr(res_ls, "method") <- "Combined global envelope"
         },
         GET = {
           attr(res_ls, "method") <- "Combined global test"
         })
  if(!is.null(attr(res_ls[[1]], "argu")))
    res_ls <- envelope_set_labs(res_ls, xlab=attr(res_ls[[1]], "xlab"),
                                ylab=substitute(italic(T(i)), list(i=attr(res_ls[[1]], "argu"))))
  else
    res_ls <- envelope_set_labs(res_ls, xlab=expression(italic(r)),
                                ylab=expression(italic(T(r))))
  attr(res_ls, "alternative") <- get_alternative(res_ls[[1]])
  attr(res_ls, "type") <- type
  attr(res_ls, "alpha") <- 1-coverage
  attr(res_ls, "M") <- attr(res_erl, "M")
  attr(res_ls, "M_alpha") <- attr(res_erl, "M_alpha")
  attr(res_ls, "p") <- attr(res_erl, "p")
  attr(res_ls, "nstep") <- 2
  class(res_ls) <- c("combined_global_envelope", class(res_ls))
  if(curve_set_is2d(curve_sets[[1]]))
    class(res_ls) <- c("combined_global_envelope2d", class(res_ls))
  res_ls
}

# Functionality for combined_central_region and combined_global_envelope_test (one-step procedure)
combined_CR_or_GET_1step <- function(curve_sets, CR_or_GET = c("CR", "GET"), coverage, ...) {
  curve_set <- combine_curve_sets(curve_sets, equalr=TRUE)
  switch(CR_or_GET,
         CR = {
           res <- individual_central_region(curve_set, coverage=coverage, ...)
         },
         GET = {
           res <- individual_global_envelope_test(curve_set, alpha=1-coverage, ...)
         })
  # Transform the envelope to a combined envelope
  nfuns <- length(curve_sets)
  nr <- curve_set_narg(curve_sets[[1]]) # all curve sets have the same
  idx <- lapply(1:nfuns, FUN = function(i) ((i-1)*nr+1):(i*nr))
  # Split the envelopes to the original groups
  res_ls <- split(res, f = rep(1:nfuns, each=nr))

  # Set unreasonable attributes of individuals sets of curves to NULL
  for(i in 1:nfuns)
    attr(res_ls[[i]], "method") <- paste0("1/", nfuns, "th of a combined global envelope test")
  anames <- c("p", "p_interval", "ties", "M", "M_alpha", "alpha")
  anames <- anames[anames %in% names(attributes(res_ls[[1]]))]
  for(name in anames) {
    for(i in 1:nfuns) attr(res_ls[[i]], name) <- NULL
  }

  mostattributes(res_ls) <- attributes(res)
  attr(res_ls, "row.names") <- NULL
  if(!is.null(names(curve_sets))) names(res_ls) <- names(curve_sets)
  switch(CR_or_GET,
         CR = {
           attr(res_ls, "method") <- "Combined global envelope"
         },
         GET = {
           attr(res_ls, "method") <- "Combined global test"
         })
  attr(res_ls, "nstep") <- 1
  class(res_ls) <- c("combined_global_envelope", "list")
  if(curve_set_is2d(curve_sets[[1]]))
    class(res_ls) <- c("combined_global_envelope2d", class(res_ls))
  res_ls
}

#' Print method for the class 'global_envelope'
#'
#' @param x A 'global_envelope' object.
#' @param ... Ignored.
#' @export
print.global_envelope <- function(x, ...) {
  printhelper_ge_base(x)
}

#' Print method for the class 'combined_global_envelope'
#'
#' @param x A 'combined_global_envelope' object
#' @param ... Ignored.
#' @export
print.combined_global_envelope <- function(x, ...) {
  printhelper_ge_combined(x)
}

#' Plot method for the class 'global_envelope'
#'
#' @param x An 'global_envelope' object
#' @param dotplot Logical. If TRUE, then instead of envelopes a dot plot is done.
#' Suitable for low dimensional test vectors.
#' Default: TRUE if the dimension is less than 10, FALSE otherwise.
#' @param sign.col The color for the observed curve when outside the global envelope
#' (significant regions). Default to "red". Setting the color to \code{NULL} corresponds
#' to no coloring.
#' @param labels A character vector of suitable length.
#' If \code{dotplot = TRUE}, then labels for the tests at x-axis.
#' @param digits The number of digits used for printing the p-value or p-interval
#' in the default main.
#' @param ... Ignored.
#'
#' @export
#' @seealso \code{\link{central_region}}, \code{\link{global_envelope_test}}
#' @examples
#' if(require("spatstat", quietly=TRUE)) {
#'   X <- unmark(spruces)
#'   \donttest{nsim <- 1999 # Number of simulations}
#'   \dontshow{nsim <- 19 # Number of simulations}
#'   env <- envelope(X, fun="Kest", nsim=nsim,
#'                   savefuns=TRUE, # save the functions
#'                   correction="translate", # edge correction for K
#'                   simulate=expression(runifpoint(ex=X))) # Simulate CSR
#'   res <- global_envelope_test(env, type="erl")
#'
#'   # Default plot
#'   plot(res)
#'   # Plots can be edited, e.g.
#'   # Remove legend
#'   plot(res) + ggplot2::theme(legend.position = "none")
#'   # Change its position
#'   plot(res) + ggplot2::theme(legend.position = "right")
#'   # Change the outside color
#'   plot(res, sign.col = "#5DC863FF")
#'   plot(res, sign.col = NULL)
#'   # Change default title and x- and y-labels
#'   plot(res) + ggplot2::labs(title="95% global envelope", x="x", y="f(x)")
#'
#'   # Prior to the plot, you can set your preferred ggplot theme by theme_set
#'   old <- ggplot2::theme_set(ggplot2::theme_bw())
#'   plot(res)
#'
#'   # Do other edits, e.g. turn off expansion with the default limits
#'   plot(res) + ggplot2::coord_cartesian(expand=FALSE)
#'
#'   # Go back to the old theme
#'   ggplot2::theme_set(old)
#'
#'   # If you are working with the R package spatstat and its envelope-function,
#'   # you can obtain global envelope plots in the style of spatstat using plot.fv:
#'   plot.fv(res)
#' }
#' @importFrom ggplot2 labs
plot.global_envelope <- function(x, dotplot = length(x$r)<10, sign.col = "red",
                                 labels = NULL, digits = 3, ...) {
  if(!is.null(x[['r']]) && !all(x[['r']][-1] - x[['r']][-length(x[['r']])] > 0))
    warning("r values non-increasing. Plot not valid.")
  if(missing(labels)) labels <- default_labels(x, labels)
  main <- env_main_default(x, digits=digits)
  d <- plotdefaultlabs(x)
  if(dotplot) {
    env_dotplot_ggplot(x, labels=labels, sign.col=sign.col) +
      labs(title=main, x=d$xlab, y=d$ylab)
  } else {
    env_ggplot(x, main=main, xlab=d$xlab, ylab=d$ylab, sign.col=sign.col)
  }
}

#' Plot method for the class 'combined_global_envelope'
#'
#' Plotting method for the class 'combined_global_envelope', i.e. combined envelopes for
#' 1d functions.
#'
#' @description This function provides plots for combined global envelopes.
#' @param x An 'combined_global_envelope' object
#' @inheritParams plot.global_envelope
#' @param labels A character vector of suitable length.
#' If \code{dotplot = TRUE} (for the level 2 test), then labels for the tests at x-axis.
#' Otherwise labels for the separate plots.
#' @param scales See \code{\link[ggplot2]{facet_wrap}}.
#' Use \code{scales = "free"} when the scales of the functions in the global envelope
#' vary. \code{scales = "fixed"} is a good choice, when you want the same y-axis for all components.
#' A sensible default based on r-values exists.
#' @param ncol The maximum number of columns for the figures.
#' Default 2 or 3, if the length of x equals 3.
#' (Relates to the number of curve_sets that have been combined.)
#' @param level 1 or 2. In the case of two-step combined tests (with several test functions),
#' two different plots are available:
#' 1 for plotting the combined global envelopes (default and most often wanted) or
#' 2 for plotting the second level test result.
#' @export
#' @seealso \code{\link{central_region}}
plot.combined_global_envelope <- function(x, labels, scales, sign.col = "red",
                                          ncol = 2 + 1*(length(x)==3),
                                          digits = 3, level = 1, ...) {
  if(!(level %in% c(1,2))) stop("Unreasonable value for level.")

  if(missing(scales)) {
    if(all(sapply(x, FUN=function(y) { all(range(y[['r']]) == range(x[[1]][['r']])) })))
      scales <- "fixed"
    else
      scales <- "free"
  }

  alt <- get_alternative(x[[1]])
  main <- env_main_default(x, digits=digits, alternative=alt)
  d <- plotdefaultlabs(x)

  if(level == 1) {
    if(missing(labels)) labels <- default_labels(x, labels)
    env_combined_ggplot(x, main=main, xlab=d$xlab, ylab=d$ylab,
                        labels=labels, scales=scales,
                        max_ncols_of_plots=ncol, sign.col=sign.col)
  }
  else {
    if(attr(x, "nstep") != 2)
      stop("level = 2 plot not available for one-step combined global envelopes.")
    if(missing(labels)) labels <- default_labels(attr(x, "level2_ge"), labels)
    env_dotplot_ggplot(attr(x, "level2_ge"), labels=labels)
  }
}

#' Central region / Global envelope
#'
#' Provides central regions or global envelopes or confidence bands
#'
#'
#' Given a \code{curve_set} (see \code{\link{create_curve_set}} for how to create such an object)
#' or an \code{envelope} object of \pkg{spatstat} or \code{fdata} object of \pkg{fda.usc},
#' the function \code{central_region} construcst a central region, i.e. a global envelope,
#' from the given set of functions (or vectors).
#'
#' Generally an envelope is a band bounded by the vectors (or functions)
#' \eqn{T_{low}}{T_lo} and \eqn{T_{hi}}{T_hi}.
#' A \eqn{100(1-\alpha)}{100(1-alpha)}\% or 100*coverage\% global envelope is a set
#' \eqn{(T_{low}, T_{hi})}{(T_lo, T_hi)} of envelope vectors
#' such that the probability that \eqn{T_i}{T_i} falls outside this envelope
#' in any of the d points of the vector \eqn{T_i}{T_i} is less or equal to \eqn{\alpha}{alpha}.
#' The global envelopes can be constructed based on different measures
#' that order the functions from the most extreme one to the least extreme one.
#' We use such orderings of the functions for which we are able to construct global envelopes
#' with intrinsic graphical interpretation.
#'
#' The type of the global envelope can be chosen with the argument \code{type} and
#' the options are given in the following.
#' Further information about the measures, on which the global envelopes are based,
#' can be found in Myllymäki and Mrkvička (2020, Section 2.).
#' \itemize{
#'  \item \code{'rank'}: The global rank envelope
#' proposed by Myllymäki et al. (2017) based on the extreme rank defined as the minimum of pointwise
#' ranks.
#'  \item \code{'erl'}: The global rank envelope based on the extreme rank
#'  length (Myllymäki et al.,2017, Mrkvička et al., 2018).
#' This envelope is constructed as the convex hull of the functions which have extreme rank
#' length measure that is larger or equal to the critical \eqn{\alpha}{alpha} level of the
#' extreme rank length measure.
#'  \item \code{'cont'}: The global rank envelope based on the continuous rank
#'  (Hahn, 2015; Mrkvička et al., 2019) based on minimum of continuous pointwise ranks.
#'  It is contructed as the convex hull in a similar way as the \code{'erl'} envelope.
#'  \item \code{'area'}: The global rank envelope based on the area rank (Mrkvička et al., 2019)
#'  which is based on area between continuous pointwise ranks and minimum pointwise ranks
#'  for those argument (r) values for which pointwise ranks achieve the minimum
#'  (it is a combination of erl and cont).
#'  It is contructed as the convex hull in a similar way as the \code{'erl'} and \code{'area'} envelopes.
#'  \item \code{'qdir'}: The directional quantile envelope based on
#'  the directional quantile maximum absolute deviation (MAD) test (Myllymäki et al., 2017, 2015),
#' which takes into account the unequal variances of the test function T(r) for
#' different distances r and is also protected against asymmetry of distribution of T(r).
#'  \item \code{'st'}: The studentised envelope based on the studentised MAD
#'  measure (Myllymäki et al., 2017, 2015),
#'  which takes into account the unequal variances of the test function T(r) for different distances r.
#'  \item \code{'unscaled'}: The unscaled envelope (Ripley, 1981),
#' which leads to envelopes with constant width. It corresponds to the classical
#' maximum deviation test without scaling. This test suffers from unequal variance
#' of T(r) over the distances r and from the asymmetry of distribution of T(r).
#' We recommend to use the other alternatives instead. This unscaled global envelope is
#' provided for reference.
#' }
#'
#' The values of the chosen measure M are determined for each curve in the \code{curve_set}, and
#' based on the chosen measure, the central region, i.e. the global envelope, is constructed
#' for the given curves.
#'
#' If a list of (suitable) objects are provided in the argument \code{curve_sets},
#' then by default (\code{nstep = 2}) the two-step combining procedure is used to
#' construct the combined global envelope as described in Myllymäki and Mrkvička (2020, Section 2.2.).
#' If \code{nstep = 1} and the lengths of the multivariate vectors in each component
#' of the list are equal, then the one-step combining procedure is used where the
#' functions are concatenated together into a one long vector (see again Myllymäki and Mrkvička, 2020, Section 2.2.).
#'
#' @references
#' Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2020) A one-way ANOVA test for functional data with graphical interpretation. Kybernetika 56 (3), 432-458. doi: 10.14736/kyb-2020-3-0432
#'
#' Mrkvička, T., Myllymäki, M., Kuronen, M. and Narisetty, N. N. (2020) New methods for multiple testing in permutation inference for the general linear model. arXiv:1906.09004 [stat.ME]
#'
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2015). Deviation test construction and power comparison for marked spatial point patterns. Spatial Statistics 11: 19-34. doi: 10.1016/j.spasta.2014.11.004
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' Myllymäki, M. and Mrkvička, T. (2020). GET: Global envelopes in R. arXiv:1911.06583 [stat.ME]
#'
#' Ripley, B.D. (1981). Spatial statistics. Wiley, New Jersey.
#'
#' @inheritParams forder
#' @param type The type of the global envelope with current options for 'rank', 'erl', 'cont', 'area',
#' 'qdir', 'st' and 'unscaled'. See details.
#' @param coverage A number between 0 and 1. The 100*coverage\% central region will be calculated.
#' @param central Either "mean" or "median". If the curve sets do not contain the component
#' \code{theo} for the theoretical central function, then the central function (used for plotting only)
#' is calculated either as the mean or median of functions provided in the curve sets.
#' @param nstep 1 or 2 for how to contruct a combined global envelope if list of curve sets
#' is provided. 2 (default) for a two-step combining procedure, 1 for one-step.
#' @param ... Ignored.
#' @return Either an object of class \code{global_envelope} and or an \code{combined_global_envelope} object.
#' The former class is obtained when a set of curves is provided, while the latter in the case
#' that \code{curve_sets} is a list of objects. The print and plot function are defined for the
#' returned objects (see examples).
#'
#' The \code{global_envelope} object is essentially a data frame containing columns
#' \itemize{
#' \item r = the vector of values of the argument r at which the test was made
#' \item lo = the lower envelope based on the simulated functions
#' \item hi = the upper envelope based on the simulated functions
#' \item central = If the \code{curve_set} (or \code{envelope} object) contains a theoretical curve,
#'       then this function is used as the central curve and returned in this component.
#'       Otherwise, the central curve is the mean or median (according to the argument \code{central})
#'       of the test functions T_i(r), i=2, ..., s+1. Used for visualization only.
#' }
#' and potentially additionally
#' \itemize{
#' \item obs = the data function, if there is only one data function in the given \code{curve_sets}.
#' Otherwise not existing.
#' }
#' Most often \code{central_region} is directly applied to functional data where all curves are observed.
#' Additionally, the returned object has some attributes, where
#' \itemize{
#'   \item M = A vector of the values of the chosen measure for all the function.
#'   If there is only one observed function, then M[1] gives the value of the measure for this.
#'   \item M_alpha = The critical value of M corresponding to the 100(1-alpha)\% global envelope
#'   (see Myllymäki and Mrkvička, 2020, Definition 1.1. IGI).
#' }
#' Further the object has some attributes for printing and plotting purposes, where
#' \code{alternative}, \code{type}, \code{ties}, \code{alpha} correspond to those in the function call
#' and \code{method} gives a name for the method.
#' Attributes of an object \code{res} can be obtained using the function
#' \code{\link[base]{attr}}, e.g. \code{attr(res, "M")} for the values of the ordering measure.
#'
#' If the given set of curves had the class \code{envelope} of \pkg{spatstat}, then the returned
#' \code{global_envelope} object has also the class \code{fv} of spatstat, whereby one can utilize
#' also the plotting functions of \pkg{spatstat}, see example in \code{\link{plot.global_envelope}}.
#' However, the \code{envelope} objects are most often used with \code{\link{global_envelope_test}}
#' and not with \code{central_region}.
#' For an \code{fv} object, also some further attributes exists as required by \code{fv} of \pkg{spatstat}.
#'
#' The \code{combined_global_envelope} is a list of \code{global_envelope} objects, where
#' the components correspond to the components of \code{curve_sets}.
#' The \code{combined_global_envelope} object constructed with \code{nstep = 2} contains,
#' in addition to some conventional ones (\code{method}, \code{alternative}, \code{type}, \code{alpha},
#' \code{M}, \code{M_alpha}, see above), the second level envelope information as the attributes
#' \itemize{
#' \item level2_ge = The second level envelope on which the envelope construction is based
#' \item level2_curve_set = The second level \code{curve_set} from which \code{level2_ge} is constructed
#' }
#'
#' In the case that the given curve sets are two-dimensional, i.e., their arguments values are two-dimensional,
#' then the returned objects have in addition to the class \code{global_envelope} or \code{combined_global_envelope},
#' the class \code{global_envelope2d} or \code{combined_global_envelope2d}, respectively. This class is assigned
#' for plotting purposes: For the 2d envelopes, also the default plots are 2d.
#' Otherwise the 1d and 2d objects are similar.
#' @export
#' @seealso \code{\link{forder}}, \code{\link{global_envelope_test}}
#' @aliases global_envelope
#' @examples
#' ## A central region of a set of functions
#' #----------------------------------------
#' if(requireNamespace("fda", quietly = TRUE)) {
#'   curve_set <- create_curve_set(list(r=as.numeric(row.names(fda::growth$hgtf)),
#'                                      obs=fda::growth$hgtf))
#'   plot(curve_set, ylab="height")
#'   cr <- central_region(curve_set, coverage=0.50, type="erl")
#'   plot(cr)
#' }
#'
#' ## Confidence bands for linear or polynomial regression
#' #------------------------------------------------------
#' # Simulate regression data according to the cubic model
#' # f(x) = 0.8x - 1.8x^2 + 1.05x^3 for x in [0,1]
#' par <- c(0,0.8,-1.8,1.05) # Parameters of the true polynomial model
#' res <- 100 # Resolution
#' x <- seq(0, 1, by=1/res); x2=x^2; x3=x^3;
#' f <- par[1] + par[2]*x + par[3]*x^2 + par[4]*x^3 # The true function
#' d <- f + rnorm(length(x), 0, 0.04) # Data
#' # Plot the true function and data
#' plot(f, type="l", ylim=range(d))
#' points(d)
#'
#' # Estimate polynomial regression model
#' reg <- lm(d ~ x + x2 + x3)
#' ftheta <- reg$fitted.values
#' resid0 <- reg$residuals
#' s0 <- sd(resid0)
#'
#' # Bootstrap regression
#' \donttest{B <- 2000 # Number of bootstrap samples}
#' \dontshow{B <- 20 # Number of bootstrap samples}
#' ftheta1 <- array(0, c(B,length(x)))
#' s1 <- array(0,B)
#' for(i in 1:B) {
#'   u <- sample(resid0, size=length(resid0), replace=TRUE)
#'   reg1 <- lm((ftheta+u) ~ x + x2 + x3)
#'   ftheta1[i,] <- reg1$fitted.values
#'   s1[i] <- sd(reg1$residuals)
#' }
#'
#' # Centering and scaling
#' meanftheta <- apply(ftheta1, 2, mean)
#' m <- array(0, c(B,length(x)))
#' for(i in 1:B) { m[i,] <- (ftheta1[i,]-meanftheta)/s1[i] }
#'
#' # Central region computation
#' boot.cset <- create_curve_set(list(r=1:length(x), obs=ftheta+s0*t(m)))
#' cr <- central_region(boot.cset, coverage=0.95, type="erl")
#'
#' # Plotting the result
#' plot(cr) + ggplot2::labs(x = expression(italic(x)), y = expression(italic(f(x)))) +
#'   ggplot2::geom_point(data = data.frame(id = 1:length(d), points = d),
#'                       ggplot2::aes(x = id, y = points)) + # data points
#'   ggplot2::geom_line(data = data.frame(id = 1:length(d), points = f),
#'                      ggplot2::aes(x = id, y = points)) # true function
central_region <- function(curve_sets, type = "erl", coverage = 0.50,
                           alternative = c("two.sided", "less", "greater"),
                           probs = c((1-coverage)/2, 1-(1-coverage)/2),
                           quantile.type = 7,
                           central = "median", nstep = 2, ...) {
  if(length(class(curve_sets)) == 1 && class(curve_sets) == "list") {
    if(length(curve_sets) > 1) { # Combined test
      if(!(nstep %in% c(1,2))) stop("Invalid number of steps (nstep) for combining. Should be 1 or 2.")
      if(nstep == 2) # Two-step combining procedure
        return(combined_CR_or_GET(curve_sets, CR_or_GET="CR", type=type, coverage=coverage,
                                  alternative=alternative,
                                  probs=probs, quantile.type=quantile.type,
                                  central=central, ...))
      else # One-step combining procedure
        return(combined_CR_or_GET_1step(curve_sets, CR_or_GET="CR", type=type, coverage=coverage,
                                        alternative=alternative,
                                        probs=probs, quantile.type=quantile.type,
                                        central=central, ...))
    }
    else if(length(curve_sets) == 1)
      curve_sets <- curve_sets[[1]]
    else
      stop("The given list of curve_sets is empty.")
  }
  # Individual test
  return(individual_central_region(curve_sets, type=type, coverage=coverage,
                                   alternative=alternative,
                                   probs=probs, quantile.type=quantile.type,
                                   central=central, ...))
}


#' Global envelope test
#'
#' Global envelope test, global envelopes and p-values
#'
#'
#' Given a \code{curve_set} (see \code{\link{create_curve_set}} for how to create such an object)
#' or an \code{envelope} object of \pkg{spatstat},
#' which contains both the data curve (or function or vector) \eqn{T_1(r)}{T_1(r)}
#' (in the component \code{obs}) and
#' the simulated curves \eqn{T_2(r),\dots,T_{s+1}(r)}{T_2(r),...,T_(s+1)(r)}
#' (in the component \code{sim_m}),
#' the function \code{global_envelope_test} performs a global envelope test.
#' The functionality of the function is rather similar to the function
#' \code{\link{central_region}}, but in addition to ordering the functions from
#' the most extreme one to the least extreme one using different measures
#' and providing the global envelopes with intrinsic
#' graphical interpretation, p-values are calculated for the test.
#' Thus, while \code{\link{central_region}} can be used to construct global
#' envelopes in a general setting, the function \code{\link{global_envelope_test}}
#' is devoted to testing as its name suggests.
#'
#' The function \code{global_envelope_test} is the main function for global envelope tests
#' (for simple hypotheses).
#' Different \code{type} of global envelope tests can be performed.
#' We use such ordering of the functions for which we are able to construct global
#' envelopes with intrinsic graphical interpretation.
#' \itemize{
#'   \item \code{'rank'}: the completely non-parametric rank envelope test (Myllymäki et al., 2017)
#'   based on minimum of pointwise ranks
#'   \item \code{'erl'}: the completely non-parametric rank envelope test based on extreme rank lengths
#'   (Myllymäki et al., 2017; Mrkvička et al., 2018) based on number of minimal pointwise ranks
#'  \item \code{'cont'}: the completely non-parametric rank envelope test based on continuous rank
#'  (Hahn, 2015; Mrkvička et al., 2019) based on minimum of continuous pointwise ranks
#'  \item \code{'area'}: the completely non-parametric rank envelope test based on area rank
#'  (Mrkvička et al., 2019) based on area between continuous pointwise ranks and minimum
#'  pointwise ranks for those argument (r) values for which pointwise ranks achieve the minimum
#'  (it is a combination of erl and cont)
#'   \item "qdir", the directional quantile envelope test, protected against unequal variance and
#'   asymmetry of T(r) for different distances r (Myllymäki et al., 2015, 2017)
#'   \item "st", the studentised envelope test, protected against unequal variance of T(r) for
#'   different distances r (Myllymäki et al., 2015, 2017)
#'   \item "unscaled", the unscaled envelope (providing a baseline) that has a contant width and
#'   that corresponds to the classical maximum deviation test (Ripley, 1981).
#' }
#' The first four types are global rank envelopes.
#' The \code{'rank'} envelope test is a completely non-parametric test,
#' which provides the 100(1-alpha)% global envelope for the chosen test function
#' T(r) on the chosen interval of distances and associated p-values.
#' The other three are modifications of \code{'rank'} to treat the ties in
#' the extreme rank ordering on which the \code{'rank'} test is based on.
#' The last three envelopes are global scaled maximum absolute difference (MAD)
#' envelope tests. The unscaled envelope test leads to envelopes with constant width over the
#' distances r. Thus, it suffers from unequal variance of T(r) over the distances r and
#' from the asymmetry of distribution of T(r). We recommend to use the other global
#' envelope tests available. The unscaled envelope is provided as a reference.
#'
#' See Myllymäki and Mrkvička (2020, Section 2.), i.e. \code{vignette("GET")}, for more detailed description of the measures and
#' the corresponding envelopes.
#'
#' @section Procedure:
#' 1) First the curves are ranked from the most extreme one to the least extreme one
#' by a measure that is specified by the argument \code{type}. The options are
#' \itemize{
#' \item 'rank': extreme ranks (Myllymäki et al., 2017)
#' \item 'erl': extreme rank lengths (Myllymäki et al., 2017; Mrkvička et al., 2018)
#' \item 'cont': continuous ranks (Hahn, 2015; Mrkvička et al., 2019)
#' \item 'area': area ranks (Mrkvička et al., 2019)
#' \item 'qdir': the directional quantile maximum absolute deviation (MAD) measure (Myllymäki et al., 2015, 2017)
#' \item 'st': the studentized MAD measure (Myllymäki et al., 2015, 2017)
#' \item 'unscaled': the unscaled MAD measure (Ripley, 1981)
#' }
#'
#' 2) Based on the measures used to rank the functions, the 100(1-alpha)\% global envelope is provided.
#' It corresponds to the 100*coverage\% central region.
#'
#' 3) P-values:
#' In the case \code{type="rank"}, based on the extreme ranks \eqn{k_i, i=1, ..., s+1}{k_i, i=1, ..., s+1},
#' the p-interval is calculated. Because the extreme ranks contain ties, there is not just
#' one p-value. The p-interval is given by the most liberal and the most conservative p-value
#' estimate. Also a single p-value is calculated.
#' By default this single p-value is the extreme rank length p-value ("erl") as specified by the argument \code{ties}.
#' If the case of other measures, a (single) p-value based on the given ordering
#' of the functions is calculated and returned in the attribute \code{p}.
#'
#' @section Number of simulations:
#' For the global \code{"rank"} envelope test, Myllymäki et al. (2017) recommended to use
#' at least 2500 simulations for testing at the significance level alpha = 0.05 for single
#' function tests, based on experiments with summary functions for point processes evaluated
#' approximately at 500 argument values.
#' In this case, the width of the p-interval associated with the extreme rank measure tended
#' to be smaller than 0.01.
#' The tests \code{'erl'}, \code{'cont'} and \code{'area'}, similarly as
#' the MAD deviation/envelope tests \code{'qdir'}, \code{'st'} and \code{'unscaled'},
#' allow in principle a lower number of simulations to be used than the test based on
#' extreme ranks (\code{'rank'}), because no ties occur for these measures.
#' If affordable, we recommend in any case some thousands of simulations for all the measures
#' to achieve a good power and repeatability of the test.
#' If the dimension of the test functions is higher, also the number of simulations should
#' preferably be higher.
#'
#' @section Tests based on several functions:
#' If a list of (suitable) objects are provided in the argument \code{curve_sets},
#' then by default (\code{nstep = 2}) the two-step combining procedure is used to
#' perform the combined global test as described in Myllymäki and Mrkvička (2020).
#' If \code{nstep = 1} and the lengths of the multivariate vectors in each component
#' of the list are equal, then the one-step combining procedure is used where the
#' functions are concatenated together into a one long vector.
#'
#' @references
#' Mrkvička, T., Myllymäki, M. and Hahn, U. (2017). Multiple Monte Carlo testing, with applications in spatial point processes. Statistics & Computing 27 (5): 1239-1255. doi: 10.1007/s11222-016-9683-9
#'
#' Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2020) A one-way ANOVA test for functional data with graphical interpretation. Kybernetika 56 (3), 432-458. doi: 10.14736/kyb-2020-3-0432
#'
#' Mrkvička, T., Myllymäki, M., Kuronen, M. and Narisetty, N. N. (2020) New methods for multiple testing in permutation inference for the general linear model. arXiv:1906.09004 [stat.ME]
#'
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2015). Deviation test construction and power comparison for marked spatial point patterns. Spatial Statistics 11: 19-34. doi: 10.1016/j.spasta.2014.11.004
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' Myllymäki, M. and Mrkvička, T. (2020). GET: Global envelopes in R. arXiv:1911.06583 [stat.ME]
#'
#' Ripley, B.D. (1981). Spatial statistics. Wiley, New Jersey.
#'
#' @inheritParams central_region
#' @param curve_sets A \code{curve_set} (see \code{\link{create_curve_set}})
#' or an \code{envelope} object of \pkg{spatstat} containing a data function and simulated functions.
#' If an envelope object is given, it must contain the summary
#' functions from the simulated patterns which can be achieved by setting
#' \code{savefuns = TRUE} when calling the \code{envelope} function.
#' Alternatively, a list of \code{curve_set} or \code{envelope} objects can be given.
#' @param alpha The significance level. The 100(1-alpha)\% global envelope will be calculated.
#' @param ties The method to obtain a unique p-value when  \code{type = 'rank'}.
#' Possible values are 'midrank', 'random', 'conservative', 'liberal' and 'erl'.
#' For 'conservative' the resulting p-value will be the highest possible.
#' For 'liberal' the p-value will be the lowest possible.
#' For 'random' the rank of the obs within the tied values is uniformly sampled so that the resulting
#' p-value is at most the conservative option and at least the liberal option.
#' For 'midrank' the mid-rank within the tied values is taken.
#' For 'erl' the extreme rank length p-value is calculated.
#' The default is 'erl'.
#' @param ... Additional parameters to be passed to \code{\link{central_region}}.
#' @return Either an object of class "global_envelope" or "combined_global_envelope",
#' similarly as the objects returned by \code{\link{central_region}}.
#'
#' The \code{global_envelope} is essentially a data frame containing columns
#' \itemize{
#' \item the values of the argument r at which the test was made, copied from the argument \code{curve_sets} with the corresponding names
#' \item obs = values of the data function, copied from the argument \code{curve_sets}
#' (unlike for central regions, \code{obs} always exists for a global envelope test)
#' \item lo = the lower envelope
#' \item hi = the upper envelope
#' \item central = a central curve as specified in the argument \code{central}.
#' }
#' Moreover, the returned object has the same attributes as the \code{global_envelope} object returned by
#' \code{\link{central_region}} and in addition
#' \itemize{
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#' }
#' and in the case that \code{type = 'rank'} also
#' \itemize{
#'   \item p_interval = The p-value interval \eqn{[p_{liberal}, p_{conservative}]}{[p_liberal, p_conservative]}.
#'   \item ties = As the argument \code{ties}.
#' }
#'
#' The \code{combined_global_envelope} is a list of \code{global_envelope} objects
#' containing the above mentioned columns and which all together form the global envelope.
#' It has the same attributes as described in \code{\link{central_region}}, and in addition also
#' the p-value \code{p}.
#' The 2d classes are attached as described in \code{\link{central_region}}.
#' @export
#' @seealso \code{\link{plot.global_envelope}}, \code{\link{central_region}},
#' \code{\link{GET.composite}}
#' @examples
#' # Goodness-of-fit testing for simple hypothesis
#' if(require("spatstat", quietly=TRUE)) {
#'   # Testing complete spatial randomness (CSR)
#'   #==========================================
#'   X <- unmark(spruces)
#'
#'   \donttest{nsim <- 1999 # Number of simulations}
#'   \dontshow{nsim <- 19 # Number of simulations}
#'
#'   # Illustration of general workflow for simple hypotheses
#'   #=======================================================
#'   # First illustrate the general workflow for the test by this example
#'   # of CSR test for a point pattern X using the empirical L-function.
#'   # Define the argument values at which the functions are evaluated
#'   obs.L <- Lest(X, correction = "translate")
#'   r <- obs.L[['r']]
#'   # The test function for the data
#'   obs <- obs.L[['trans']] - r
#'   # Prepare simulations and calculate test functions for them at same r as 'obs'
#'   sim <- matrix(nrow = length(r), ncol = nsim)
#'   for(i in 1:nsim) {
#'     sim.X <- runifpoint(ex = X) # simulation under CSR
#'     sim[, i] <- Lest(sim.X, correction = "translate", r = r)[['trans']] - r
#'   }
#'   # Create a curve_set containing argument values, observed and simulated functions
#'   cset <- create_curve_set(list(r = r, obs = obs, sim_m = sim))
#'   # Perform the test
#'   res <- global_envelope_test(cset, type="erl")
#'   plot(res, ylab = expression(italic(hat(L)(r)-r)))
#'
#'   # Simple hypothesis for a point pattern utilizing the spatstat package
#'   #=====================================================================
#'   # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#'   env <- envelope(X, fun="Lest", nsim=nsim,
#'                   savefuns=TRUE, # save the functions
#'                   correction="translate", # edge correction for L
#'                   transform = expression(.-r), # centering
#'                   simulate=expression(runifpoint(ex=X))) # Simulate CSR
#'   # The rank envelope test (ERL)
#'   res <- global_envelope_test(env, type="erl")
#'   # Plot the result
#'   plot(res)
#'
#'   ## Advanced use:
#'   # Choose the interval of distances [r_min, r_max] (at the same time create a curve_set from 'env')
#'   cset <- crop_curves(env, r_min=1, r_max=7)
#'   # Do the rank envelope test (erl)
#'   res <- global_envelope_test(cset, type="erl")
#'   plot(res, ylab=expression(italic(L(r)-r)))
#'
#'   \donttest{
#'   # Random labeling test
#'   #=====================
#'   mpp <- spruces
#'   # 1) Perform simulations under the random labelling hypothesis and calculate
#'   # the test function T(r) for the data pattern (mpp) and each simulation.
#'   # The command below specifies that the test function is T(r) = \hat{L}_mm(r),
#'   # which is an estimator of the mark-weighted L function, L_mm(r),
#'   # with translational edge correction.
#'   nsim <- 1999 # Number of simulations
#'   env <- envelope(mpp, fun=Kmark, nsim = nsim, f=function(m1, m2) { m1*m2 },
#'                   correction="translate", returnL=TRUE,
#'                   simulate=expression(rlabel(mpp, permute=TRUE)), # Permute the marks
#'                   savefuns=TRUE) # Save the functions
#'   # 2)
#'   # Crop curves to desired r-interval
#'   curve_set <- crop_curves(env, r_min=1.5, r_max=9.5)
#'   # Center the functions, i.e. take \hat{L}_mm(r)-T_0(r).
#'   # Below T_0(r) = \hat{L}(r) is the mean of simulated functions.
#'   # (This is only for visualization, does not affect the test result.)
#'   curve_set <- residual(curve_set)
#'   # 3) Do the rank envelope test
#'   res <- global_envelope_test(curve_set)
#'   # 4) Plot the test result
#'   plot(res, ylab=expression(italic(L[mm](r)-L(r))))
#'
#'   # A combined global envelope test
#'   #================================
#'   # As an example test CSR of the saplings point pattern by means of
#'   # L, F, G and J functions.
#'   data(saplings)
#'   X <- as.ppp(saplings, W=square(75))
#'
#'   \donttest{nsim <- 499 # Number of simulations}
#'   \dontshow{nsim <- 19 # Number of simulations}
#'   # Specify distances for different test functions
#'   n <- 500 # the number of r-values
#'   rmin <- 0; rmax <- 20; rstep <- (rmax-rmin)/n
#'   rminJ <- 0; rmaxJ <- 8; rstepJ <- (rmaxJ-rminJ)/n
#'   r <- seq(0, rmax, by=rstep)    # r-distances for Lest
#'   rJ <- seq(0, rmaxJ, by=rstepJ) # r-distances for Fest, Gest, Jest
#'
#'   # Perform simulations of CSR and calculate the L-functions
#'   env_L <- envelope(X, nsim=nsim,
#'    simulate=expression(runifpoint(ex=X)),
#'    fun="Lest", correction="translate",
#'    transform=expression(.-r), # Take the L(r)-r function instead of L(r)
#'    r=r,                       # Specify the distance vector
#'    savefuns=TRUE,             # Save the estimated functions
#'    savepatterns=TRUE)         # Save the simulated patterns
#'   # Take the simulations from the returned object
#'   simulations <- attr(env_L, "simpatterns")
#'   # Then calculate the other test functions F, G, J for each simulated pattern
#'   env_F <- envelope(X, nsim=nsim,
#'                     simulate=simulations,
#'                     fun="Fest", correction="Kaplan", r=rJ,
#'                     savefuns=TRUE)
#'   env_G <- envelope(X, nsim=nsim,
#'                     simulate=simulations,
#'                     fun="Gest", correction="km", r=rJ,
#'                     savefuns=TRUE)
#'   env_J <- envelope(X, nsim=nsim,
#'                     simulate=simulations,
#'                     fun="Jest", correction="none", r=rJ,
#'                     savefuns=TRUE)
#'
#'   # Crop the curves to the desired r-interval I
#'   curve_set_L <- crop_curves(env_L, r_min=rmin, r_max=rmax)
#'   curve_set_F <- crop_curves(env_F, r_min=rminJ, r_max=rmaxJ)
#'   curve_set_G <- crop_curves(env_G, r_min=rminJ, r_max=rmaxJ)
#'   curve_set_J <- crop_curves(env_J, r_min=rminJ, r_max=rmaxJ)
#'
#'   res <- global_envelope_test(curve_sets=list(L = curve_set_L, F = curve_set_F,
#'                                               G = curve_set_G, J = curve_set_J))
#'   plot(res)
#'   plot(res, labels=c("L(r)-r", "F(r)", "G(r)", "J(r)"))
#'   }
#' }
#'
#' # A test based on a low dimensional random vector
#' #================================================
#' # Let us generate some example data.
#' X <- matrix(c(-1.6,1.6),1,2) # data pattern X=(X_1,X_2)
#' if(requireNamespace("mvtnorm", quietly=TRUE)) {
#'   Y <- mvtnorm::rmvnorm(200,c(0,0),matrix(c(1,0.5,0.5,1),2,2)) # simulations
#'   plot(Y, xlim=c(min(X[,1],Y[,1]), max(X[,1],Y[,1])), ylim=c(min(X[,2],Y[,2]), max(X[,2],Y[,2])))
#'   points(X, col=2)
#'
#'   # Test the null hypothesis is that X is from the distribution of Y's (or if it is an outlier).
#'
#'   # Case 1. The test vector is (X_1, X_2)
#'   cset1 <- create_curve_set(list(r=1:2, obs=as.vector(X), sim_m=t(Y)))
#'   res1 <- global_envelope_test(cset1)
#'   plot(res1)
#'
#'   # Case 2. The test vector is (X_1, X_2, (X_1-mean(Y_1))*(X_2-mean(Y_2))).
#'   t3 <- function(x, y) { (x[,1]-mean(y[,1]))*(x[,2]-mean(y[,2])) }
#'   cset2 <- create_curve_set(list(r=1:3, obs=c(X[,1],X[,2],t3(X,Y)), sim_m=rbind(t(Y), t3(Y,Y))))
#'   res2 <- global_envelope_test(cset2)
#'   plot(res2)
#' }
global_envelope_test <- function(curve_sets, type = "erl", alpha = 0.05,
                          alternative = c("two.sided", "less", "greater"),
                          ties = "erl", probs = c(0.025, 0.975), quantile.type=7,
                          central = "mean", nstep = 2, ...) {
  if(length(class(curve_sets)) == 1 && class(curve_sets) == "list") {
    if(length(curve_sets) > 1) { # Combined test
      if(!(nstep %in% c(1,2))) stop("Invalid number of steps (nstep) for combining. Should be 1 or 2.")
      if(nstep == 2) # Two-step combining procedure
        return(combined_CR_or_GET(curve_sets, CR_or_GET="GET", type=type, coverage=1-alpha,
                                  alternative=alternative,
                                  probs=probs, quantile.type=quantile.type,
                                  central=central, ...))
      else # One-step combining procedure
        return(combined_CR_or_GET_1step(curve_sets, CR_or_GET="GET", type=type, coverage=1-alpha,
                                        alternative=alternative,
                                        probs=probs, quantile.type=quantile.type,
                                        central=central, ...))
    }
    else if(length(curve_sets) == 1)
      curve_sets <- curve_sets[[1]]
    else
      stop("The given list of curve_sets is empty.")
  }
  return(individual_global_envelope_test(curve_sets, type=type, alpha=alpha,
                                         alternative=alternative, ties=ties,
                                         probs=probs, quantile.type=quantile.type,
                                         central=central, ...))
}

#' The rank envelope test
#'
#' The rank envelope test, p-values and global envelopes.
#' The test corresponds to the global envelope test that can be carriet out by
#' \code{\link{global_envelope_test}} by specifying the \code{type} for which the options
#' \code{"rank"}, \code{"erl"}, \code{"cont"} and \code{"area"} are available. The last
#' three are modifications of the first one to treat the ties in the extreme rank ordering
#' used in \code{"rank"}. This function is kept for historical reasons.
#'
#' The \code{"rank"} envelope test is a completely non-parametric test, which provides
#' the 100(1-alpha)\% global envelope for the chosen test function T(r) on
#' the chosen interval of distances and associated p-values.
#' The other three types are solutions to break the ties in the extreme ranks
#' on which the \code{"rank"} envelope test is based on.
#'
#' Note: The method to break ties for the global \code{type = "rank"} envelope
#' (Myllymäki et al., 2017) can be done by the argument \code{ties} with default
#' to \code{ties = "erl"} corresponding to the extreme rank length breaking of ties.
#' In this case the global envelope corresponds to the extreme rank measure.
#' If instead choosing \code{type} to be \code{"erl"}, \code{"cont"} or \code{"area"},
#' then the global envelope corresponds to these measures.
#'
#' @section Number of simulations:
#' The global \code{"erl"}, \code{"cont"}, \code{"area"} envelope tests allow 
#' in principle a lower number of simulations to be used than the global \code{"rank"} test
#' based on extreme ranks.
#' However, if feasible, we recommend some thousands of simulations in any case to achieve
#' a good power and repeatability of the test.
#' For the global \code{"rank"} envelope test, Myllymäki et al. (2017) recommended to use
#' at least 2500 simulations for testing at the significance level alpha = 0.05 for single
#' function tests, experimented with summary functions for point processes.
#'
#' @references
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' Mrkvička, T., Myllymäki, M. and Hahn, U. (2017). Multiple Monte Carlo testing, with applications in spatial point processes. Statistics & Computing 27 (5): 1239-1255. doi: 10.1007/s11222-016-9683-9
#'
#' Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2020) A one-way ANOVA test for functional data with graphical interpretation. Kybernetika 56 (3), 432-458. doi: 10.14736/kyb-2020-3-0432
#'
#' @param curve_set A curve_set (see \code{\link{create_curve_set}}) or an \code{envelope}
#'  object of \pkg{spatstat}. If an envelope object is given, it must contain the summary
#'  functions from the simulated patterns which can be achieved by setting
#'  savefuns = TRUE when calling the function of \pkg{spatstat}.
#' @param type The type of the global envelope with current options for "rank", "erl", "cont" and "area".
#' If "rank", the global rank envelope accompanied by the p-interval is given (Myllymäki et al., 2017).
#' If "erl", the global rank envelope based on extreme rank lengths accompanied by the extreme rank
#' length p-value is given (Myllymäki et al., 2017, Mrkvička et al., 2018). See details and additional
#' sections thereafter.
#' @param ... Additional parameters to be passed to \code{\link{global_envelope_test}}.
#' @return An object of class \code{global_envelope} of \code{combined_global_envelope}
#' which can be printed and plotted directly. See \code{\link{global_envelope_test}} for more details.
#' @export
#' @seealso \code{\link{global_envelope_test}}
#' @examples
#' # See ?global_envelope_test for more examples
#'
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' if(require("spatstat", quietly=TRUE)) {
#'   X <- unmark(spruces)
#'   \donttest{nsim <- 2499 # Number of simulations}
#'   \dontshow{nsim <- 19 # Number of simulations for testing}
#'   # Generate nsim simulations under CSR, calculate centred L-function for the data and simulations
#'   env <- envelope(X, fun="Lest", nsim=nsim, savefuns=TRUE,
#'                   correction="translate", transform = expression(.-r),
#'                   simulate=expression(runifpoint(ex=X)))
#'   # The rank envelope test
#'   res <- rank_envelope(env)
#'   # Plot the result.
#'   plot(res)
#'
#'   ## Advanced use:
#'   # Choose the interval of distances [r_min, r_max] (at the same time create a curve_set from 'env')
#'   curve_set <- crop_curves(env, r_min=1, r_max=7)
#'   # Do the rank envelope test
#'   res <- rank_envelope(curve_set); plot(res)
#' }
rank_envelope <- function(curve_set, type = "rank", ...) {
  if(!(type %in% c("rank", "erl", "cont", "area"))) stop("No such type for the global rank envelope.")
  global_envelope_test(curve_set, type=type, ...)
}

#' Global scaled maximum absolute difference (MAD) envelope tests
#'
#' Performs the global scaled MAD envelope tests, either directional quantile or studentised,
#' or the unscaled MAD envelope test. These tests correspond to calling the
#' function \code{\link{global_envelope_test}} with \code{type="qdir"}, \code{type = "st"} and
#' \code{type="unscaled"}, respectively. The functions \code{qdir_envelope}, \code{st_envelope} and
#' \code{unscaled_envelope} have been kept for historical reasons;
#' preferably use \code{\link{global_envelope_test}} with the suitable \code{type} argument.
#'
#' The directional quantile envelope test (Myllymäki et al., 2015, 2017)
#' takes into account the unequal variances of the test function T(r)
#' for different distances r and is also protected against asymmetry of T(r).
#'
#' @references
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2015). Deviation test construction and power comparison for marked spatial point patterns. Spatial Statistics 11: 19-34. doi: 10.1016/j.spasta.2014.11.004
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' @inheritParams rank_envelope
#' @return An object of class \code{global_envelope} of \code{combined_global_envelope}
#' which can be printed and plotted directly. See \code{\link{global_envelope_test}} for more details.
#' @export
#' @name qdir_envelope
#' @seealso \code{\link{global_envelope_test}}
#' @examples
#' # See more examples in ?global_envelope_test
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' if(require("spatstat", quietly=TRUE)) {
#'   X <- spruces
#'   \donttest{nsim <- 999 # Number of simulations}
#'   \dontshow{nsim <- 19 # Number of simulations for testing}
#'   ## Test for complete spatial randomness (CSR)
#'   # Generate nsim simulations under CSR, calculate centred L-function for the data and simulations
#'   env <- envelope(X, fun="Lest", nsim=nsim, savefuns=TRUE,
#'                   correction="translate", transform = expression(.-r),
#'                   simulate=expression(runifpoint(ex=X)))
#'   res_qdir <- qdir_envelope(env) # The directional quantile envelope test
#'   plot(res_qdir)
#'
#'   ## Advanced use:
#'   # Create a curve set, choosing the interval of distances [r_min, r_max]
#'   curve_set <- crop_curves(env, r_min=1, r_max=8)
#'   # The directional quantile envelope test
#'   res_qdir <- qdir_envelope(curve_set); plot(res_qdir)
#'   # The studentised envelope test
#'   res_st <- st_envelope(curve_set); plot(res_st)
#'   # The unscaled envelope test
#'   res_unscaled <- unscaled_envelope(curve_set); plot(res_unscaled)
#' }
qdir_envelope <- function(curve_set, ...) {
  args <- list(...)
  if("type" %in% names(args)) warning("type is hardcoded to be qdir here. No other options.")
  global_envelope_test(curve_set, type="qdir", ...)
}

#' Studentised envelope test
#'
#' @details The studentised envelope test (Myllymäki et al., 2015, 2017)
#' takes into account the unequal variances of the test function T(r)
#' for different distances r.
#'
#' @export
#' @rdname qdir_envelope
st_envelope <- function(curve_set, ...) {
  args <- list(...)
  if("type" %in% names(args)) warning("type is hardcoded to be st here. No other options.")
  global_envelope_test(curve_set, type="st", ...)
}

#' Unscaled envelope test
#'
#' @details The unscaled envelope test (Ripley, 1981) corresponds to the classical maximum
#' deviation test without scaling, and leads to envelopes with constant width over the distances r.
#' Thus, it suffers from unequal variance of T(r) over the distances r and from the asymmetry of
#' distribution of T(r). We recommend to use the other global envelope tests available,
#' see \code{\link{global_envelope_test}} for full list of alternatives.
#'
#' @references
#' Ripley, B.D. (1981). Spatial statistics. Wiley, New Jersey.
#' @export
#' @rdname qdir_envelope
unscaled_envelope <- function(curve_set, ...) {
  args <- list(...)
  if("type" %in% names(args)) warning("type is hardcoded to be unscaled here. No other options.")
  global_envelope_test(curve_set, type="unscaled", ...)
}
