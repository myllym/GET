# Functionality for central regions based on a curve set
#' @importFrom spatstat fv
individual_central_region <- function(curve_set, type = "erl", coverage = 0.50,
                                      alternative = c("two.sided", "less", "greater"),
                                      probs = c((1-coverage)/2, 1-(1-coverage)/2),
                                      central = "median", ...) {
  if(!is.numeric(coverage) || (coverage < 0 | coverage > 1)) stop("Unreasonable value of coverage.\n")
  alpha <- 1 - coverage
  if(!(type %in% c("rank", "erl", "cont", "area", "qdir", "st", "unscaled")))
    stop("No such type for global envelope.\n")
  alternative <- match.arg(alternative)
  if(type %in% c("qdir", "st", "unscaled") && alternative != "two.sided") {
    warning("For qdir, st and unscaled envelopes only the two.sided alternative is valid.\n")
    alternative <- "two.sided"
  }
  check_probs(probs)
  if(!(central %in% c("mean", "median"))) {
    central <- "median"
    warning("Invalid option fiven for central. Using central = median.\n")
  }
  picked_attr <- pick_attributes(curve_set, alternative=alternative, type=type) # saving for attributes / plotting purposes
  curve_set <- convert_envelope(curve_set)

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
                     alternative=alternative, probs=probs, ...)

  data_and_sim_curves <- data_and_sim_curves(curve_set) # all the functions
  Nfunc <- length(distance) # Number of functions
  nr <- length(curve_set[['r']])
  # Define the central curve T_0
  T_0 <- get_T_0(curve_set)

  #-- Global envelopes
  switch(type,
         rank = {
           #-- the 100(1-alpha)% global rank envelope
           distancesorted <- sort(distance, decreasing=TRUE)
           kalpha <- distancesorted[floor((1-alpha)*(Nfunc))]
           LB <- array(0, nr)
           UB <- array(0, nr)
           for(i in 1:nr){
             Hod <- sort(data_and_sim_curves[,i])
             LB[i]<- Hod[kalpha]
             UB[i]<- Hod[Nfunc-kalpha+1]
           }
         },
         erl = {
           #-- the 100(1-alpha)% global ERL envelope
           distance_lexo_sorted <- sort(distance, decreasing=TRUE)
           kalpha <- distance_lexo_sorted[floor((1-alpha)*Nfunc)]
           curves_for_envelope <- data_and_sim_curves[which(distance >= kalpha),]
           LB <- apply(curves_for_envelope, MARGIN=2, FUN=min)
           UB <- apply(curves_for_envelope, MARGIN=2, FUN=max)
         },
         cont = {
           #-- the 100(1-alpha)% global 'area' envelope (determined similarly as ERL from 'distance')
           distance_cont_sorted <- sort(distance, decreasing=TRUE)
           kalpha <- distance_cont_sorted[floor((1-alpha)*Nfunc)]
           curves_for_envelope <- data_and_sim_curves[which(distance >= kalpha),]
           LB <- apply(curves_for_envelope, MARGIN=2, FUN=min)
           UB <- apply(curves_for_envelope, MARGIN=2, FUN=max)
         },
         area = {
           #-- the 100(1-alpha)% global 'area' envelope (determined similarly as ERL from 'distance')
           distance_area_sorted <- sort(distance, decreasing=TRUE)
           kalpha <- distance_area_sorted[floor((1-alpha)*Nfunc)]
           curves_for_envelope <- data_and_sim_curves[which(distance >= kalpha),]
           LB <- apply(curves_for_envelope, MARGIN=2, FUN=min)
           UB <- apply(curves_for_envelope, MARGIN=2, FUN=max)
         },
         qdir = {
           curve_set_res <- residual(curve_set, use_theo=TRUE)
           quant_m <- curve_set_quant(curve_set_res, probs=probs)
           #-- the 100(1-alpha)% global directional quantile envelope
           distancesorted <- sort(distance)
           kalpha <- distancesorted[floor((1-alpha)*Nfunc)]
           LB <- T_0 - kalpha*abs(quant_m[1,])
           UB <- T_0 + kalpha*abs(quant_m[2,])
         },
         st = {
           sdX <- curve_set_sd(curve_set)
           #-- calculate the 100(1-alpha)% global studentized envelope
           distancesorted <- sort(distance)
           kalpha <- distancesorted[floor((1-alpha)*Nfunc)]
           LB <- T_0 - kalpha*sdX
           UB <- T_0 + kalpha*sdX
         },
         unscaled = {
           #-- calculate the 100(1-alpha)% global unscaled envelope
           distancesorted <- sort(distance)
           kalpha <- distancesorted[floor((1-alpha)*Nfunc)]
           LB <- T_0 - kalpha
           UB <- T_0 + kalpha
         })

  switch(alternative,
         "two.sided" = {},
         "less" = { UB <- Inf },
         "greater" = { LB <- -Inf })

  if(is.vector(curve_set[['obs']])) {
    df <- data.frame(r=curve_set[['r']], obs=curve_set[['obs']], central=T_0, lo=LB, hi=UB)
    picked_attr$einfo$nsim <- Nfunc-1
  }
  else {
    df <- data.frame(r=curve_set[['r']], central=T_0, lo=LB, hi=UB)
    picked_attr$einfo$nsim <- Nfunc
  }
  res <- spatstat::fv(x=df, argu = picked_attr[['argu']],
            ylab = picked_attr[['ylab']], valu = "central", fmla = ". ~ r",
            alim = c(min(curve_set[['r']]), max(curve_set[['r']])),
            labl = picked_attr[['labl']], desc = picked_attr[['desc']],
            unitname = NULL, fname = picked_attr[['fname']], yexp = picked_attr[['yexp']])
  attr(res, "shade") <- c("lo", "hi")
  if(type == "st") picked_attr$einfo$nSD <- kalpha
  if(type == "rank") picked_attr$einfo$nrank <- kalpha
  attr(res, "einfo") <- picked_attr[['einfo']]
  attr(res, "xlab") <- picked_attr[['xlab']]
  attr(res, "xexp") <- picked_attr[['xexp']]
  # Extra for global envelopes
  class(res) <- c("global_envelope", class(res))
  attr(res, "method") <- "Global envelope"
  attr(res, "type") <- type
  attr(res, "k_alpha") <- kalpha
  attr(res, "alpha") <- 1 - coverage
  attr(res, "k") <- distance
  attr(res, "call") <- match.call()
  res
}

# Functionality for global envelope tests based on a curve set (individual central region + p-values)
individual_global_envelope_test <- function(curve_set, type = "erl", alpha = 0.05,
                                            alternative = c("two.sided", "less", "greater"),
                                            ties = "erl", probs = c(0.025, 0.975),
                                            central = "mean", ...) {
  alternative <- match.arg(alternative)
  if(!is.numeric(alpha) || (alpha < 0 | alpha > 1)) stop("Unreasonable value of alpha.\n")
  res <- individual_central_region(curve_set, type=type, coverage=1-alpha,
                                   alternative=alternative, probs=probs,
                                   central=central, ...)
  # The type of the p-value
  possible_ties <- c('midrank', 'random', 'conservative', 'liberal', 'erl')
  if(!(ties %in% possible_ties)) stop("Unreasonable ties argument!\n")

  # Measures for functional ordering
  distance <- attr(res, "k")

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
  if(ntests < 1) stop("Only one curve_set, no combining to be done.\n")
  check_curve_set_dimensions(curve_sets)
  CR_or_GET <- match.arg(CR_or_GET)

  # 1) First stage: Calculate the functional orderings individually for each curve_set
  res_ls <- lapply(curve_sets, FUN = function(x) { individual_central_region(x, ...) })
  type <- attr(res_ls[[1]], "type")

  # 2) Second stage: ERL central region/test
  # Create a curve_set for the ERL test
  k_ls <- lapply(res_ls, FUN = function(x) attr(x, "k"))
  k_mat <- do.call(rbind, k_ls, quote=FALSE)
  Nfunc <- ncol(k_mat)
  # Construct the one-sided ERL central region
  if(type %in% c("qdir", "st", "unscaled")) alt2 <- "greater"
  else alt2 <- "less"
  switch(CR_or_GET,
         CR = {
           curve_set_u <- create_curve_set(list(r=1:ntests, obs=k_mat, is_residual=FALSE))
           res_erl <- individual_central_region(curve_set_u, type="erl", coverage=coverage, alternative=alt2)
         },
         GET = {
           curve_set_u <- create_curve_set(list(r=1:ntests, obs=k_mat[,1], sim_m=k_mat[,-1], is_residual=FALSE))
           res_erl <- individual_global_envelope_test(curve_set_u, type="erl", alpha=1-coverage, alternative=alt2)
         }
  )

  # 3) The 100(1-alpha)% global combined ERL envelope
  distance_lexo_sorted <- sort(attr(res_erl, "k"), decreasing=TRUE)
  kalpha <- distance_lexo_sorted[floor(coverage*Nfunc)]
  # Indices of the curves from which to calculate the convex hull
  curves_for_envelope_ind <- which(attr(res_erl, "k") >= kalpha)
  # Curves
  curve_sets <- lapply(curve_sets, FUN=convert_envelope)
  data_and_sim_curves_l <- lapply(curve_sets, function(x) { data_and_sim_curves(x) })
  # Curves from which to calculate the convex hull
  curves_for_envelope_l <- lapply(data_and_sim_curves_l, function(x) { x[curves_for_envelope_ind,] })
  # Bounding curves
  LB <- lapply(curves_for_envelope_l, FUN = function(x) { apply(x, MARGIN=2, FUN=min) })
  UB <- lapply(curves_for_envelope_l, FUN = function(x) { apply(x, MARGIN=2, FUN=max) })
  # Update the bounding curves (lo, hi) and kalpha to the first level central regions
  for(i in 1:length(curve_sets)) {
    if(attr(res_ls[[i]], "einfo")$alternative != "greater") res_ls[[i]]$lo <- LB[[i]]
    if(attr(res_ls[[i]], "einfo")$alternative != "less") res_ls[[i]]$hi <- UB[[i]]
    attr(res_ls[[i]], "alpha") <- NA
    attr(res_ls[[i]], "k_alpha") <- NULL
  }
  if(!is.null(names(curve_sets))) names(res_ls) <- names(curve_sets)

  # Return
  res <- res_ls
  attr(res, "level2_ge") <- res_erl
  attr(res, "level2_curve_set") <- curve_set_u
  attr(res, "method") <- "Combined global envelope (two-step)"
  class(res) <- c("combined_global_envelope", class(res))
  res
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
  nr <- length(curve_sets[[1]]$r) # all curve sets have the same
  idx <- lapply(1:nfuns, FUN = function(i) ((i-1)*nr+1):(i*nr))
  # Split the envelopes to the original groups
  res_ls <- split(res, f = rep(1:nfuns, each=nr))
  res_ls <- lapply(res_ls, FUN = function(x) { class(x) <- c("global_envelope", class(x)); x })
  # Create empty "level2_ge" attribute containing the test information
  attr(res_ls, "level2_ge") <- data.frame(r=1, obs=attr(res, "k")[1],
                                          central=mean(attr(res, "k")),
                                          lo=attr(res, "k_alpha"), hi=NA)
  mostattributes(attr(res_ls, "level2_ge")) <- attributes(res)
  # Set unreasonable attributes of individuals sets of curves to NULL
  anames <- c("p", "p_interval", "ties", "k", "k_alpha", "method")
  anames <- anames[anames %in% names(attributes(res_ls[[1]]))]
  for(name in anames) {
    for(i in 1:length(res_ls)) attr(res_ls[[i]], name) <- NULL
  }
  attr(res_ls[[i]], "alpha") <- NA
  if(!is.null(curve_sets)) names(res_ls) <- names(curve_sets)
  attr(res_ls, "method") <- "Combined global envelope (one-step)"
  class(res_ls) <- c("combined_global_envelope", class(res_ls))
  res_ls
}

# Helper function for printing object with attributes "alpha", "type", "method"
# and optionally "p", "p_interval", "ties", "alpha_star"
GEprinthelper <- function(x, ...) {
  if(is.null(attr(x, "p"))) { # The case of a central region
    if(inherits(x, c("fboxplot", "combined_fboxplot")))
      cat(attr(x, "method"), " based on ", sep="")
    cat(100*(1-attr(x, "alpha")), "% central region (", attr(x, "type"), "). \n",
        " Plot the object instead.\n", sep="")
  }
  else { # The case of a global envelope test
    cat(attr(x, "method"), " (", attr(x, "type"), ")\n",
        " p-value of the test: ", attr(x, "p"), sep="")
    if(!is.null(attr(x, "ties"))) cat(" (ties method: ", attr(x, "ties"), ")\n", sep="")
    else cat("\n")
    if(!is.null(attr(x, "p_interval")))
      cat(" p-interval         : (", attr(x, "p_interval")[1], ", ", attr(x, "p_interval")[2],")\n", sep="")
  }
  if(!is.null(attr(x, "alpha_star"))) cat(paste("The adjusted level of the test: ", attr(x, "alpha_star"), "\n", sep=""))
}

#' Print method for the class 'global_envelope'
#' @usage \method{print}{global_envelope}(x, ...)
#'
#' @param x an 'global_envelope' object
#' @param ... Ignored.
#'
#' @method print global_envelope
#' @export
print.global_envelope <- function(x, ...) {
  GEprinthelper(x)
}

#' Print method for the class 'global_envelope'
#' @usage \method{print}{combined_global_envelope}(x, ...)
#'
#' @param x an 'combined_global_envelope' object
#' @param ... Ignored.
#'
#' @method print combined_global_envelope
#' @export
print.combined_global_envelope <- function(x, ...) {
  GEprinthelper(attr(x, "level2_ge"))
}

#' Plot method for the class 'global_envelope'
#' @usage \method{plot}{global_envelope}(x, plot_style = c("ggplot2", "fv", "basic"),
#' dotplot = length(x$r)<10,
#' main, ylim, xlab, ylab,
#' color_outside = TRUE, env.col = 1, base_size = 15,
#' labels = NULL, add = FALSE, digits = 3, legend = TRUE, ...)
#' @param x An 'global_envelope' object
#' @param plot_style One of the following "basic", "fv" or "ggplot2".
#' The option "basic" (default) offers a very basic global envelope plot.
#' The option "fv" utilizes the plot routines of the function value table \code{\link[spatstat]{fv.object}}.
#' For "ggplot2", a plot with a coloured envelope ribbon is provided. Requires R library ggplot2.
#' The option "fv" is currently only available for tests with one test function, whereas the other true allow
#' also tests with several tests functions.
#' @param dotplot Logical. If TRUE, then instead of envelopes a dot plot is done.
#' Suitable for low dimensional test vectors. Only applicable if \code{plot_style} is "basic".
#' Default: TRUE if the dimension is less than 10, FALSE otherwise.
#' @param main See \code{\link{plot.default}}. A sensible default exists.
#' @param ylim See \code{\link{plot.default}}. A sensible default exists.
#' @param xlab See \code{\link{plot.default}}. A sensible default exists.
#' @param ylab See \code{\link{plot.default}}. A sensible default exists.
#' @param color_outside Logical. Whether to color the places where the data function goes
#' outside the envelope. Currently red color is used. Relevant only for \code{plot_style = "basic"}.
#' @param env.col The color for the envelope lines (or dotplot arrows). Default 1 (black).
#' @param base_size Base font size, to be passed to theme style when \code{plot_style = "ggplot2"}.
#' @param labels A character vector of suitable length.
#' If \code{dotplot = TRUE}, then labels for the tests at x-axis.
#' @param add Whether to add the plot to an existing plot (TRUE) or to draw a new plot (FALSE).
#' Not available for \code{plot_style = "ggplot2"}.
#' @param digits The number of digits used for printing the p-value or p-interval in the main,
#' if using the default main.
#' @param legend Logical. If FALSE, then the legend is removed from the "ggplot2" style plot.
#' @param ... Additional parameters to be passed to \code{\link{plot}} or \code{\link{lines}}.
#'
#' @method plot global_envelope
#' @export
#' @seealso \code{\link{central_region}}
plot.global_envelope <- function(x, plot_style = c("ggplot2", "fv", "basic"),
                                 dotplot = length(x$r)<10,
                                 main, ylim, xlab, ylab,
                                 color_outside = TRUE, env.col = 1, base_size = 15,
                                 labels = NULL, add = FALSE, digits = 3, legend = TRUE, ...) {
  plot_style <- match.arg(plot_style)
  if(dotplot) plot_style <- "basic"
  # main
  if(missing('main')) {
      main <- env_main_default(x, digits=digits)
  }
  # ylim
  if(missing('ylim')) {
    ylim <- env_ylim_default(x, plot_style == "ggplot2")
  }
  # ylab, ylab, labels
  if(missing('xlab')) {
    if(attr(x, "xlab") == attr(x, "argu")) xlab <- substitute(italic(i), list(i=attr(x, "xexp")))
    else xlab <- substitute(paste(i, " (", italic(j), ")", sep=""), list(i=attr(x, "xexp"), j=attr(x, "argu")))
  }
  if(missing('ylab')) {
    if(inherits(attr(x, "yexp"), "character")) ylab <- attr(x, "yexp")
    else ylab <- substitute(italic(i), list(i=attr(x, "yexp")))
  }
  if(is.null(labels)) if(!is.null(attr(x, "labels"))) labels <- attr(x, "labels")

  switch(plot_style,
         basic = {
           if(dotplot) {
             if(length(labels) != length(x$r)) labels <- NULL
             env_dotplot(x, main=main, ylim=ylim, xlab=xlab, ylab=ylab,
                         color_outside=color_outside, labels=labels,
                         add=add, arrows.col=env.col, ...)
           }
           else {
             env_basic_plot(x, main=main, ylim=ylim, xlab=xlab, ylab=ylab,
                            color_outside=color_outside,
                            add=add, env.col=env.col, ...)
           }
         },
         fv = {
           spatstat::plot.fv(x, main=main, ylim=ylim, xlab=xlab, ylab=ylab, add=add, ...)
         },
         ggplot2 = {
           env_ggplot(x, base_size=base_size, main=main, ylim=ylim, xlab=xlab, ylab=ylab,
                      labels=labels, legend=legend, color_outside=color_outside, ...)
         })
}

#' Plot method for the class 'combined_global_envelope'
#' @usage \method{plot}{combined_global_envelope}(x, plot_style = c("ggplot2", "fv", "basic"),
#'    main, ylim, xlab, ylab,
#'    color_outside = TRUE, env.col = 1, base_size = 15,
#'    labels = NULL, add = FALSE, digits=3,
#'    level = 1, max_ncols_of_plots = 2, nticks = 5,
#'    legend = TRUE, ...)
#' @param x An 'combined_global_envelope' object
#' @inheritParams plot.global_envelope
#' @param labels A character vector of suitable length.
#' If \code{dotplot = TRUE} (for the level 2 test), then labels for the tests at x-axis.
#' Otherwise labels for the separate plots.
#' @param level 1 or 2. In the case of combined tests (with several test functions), two different plots are available:
#' 1 for plotting the combined global envelopes (default and most often wanted) or
#' 2 for plotting the second level test result.
#' @param max_ncols_of_plots The maximum number of columns for the figures. Default 2.
#' (Relates to the number of curve_sets that have been combined.)
#' @param nticks The number of ticks on the xaxis.
#' @method plot combined_global_envelope
#' @export
#' @seealso \code{\link{central_region}}
plot.combined_global_envelope <- function(x, plot_style = c("ggplot2", "fv", "basic"),
                                 main, ylim, xlab, ylab,
                                 color_outside = TRUE, env.col = 1, base_size = 15,
                                 labels = NULL, add = FALSE, digits = 3,
                                 level = 1, max_ncols_of_plots = 2, nticks = 5,
                                 legend = TRUE, ...) {
  plot_style <- match.arg(plot_style)
  if(!(level %in% c(1,2))) stop("Unreasonable value for level.\n")
  # main
  if(missing('main')) {
    alt <- attr(x[[1]], "einfo")$alternative
    main <- env_main_default(attr(x, "level2_ge"), digits=digits, alternative=alt)
  }
  # ylim
  if(missing('ylim'))
    if(level == 1) ylim <- env_ylim_default(x, use_ggplot2 = plot_style == "ggplot2")
    else ylim <- env_ylim_default(attr(x, "level2_ge"), use_ggplot2 = FALSE)
  # ylab, ylab, labels
  if(missing('xlab')) {
    if(plot_style == "ggplot2") {
      xlab <- substitute(italic(i), list(i=attr(attr(x, "level2_ge"), "xexp")))
    }
    else {
      if(attr(x[[1]], "xlab") == attr(x[[1]], "argu")) xlab <- lapply(x, function(y) { substitute(italic(i), list(i=attr(y, "xexp"))) })
      else xlab <- lapply(x, function(y) { substitute(paste(i, " (", italic(j), ")", sep=""), list(i=attr(y, "xexp"), j=attr(y, "argu"))) })
    }
  }
  if(missing('ylab')) {
    if(plot_style == "ggplot2") {
      ylab <- substitute(italic(i), list(i=attr(attr(x, "level2_ge"), "yexp")))
    }
    else {
      if(inherits(attr(x[[1]], "yexp"), "character")) ylab <- lapply(x, function(y) attr(y, "yexp"))
      else ylab <- lapply(x, function(y) { substitute(italic(i), list(i=attr(y, "yexp"))) })
    }
  }
  if(is.null(labels)) {
    if(!is.null(attr(x, "labels"))) labels <- attr(x, "labels")
    else {
      if(!is.null(names(x))) labels <- names(x)
      else {
        if(plot_style == "ggplot2") {
          labels <- sapply(x, function(y) attr(y, "ylab"), simplify=TRUE)
          if(all(sapply(labels, FUN=identical, y=labels[[1]]))) labels <- NULL
        }
      }
    }
  }

  if(level == 1) {
    switch(plot_style,
           basic = {
             env_basic_plot(x, main=labels, ylim=ylim, xlab=xlab, ylab=ylab,
                            color_outside=color_outside,
                            max_ncols_of_plots=max_ncols_of_plots, add=add, env.col=env.col,
                            nticks=nticks, ...)
           },
           fv = {
             n_of_plots <- length(x)
             ncols_of_plots <- min(n_of_plots, max_ncols_of_plots)
             nrows_of_plots <- ceiling(n_of_plots / ncols_of_plots)
             par(mfrow=c(nrows_of_plots, ncols_of_plots))
             if(is.vector(ylim) & length(ylim)==2) ylim <- rep(list(ylim), times=n_of_plots)
             for(i in 1:length(x))
               spatstat::plot.fv(x[[i]],
                                 main=labels[i], ylim=ylim[[i]], xlab=xlab[[i]], ylab=ylab[[i]], add=FALSE, ...)
           },
           ggplot2 = {
             env_ggplot(x, base_size=base_size,
                        main=main, ylim=ylim, xlab=xlab, ylab=ylab,
                        max_ncols_of_plots=max_ncols_of_plots,
                        labels=labels, nticks=nticks, legend=legend, ...)
           })
  }
  else {
     plot.global_envelope(attr(x, "level2_ge"), dotplot = TRUE,
                          main=main, ylim=ylim, xlab=xlab, ylab=ylab,
                          color_outside=color_outside, env.col=env.col, base_size=base_size,
                          labels=labels, add=add, digits=digits, ...)
  }
}

#' Central region / Global envelope
#'
#' Provides central regions or global envelopes or confidence bands
#'
#'
#' Given a \code{curve_set} (see \code{\link{create_curve_set}} for how to create such an object)
#' or an \code{\link[spatstat]{envelope}} object, the function \code{central_region}
#' construcst a central region, i.e. a global envelope, from the given set of functions (or vectors).
#' There are two options for the functions that the \code{curve_set} can contain:
#' \itemize{
#'  \item If the component \code{obs} of the \code{curve_set} is a matrix,
#' then it is assumed that all the functions are data/observed. In this case,
#' the component \code{sim_m} of the \code{curve_set} (which can be then NULL)
#' is ignored and the central region constructed from the functions given in \code{obs}.
#'  \item If the component \code{obs} is a vector, then \code{sim_m} should be provided as well
#' and it is assumed to contain simulated functions (obtained, e.g., from some model or by permutation).
#' The central region is calculated from all the functions.
#' }
#' Thus the \code{curve_set} contains functions (or vectors)
#' \eqn{T_1(r),\dots,T_s(r)}{T_1(r),...,T_s(r)}.
#' In the case of one observed function only,
#' the data function is considered to be \eqn{T_1(r)}{T_1(r)}.
#'
#' Generally an envelope is a band bounded by the vectors (or functions)
#' \eqn{T_{\text{low}}}{T_lo} and \eqn{T_{\text{hi}}}{T_hi}.
#' A \eqn{100(1-\alpha)}{100(1-alpha)}\% or 100*coverage\% global envelope is a set
#' \eqn{(T_{\text{low}}, T_{\text{hi}})}{(T_lo, T_hi)} of envelope vectors
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
#' can be found in \code{\link{forder}}.
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
#' We note that the global envelopes \code{'rank'}, \code{'erl'}, \code{'cont'} and \code{'area'}
#' are completely non-parametric tests and thus protected against the unequal variances
#' of the test function T(r) for different distances r and also against asymmetry of the distribution
#' of T(r).
#'
#' For each curve in the curve_set, both the data curve and the simulations,
#' an above mention measure k is determined. The measure values
#' \eqn{k_1, k_2, ..., k_s}{k_1, k_2, ..., k_s}
#' are returned in the attribute 'k' (in a case of one observed curve only, k[1] is its value).
#' Based on the chosen measure, the central region, i.e. the global envelope, is constructed
#' on the chosen interval of argument values (the functions in the \code{curve_set} are assumed
#' to be given on this interval only).
#'
#' @references
#' Mrkvička, T., Hahn, U. and Myllymäki, M. (2018). A one-way ANOVA test for functional data with graphical interpretation. arXiv:1612.03608 [stat.ME]
#'
#' Mrkvička, T., Myllymäki, M. and Narisetty, N. N. (2019) New methods for multiple testing in permutation inference for the general linear model.
#'
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2015). Deviation test construction and power comparison for marked spatial point patterns. Spatial Statistics 11: 19-34. doi: 10.1016/j.spasta.2014.11.004
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
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
#' @param ... Additional parameters to be passed to \code{\link{forder}}.
#' @return An object of class "global_envelope" and "fv" (see \code{\link[spatstat]{fv.object}}),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = the vector of values of the argument r at which the test was made
#' \item obs = the data function, if there is only one data function. Otherwise not existing.
#' \item lo = the lower envelope based on the simulated functions
#' \item hi = the upper envelope based on the simulated functions
#' \item central = If the curve_set (or envelope object) contains a component 'theo',
#'       then this function is used as the central curve and returned in this component.
#'       Otherwise, the central_curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'       Used for visualization only.
#' }
#' Additionally, the return value has attributes
#' \itemize{
#'   \item method = The name of the envelope test used for plotting purposes ("Global envelope")
#'   \item alternative = The alternative specified in the function call.
#'   \item ties = As the argument \code{ties}.
#'   \item k_alpha = The value of k corresponding to the 100(1-alpha)\% global envelope.
#'   \item k = The values of the chosen measure for all the functions. If there is only one
#'   observed function, then k[1] will give the value of the measure for this.
#'   \item call = The call of the function.
#' }
#' and a punch of attributes for the "fv" object type, see \code{\link[spatstat]{fv}}.
#' Attributes of an object \code{res} can be obtained using the function
#' \code{\link[base]{attr}}, e.g. \code{attr(res, "k")} for the values of the ordering measure.
#' @export
#' @seealso \code{\link{global_envelope_test}}
#' @aliases global_envelope
#' @examples
#' ## A central region of a set of functions
#' #----------------------------------------
#' if(requireNamespace("fda", quietly = TRUE)) {
#'   curve_set <- create_curve_set(list(r=as.numeric(row.names(fda::growth$hgtf)),
#'                                      obs=fda::growth$hgtf))
#'   plot(curve_set, ylab="height")
#'   cr <- central_region(curve_set, coverage=0.50, type="erl")
#'   plot(cr, main="50% central region")
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
#' boot.cset <- create_curve_set(list(r=1:length(x), obs=t(m)))
#' cr <- central_region(boot.cset, coverage=0.95, type="erl")
#' CB.lo <- ftheta+s0*cr$lo
#' CB.hi <- ftheta+s0*cr$hi
#'
#' # Plotting the result
#' plot(d, ylab="f(x)", xaxt="n", xlab="x", main="95% central region")
#' axis(1, at=(0:5)*20, labels=(0:5)/5)
#' lines(ftheta)
#' lines(CB.lo, lty=2)
#' lines(CB.hi, lty=2)
central_region <- function(curve_sets, type = "erl", coverage = 0.50,
                           alternative = c("two.sided", "less", "greater"),
                           probs = c((1-coverage)/2, 1-(1-coverage)/2),
                           central = "median", nstep = 2, ...) {
  if(class(curve_sets)[1] == "list" & length(curve_sets) == 1) curve_sets <- curve_sets[[1]]
  if(class(curve_sets)[1] == "list") {
    if(!(nstep %in% c(1,2))) stop("Invalid number of steps (nstep) for combining. Should be 1 or 2.")
    if(nstep == 2) {
      res <- combined_CR_or_GET(curve_sets, CR_or_GET="CR", type=type, coverage=coverage,
                                alternative=alternative, probs=probs,
                                central=central, ...)
    }
    else { # One-step combining procedure
      res <- combined_CR_or_GET_1step(curve_sets, CR_or_GET="CR", type=type, coverage=coverage,
                                      alternative=alternative, probs=probs,
                                      central=central, ...)
    }
  }
  else {
    res <- individual_central_region(curve_sets, type=type, coverage=coverage,
                                     alternative=alternative, probs=probs,
                                     central=central, ...)
  }
  res
}

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
#'   curve_set <- create_curve_set(list(r=as.numeric(row.names(fda::growth$hgtf)),
#'                                      obs=fda::growth$hgtf))
#'   plot(curve_set, ylab="height")
#'   bp <- fBoxplot(curve_set, coverage=0.50, type="area", factor=1)
#'   plot(bp)
#' }
fBoxplot <- function(curve_sets, factor = 1.5, ...) {
  res <- central_region(curve_sets, ...)
  if(inherits(res, "combined_global_envelope")) {
    dist <- factor * (attr(res, "level2_ge")$hi - attr(res, "level2_ge")$lo)
    attr(attr(res, "level2_ge"), "whisker.lo") <- attr(res, "level2_ge")$lo - dist
    attr(attr(res, "level2_ge"), "whisker.hi") <- attr(res, "level2_ge")$hi + dist
    attr(attr(res, "level2_ge"), "method") <- "Functional boxplot"
    class(attr(res, "level2_ge")) <- c("fboxplot", class(attr(res, "level2_ge")))
    for(i in 1:length(res)) {
      dist <- factor * (res[[i]]$hi - res[[i]]$lo)
      attr(res[[i]], "whisker.lo") <- res[[i]]$lo - dist
      attr(res[[i]], "whisker.hi") <- res[[i]]$hi + dist
      attr(res[[i]], "method") <- "Functional boxplot"
      class(res[[i]]) <- c("fboxplot", class(res[[i]]))
    }
    attr(res, "curve_sets") <- curve_sets
    attr(res, "factor") <- factor
    attr(res, "method") <- "Combined functional boxplot"
    attr(res, "call") <- match.call()
    class(res) <- c("combined_fboxplot", class(res))
  }
  else {
    dist <- factor * (res$hi - res$lo)
    attr(res, "whisker.lo") <- res$lo - dist
    attr(res, "whisker.hi") <- res$hi + dist
    attr(res, "method") <- "Functional boxplot"
    attr(res, "curve_sets") <- curve_sets
    attr(res, "factor") <- factor
    attr(res, "call") <- match.call()
    class(res) <- c("fboxplot", class(res))
  }
  res
}

#' Plot method for the class 'fboxplot'
#' @usage \method{plot}{fboxplot}(x, plot_style = c("ggplot2", "fv", "basic"),
#'    dotplot = length(x$r)<10,
#'    outliers = TRUE, bp.col = 2, cr.col = 1, ...)
#'
#' @param x an 'fboxplot' object
#' @inheritParams plot.global_envelope
#' @param outliers Logical. If TRUE, then the functions outside the functional boxplot are drawn.
#' @param bp.col The color for the boxplot bounds. Default 2 (red).
#' @param cr.col The color for the central region bounds.
#' @param ... Additional arguments to be passed to \code{\link{plot.global_envelope}}.
#'
#' @method plot fboxplot
#' @export
plot.fboxplot <- function(x, plot_style = c("ggplot2", "fv", "basic"),
                          dotplot = length(x$r)<10,
                          outliers = TRUE, bp.col = 2, cr.col = 1, ...) {
  plot_style <- match.arg(plot_style)
  if(outliers) curve_sets <- attr(x, "curve_sets") else curve_sets <- NULL
  cr <- x
  x$lo <- attr(x, "whisker.lo"); x$hi <- attr(x, "whisker.hi") # Functional boxplot

  if(retick_xaxis(x)$retick_xaxis) {
    if(plot_style == "fv") {
      warning("The plot style fv not available for the case where r distances are not increasing.\n Setting plot_style to basic.\n")
      plot_style <- "basic"
    }
  }
  switch(plot_style,
         basic = {
           if(dotplot) {
             # Functional boxplot
             plot.global_envelope(x, plot_style=plot_style, env.col=bp.col, dotplot=TRUE, ..., curve_sets=curve_sets)
             # Central region
             plot.global_envelope(cr, plot_style=plot_style, env.col=cr.col, dotplot=TRUE, add=TRUE, ..., curve_sets=NULL)
           }
           else {
             plot.global_envelope(x, plot_style=plot_style, env.col=bp.col, dotplot=FALSE, ..., curve_sets=curve_sets)
             # Central region
             lines(cr$r, cr$lo, lty=2, col=cr.col)
             lines(cr$r, cr$hi, lty=2, col=cr.col)
           }
         },
         fv = {
           plot.global_envelope(x, plot_style=plot_style, env.col=bp.col, ..., curve_sets=NULL)
           # Outliers
           for(i in 1:ncol(attr(x, "curve_set")$obs)) {
             if(any(curve_sets$obs[,i] < x$lo | curve_sets$obs[,i] > x$hi))
               lines(curve_sets$r, curve_sets$obs[,i], col=grey(0.5))
           }
           # Central region
           lines(cr$r, cr$lo, lty=2, col=cr.col)
           lines(cr$r, cr$hi, lty=2, col=cr.col)
         },
         ggplot2 = {
           plot.global_envelope(x, plot_style=plot_style, env.col=bp.col, ..., curve_sets=curve_sets, x2=cr)
         })
}

#' Plot method for the class 'combined_fboxplot'
#' @usage \method{plot}{combined_fboxplot}(x, plot_style = c("ggplot2", "fv", "basic"), level = 1,
#'    outliers = TRUE, bp.col = 2, cr.col = 1, ...)
#'
#' @param x an 'combined_fboxplot' object
#' @inheritParams plot.combined_global_envelope
#' @inheritParams plot.fboxplot
#' @param ... Additional arguments to be passed to \code{\link{plot.combined_global_envelope}}.
#'
#' @method plot combined_fboxplot
#' @export
plot.combined_fboxplot <- function(x, plot_style = c("ggplot2", "fv", "basic"), level = 1,
                          outliers = TRUE, bp.col = 2, cr.col = 1, ...) {
  plot_style <- match.arg(plot_style)
  if(!(level %in% c(1,2))) stop("Unreasonable value for level.\n")
  if(outliers) curve_sets <- attr(x, "curve_sets") else curve_sets <- NULL
  cr <- x
  attr(x, "level2_ge")$lo <- attr(attr(x, "level2_ge"), "whisker.lo")
  attr(x, "level2_ge")$hi <- attr(attr(x, "level2_ge"), "whisker.hi") # Functional boxplot

  # Combined test, level 1 plots
  if(level == 1) {
    # Set also first level bounds of x to whiskers
    for(i in 1:length(x)) {
      x[[i]]$lo <- attr(x[[i]], "whisker.lo")
      x[[i]]$hi <- attr(x[[i]], "whisker.hi")
    }
    switch(plot_style,
           basic = {
             plot.combined_global_envelope(x, plot_style=plot_style, env.col=bp.col, ..., curve_sets=curve_sets)
           },
           fv = {
             plot.combined_global_envelope(x, plot_style=plot_style, env.col=bp.col, ...)
             if(outliers) cat("For fv style & combined boxplots, plotting outliers is not implemented.\n")
           },
           ggplot2 = {
             plot.combined_global_envelope(x, plot_style=plot_style, env.col=bp.col, ..., curve_sets=curve_sets, x2=cr)
           })
  }
  else {
    plot.fboxplot(x, plot_style=plot_style, outliers=outliers, bp.col=bp.col, cr.col=cr.col, ...)
  }
}


#' Global envelope test
#'
#' Global envelope test, p-values and global envelopes
#'
#'
#' Given a \code{curve_set} (see \code{\link{create_curve_set}} for how to create such an object)
#' or an \code{\link[spatstat]{envelope}} object,
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
#' using single functions (for simple hypotheses).
#' Different \code{type} of global envelope tests can be performed.
#' We use such ordering of the functions for which we are able to construct global envelopes
#' with intrinsic graphical interpretation.
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
#' See \code{\link{forder}} and \code{\link{central_region}} and the references
#' for more detailed description of the measures and the corresponding envelopes.
#'
#' @section Ranking of the curves:
#' The options for measures to order the functions from the most extreme one to the least extreme one
#' are given by the argument \code{type}: 'rank', 'erl', 'cont', 'area', 'qdir', 'st', 'unscaled'.
#' The options are
#' \itemize{
#' \item 'rank': extreme ranks (Myllymäki et al., 2017)
#' \item 'erl': extreme rank lengths (Myllymäki et al., 2017; Mrkvička et al., 2018)
#' \item 'cont': continuous ranks (Hahn, 2015; Mrkvička et al., 2019)
#' \item 'area': area ranks (Mrkvička et al., 2019)
#' \item 'qdir': the directional quantile maximum absolute deviation (MAD) measure (Myllymäki et al., 2015, 2017)
#' \item 'st': the studentized MAD measure (Myllymäki et al., 2015, 2017)
#' \item 'unscaled': the unscaled MAD measure (Ripley, 1981)
#' }
#' See more detailed description of the envelopes and measures in \code{\link{central_region}}
#' and \code{\link{forder}}.
#'
#' @section Global envelope:
#' Based on the measures used to rank the functions, the 100(1-alpha)\% global envelope is provided.
#' It corresponds to the 100*coverage\% central region.
#'
#' @section P-values:
#' In the case \code{type="rank"}, based on the extreme ranks k_i, i=1, ..., s+1,
#' the p-interval is calculated. Because the extreme ranks contain ties, there is not just
#' one p-value. The p-interval is given by the most liberal and the most conservative p-value
#' estimate. Also a single p-value is calculated.
#' By default this single p-value is the extreme rank length p-value ("erl"),
#' but another option can be used by specifying \code{ties} argument.
#'
#' If the case \code{type = "erl"}, the (single) p-value based on the extreme rank length ordering
#' of the functions is calculated and returned in the attribute \code{p}.
#' The same is done for other measures, the p-value always being correspondent to the chosen measure.
#'
#' @section Number of simulations:
#' The tests \code{'erl'}, \code{'cont'} and \code{'area'}, similarly as
#' the MAD deviation/envelope tests \code{'qdir'}, \code{'st'} and \code{'unscaled'},
#' allow in principle a lower number of simulations to be used than the test based on
#' extreme ranks (\code{'rank'}), because no ties occur for these measures.
#' However, if affordable, we recommend some thousands of simulations in any case
#' to achieve a good power and repeatability of the test.
#'
#' @references
#' Mrkvička, T., Myllymäki, M. and Hahn, U. (2017). Multiple Monte Carlo testing, with applications in spatial point processes. Statistics & Computing 27 (5): 1239-1255. doi: 10.1007/s11222-016-9683-9
#'
#' Mrkvička, T., Myllymäki, M., Jilek, M., and Hahn, U. (2018). A one-way ANOVA test for functional data with graphical interpretation. arXiv:1612.03608 [stat.ME]
#'
#' Mrkvička, T., Myllymäki, M. and Narisetty, N. N. (2019) New methods for multiple testing in permutation inference for the general linear model.
#'
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2015). Deviation test construction and power comparison for marked spatial point patterns. Spatial Statistics 11: 19-34. doi: 10.1016/j.spasta.2014.11.004
#'
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' Ripley, B.D. (1981). Spatial statistics. Wiley, New Jersey.
#'
#' @inheritParams central_region
#' @param curve_sets A \code{curve_set} (see \code{\link{create_curve_set}})
#' or an \code{\link[spatstat]{envelope}} object containing a data function and simulated functions.
#' If an envelope object is given, it must contain the summary
#' functions from the simulated patterns which can be achieved by setting
#' \code{savefuns = TRUE} when calling \code{\link[spatstat]{envelope}}.
#' Alternatively, a list of \code{curve_set} or \code{\link[spatstat]{envelope}} objects can be given.
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
#' @return An object of class "global_envelope" and "fv"
#' (see \code{\link[spatstat]{fv.object}}), which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = the vector of values of the argument r at which the test was made
#' \item obs = values of the data function
#' \item lo = the lower envelope based on the simulated functions
#' \item hi = the upper envelope based on the simulated functions
#' \item central = If the curve_set (or envelope object) contains a component 'theo',
#'       then this function is used as the central curve and returned in this component.
#'       Otherwise, the central_curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'       Used for visualization only.
#' }
#' Moreover, the return value has the same attributes as the object returned by
#' \code{\link{central_region}} and in addition
#' \itemize{
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#' }
#' and in the case that \code{type = 'rank'} also
#' \itemize{
#'   \item p_interval = The p-value interval [p_liberal, p_conservative].
#'   \item ties = As the argument \code{ties}.
#' }
#' @export
#' @seealso \code{\link{plot.global_envelope}}, \code{\link{central_region}},
#' \code{\link{global_envelope_test2d}}
#' @aliases GET
#' @examples
#' if(require(spatstat, quietly=TRUE)) {
#'   ## Testing complete spatial randomness (CSR)
#'   #-------------------------------------------
#'   \donttest{nsim <- 1999 # Number of simulations}
#'   \dontshow{nsim <- 19 # Number of simulations}
#'   pp <- unmark(spruces)
#'   # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#'   env <- envelope(pp, fun="Lest", nsim=nsim,
#'                   savefuns=TRUE, # save the functions
#'                   correction="translate", # edge correction for L
#'                   simulate=expression(runifpoint(ex=pp))) # Simulate CSR
#'   # The rank envelope test (extreme rank length (ERL) breaking of ties)
#'   res <- global_envelope_test(env)
#'   # Plot the result.
#'   # - The central curve is now obtained from env[['theo']], which is the
#'   # value of the L-function under the null hypothesis (L(r) = r).
#'   # - Three different plot styles are provided:
#'   # a) ggplot2 style (requires ggplot2)
#'   plot(res, plot_style="ggplot2")
#'   # b) spatstat's style (requires spatstat)
#'   plot(res, plot_style="fv")
#'   # c) a basic style
#'   plot(res, plot_style="basic")
#'
#'   ## Advanced use:
#'   # Choose the interval of distances [r_min, r_max] (at the same time create a curve_set from 'env')
#'   curve_set <- crop_curves(env, r_min=1, r_max=7)
#'   # For better visualisation, take the L(r)-r function
#'   curve_set <- residual(curve_set, use_theo=TRUE)
#'   # Do the rank envelope test (erl)
#'   res <- global_envelope_test(curve_set); plot(res, ylab=expression(italic(L(r)-r)))
#'
#'   ## Random labeling test
#'   #----------------------
#'   mpp <- spruces
#'   # 1) Perform simulations under the random labelling hypothesis and calculate
#'   # the test function T(r) for the data pattern (mpp) and each simulation.
#'   # The command below specifies that the test function is T(r) = \hat{L}_mm(r),
#'   # which is an estimator of the mark-weighted L function, L_mm(r),
#'   # with translational edge correction.
#'   \donttest{nsim <- 1999 # Number of simulations}
#'   \dontshow{nsim <- 19 # Number of simulations}
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
#'   plot(res, ylab=expression(italic(L[m](r)-L(r))))
#'
#'   ## Goodness-of-fit test (typically conservative, see dg.global_envelope for adjusted tests)
#'   #-----------------------------------------------
#'   pp <- unmark(spruces)
#'   # Minimum distance between points in the pattern
#'   min(nndist(pp))
#'   # Fit a model
#'   fittedmodel <- ppm(pp, interaction=Hardcore(hc=1)) # Hardcore process
#'
#'   # Simulating Gibbs process by 'envelope' is slow, because it uses the MCMC algorithm
#'   #env <- envelope(fittedmodel, fun="Jest", nsim=999, savefuns=TRUE,
#'   #                correction="none", r=seq(0, 4, length=500))
#'
#'   # Using direct algorihm can be faster, because the perfect simulation is used here.
#'   simulations <- NULL
#'   \donttest{nsim <- 999 # Number of simulations}
#'   \dontshow{nsim <- 19 # Number of simulations}
#'   for(j in 1:nsim) {
#'      simulations[[j]] <- rHardcore(beta=exp(fittedmodel$coef[1]),
#'                                    R=fittedmodel$interaction$par$hc,
#'                                    W=pp$window)
#'      if(j%%10==0) cat(j, "...", sep="")
#'   }
#'   env <- envelope(pp, simulate=simulations, fun="Jest", nsim=length(simulations),
#'                   savefuns=TRUE, correction="none", r=seq(0, 4, length=500))
#'   curve_set <- crop_curves(env, r_min=1, r_max=3.5)
#'   res <- global_envelope_test(curve_set, type="erl"); plot(res, ylab=expression(italic(J(r))))
#'
#'   # A combined global envelope test
#'   #--------------------------------
#'   # As an example test CSR of the saplings point pattern by means of
#'   # L, F, G and J functions.
#'   data(saplings)
#'   X <- saplings
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
#'   res <- global_envelope_test(curve_sets=list(curve_set_L, curve_set_F,
#'                                               curve_set_G, curve_set_J))
#'   plot(res, labels=c("L(r)-r", "F(r)", "G(r)", "J(r)"))
#' }
#'
#' ## A test based on a low dimensional random vector
#' #-------------------------------------------------
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
#'
global_envelope_test <- function(curve_sets, type = "erl", alpha = 0.05,
                          alternative = c("two.sided", "less", "greater"),
                          ties = "erl", probs = c(0.025, 0.975),
                          central = "mean", nstep = 2, ...) {
  if(class(curve_sets)[1] == "list" & length(curve_sets) == 1) curve_sets <- curve_sets[[1]]
  if(class(curve_sets)[1] == "list") {
    if(!(nstep %in% c(1,2))) stop("Invalid number of steps (nstep) for combining. Should be 1 or 2.")
    if(nstep == 2) {
      res <- combined_CR_or_GET(curve_sets, CR_or_GET="GET", type=type, coverage=1-alpha,
                                alternative=alternative, probs=probs,
                                central=central, ...)
    }
    else { # One-step combining procedure
      res <- combined_CR_or_GET_1step(curve_sets, CR_or_GET="GET", type=type, coverage=1-alpha,
                                alternative=alternative, probs=probs,
                                central=central, ...)
    }
  }
  else {
    res <- individual_global_envelope_test(curve_sets, type=type, alpha=alpha,
                                           alternative=alternative,
                                           ties=ties, probs=probs,
                                           central=central, ...)
  }
  res
}


#' The rank envelope test
#'
#' The rank envelope test, p-values and global envelopes
#'
#'
#' The rank envelope test is a completely non-parametric test, which provides
#' the 100(1-alpha)\% global envelope for the chosen test function T(r) on
#' the chosen interval of distances and associated p-values.
#'
#' The test corresponds to the global envelope test that can be carriet out by
#' \code{\link{global_envelope_test}} by specifying the \code{type} for which the options
#' \code{"rank"}, \code{"erl"}, \code{"cont"} and \code{"area"} are available. The last
#' three are modifications of the first one to treat the ties in the extreme rank ordering
#' used in \code{"rank"}.
#'
#' Note: Earlier it was possible to specify to the extreme rank lengths breaking of ties for the rank
#' envelope with specifying the argument \code{lexo = TRUE}. This is obsolete now. The same can be done
#' by choosing \code{type = "rank"} and \code{ties = "erl"}, which is in fact the default of this
#' \code{rank_envelope} function.
#'
#' @section Global envelope:
#' The 100(1-alpha)\% global envelope is provided in addition to the p-values.
#' If \code{type = "rank"} then the envelope is the global rank envelope proposed by
#' Myllymäki et al. (2017).
#' If \code{type = "erl"} then the envelope is the global rank envelope based on the
#' extreme rank length ordering. This envelope is constructed as the convex hull of
#' the functions which have extreme rank length measure \eqn{R_i^{\text{erl}}}{Rerl_i}
#' that is larger or equal to the critical \eqn{\alpha}{alpha} level of the extreme rank
#' length measure (Mrkvička et al., 2018).
#'
#' @section Number of simulations:
#' The extreme rank length ordering test (\code{type = "erl"}) allows in principle a lower numbe
#' of simulations to be used than the test based on extreme ranks (\code{type = "rank"}).
#' However, we recommend some thousands of simulations in any case to achieve a good power
#' and repeatability of the test.
#'
#' @references
#' Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for spatial point patterns. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 79: 381–404. doi: 10.1111/rssb.12172
#'
#' Mrkvička, T., Myllymäki, M. and Hahn, U. (2017). Multiple Monte Carlo testing, with applications in spatial point processes. Statistics & Computing 27 (5): 1239-1255. doi: 10.1007/s11222-016-9683-9
#'
#' Mrkvička, T., Myllymäki, M., Jilek, M. and Hahn, U. (2018). A one-way ANOVA test for functional data with graphical interpretation. arXiv:1612.03608 [stat.ME]
#'
#' @param curve_set A curve_set (see \code{\link{create_curve_set}}) or an \code{\link[spatstat]{envelope}}
#'  object. If an envelope object is given, it must contain the summary
#'  functions from the simulated patterns which can be achieved by setting
#'  savefuns = TRUE when calling \code{\link[spatstat]{envelope}}.
#' @param type The type of the global envelope with current options for "rank", "erl", "cont" and "area".
#' If "rank", the global rank envelope accompanied by the p-interval is given (Myllymäki et al., 2017).
#' If "erl", the global rank envelope based on extreme rank lengths accompanied by the extreme rank
#' length p-value is given (Myllymäki et al., 2017, Mrkvička et al., 2018). See details and additional
#' sections thereafter.
#' @param ... Additional parameters to be passed to \code{\link{global_envelope_test}}.
#' @return An object of class "global_envelope" and "fv" (see \code{\link[spatstat]{fv.object}}),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = the vector of values of the argument r at which the test was made
#' \item obs = values of the test function for the data point pattern
#' \item lo = the lower envelope based on the simulated functions
#' \item hi = the upper envelope based on the simulated functions
#' \item central = If the curve_set (or envelope object) contains a component 'theo',
#'       then this function is used as the central curve and returned in this component.
#'       Otherwise, the central_curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'       Used for visualization only.
#' }
#' Additionally, the return value has attributes
#' \itemize{
#'   \item method = The name of the envelope test ("Rank envelope test" for the rank envelope test)
#'   \item alternative = The alternative specified in the function call.
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#'   \item p_interval = The p-value interval [p_liberal, p_conservative].
#'   \item ties = As the argument \code{ties}.
#'   \item k_alpha = The value of k corresponding to the 100(1-alpha)\% global envelope.
#'   \item k = Global rank values (for type="rank") or extreme rank lengths (for type="erl").
#'   k[1] is the value for the data pattern.
#'   \item call = The call of the function.
#' }
#' and a punch of attributes for the "fv" object type, see \code{\link[spatstat]{fv}}.
#' @export
#' @seealso \code{\link{global_envelope_test}}, \code{\link{plot.global_envelope}}
#' @examples
#' # See ?global_envelope_test for more examples
#'
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' if(require(spatstat, quietly=TRUE)) {
#'   pp <- unmark(spruces)
#'   \donttest{nsim <- 2499 # Number of simulations}
#'   \dontshow{nsim <- 19 # Number of simulations for testing}
#'   # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#'   env <- envelope(pp, fun="Lest", nsim=nsim, savefuns=TRUE, correction="translate",
#'                   simulate=expression(runifpoint(ex=pp)))
#'   # The rank envelope test
#'   res <- rank_envelope(env)
#'   # Plot the result.
#'   # - The central curve is now obtained from env[['theo']], which is the
#'   # value of the L-function under the null hypothesis (L(r) = r).
#'   plot(res)
#'   # or (requires R library ggplot2)
#'   plot(res, plot_style="ggplot2")
#'
#'   ## Advanced use:
#'   # Choose the interval of distances [r_min, r_max] (at the same time create a curve_set from 'env')
#'   curve_set <- crop_curves(env, r_min=1, r_max=7)
#'   # For better visualisation, take the L(r)-r function
#'   curve_set <- residual(curve_set, use_theo=TRUE)
#'   # Do the rank envelope test
#'   res <- rank_envelope(curve_set); plot(res, plot_style="ggplot2")
#' }
rank_envelope <- function(curve_set, type = "rank", ...) {
  if(!(type %in% c("rank", "erl", "cont", "area"))) stop("No such type for the global rank envelope.\n")
  global_envelope_test(curve_set, type=type, ...)
}

#' Global scaled maximum absolute difference (MAD) envelope tests
#'
#' Performs the global scaled MAD envelope tests, either directional quantile or studentised,
#' or the unscaled MAD envelope test. These tests correspond to calling the
#' function \code{\link{global_envelope_test}} with \code{"qdir"}, \code{type = "st"} and
#' \code{"unscaled"}, respectively. The functions \code{qdir_envelope}, \code{st_envelope} and
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
#' @return An object of class "global_envelope" and "fv" (see \code{\link[spatstat]{fv.object}}),
#' which can be printed and plotted directly.
#'
#' Essentially a data frame containing columns
#' \itemize{
#' \item r = the vector of values of the argument r at which the test was made
#' \item obs = values of the test function for the data point pattern
#' \item lo = the lower envelope based on the simulated functions
#' \item hi = the upper envelope based on the simulated functions
#' \item central = If the curve_set (or envelope object) contains a component 'theo',
#'       then this function is used as the central curve and returned in this component.
#'       Otherwise, the central curve is the mean of the test functions T_i(r), i=2, ..., s+1.
#'       Used for visualization only.
#' }
#' Additionally, the return value has attributes
#' \itemize{
#'   \item method = The name of the method ("Global envelope test")
#'   \item alternative = "two-sided
#'   \item p = A point estimate for the p-value (default is the mid-rank p-value).
#'   \item u_alpha = The value of u corresponding to the 100(1-alpha)\% global envelope.
#'   \item u = Deviation values (u[1] is the value for the data pattern).
#'   \item call = The call of the function.
#' }
#' and a punch of attributes for the "fv" object type.
#' @export
#' @name qdir_envelope
#' @seealso \code{\link{global_envelope_test}}, \code{\link{plot.global_envelope}},
#' \code{\link{global_envelope_test2d}}, \code{\link{dg.global_envelope_test}}
#' @examples
#' # See more examples in ?global_envelope_test
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' if(require("spatstat", quietly=TRUE)) {
#'   pp <- spruces
#'   \donttest{nsim <- 999 # Number of simulations}
#'   \dontshow{nsim <- 19 # Number of simulations for testing}
#'   ## Test for complete spatial randomness (CSR)
#'   # Generate nsim simulations under CSR, calculate L-function for the data and simulations
#'   env <- envelope(pp, fun="Lest", nsim=nsim, savefuns=TRUE, correction="translate",
#'                   simulate=expression(runifpoint(ex=pp)))
#'   res_qdir <- qdir_envelope(env) # The directional quantile envelope test
#'   plot(res_qdir)
#'   # or (requires R library ggplot2)
#'   plot(res_qdir, plot_style="ggplot2")
#'
#'   ## Advanced use:
#'   # Create a curve set, choosing the interval of distances [r_min, r_max]
#'   curve_set <- crop_curves(env, r_min=1, r_max=8)
#'   # For better visualisation, take the L(r)-r function
#'   curve_set <- residual(curve_set, use_theo=TRUE)
#'   # The directional quantile envelope test
#'   res_qdir <- qdir_envelope(curve_set); plot(res_qdir, plot_style="ggplot2")
#'   # The studentised envelope test
#'   res_st <- st_envelope(curve_set); plot(res_st, plot_style="ggplot2")
#'   # The unscaled envelope test
#'   res_unscaled <- unscaled_envelope(curve_set); plot(res_unscaled, plot_style="ggplot2")
#' }
qdir_envelope <- function(curve_set, ...) {
  args <- list(...)
  if("type" %in% names(args)) warning("type is hardcoded to be qdir here. No other options.\n")
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
  if("type" %in% names(args)) warning("type is hardcoded to be st here. No other options.\n")
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
  if("type" %in% names(args)) warning("type is hardcoded to be unscaled here. No other options.\n")
  global_envelope_test(curve_set, type="unscaled", ...)
}
