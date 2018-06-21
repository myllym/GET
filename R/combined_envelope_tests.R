# Functionality for combined_central_region and combined_global_envelope_test
combined_CR_or_GET <- function(curve_sets, CR_or_GET = c("CR", "GET"), coverage, ...) {
  ntests <- length(curve_sets)
  curve_sets <- check_curve_set_dimensions(curve_sets)

  # 1) First stage: Calculate the functional orderings individually for each curve_set
  # switch(CR_or_GET,
  #        CR = {
           res_ls <- lapply(curve_sets, FUN = function(x) { central_region(x, savedevs=TRUE, ...) })
         # },
         # GET = {
         #   res_ls <- lapply(curve_sets, FUN = function(x) { global_envelope_test(x, savedevs=TRUE, ...) })
         # })
  type <- attr(res_ls[[1]], "type")

  # 2) Second stage: ERL central region
  # Create a curve_set for the ERL test
  k_ls <- lapply(res_ls, FUN = function(x) attr(x, "k"))
  k_mat <- do.call(rbind, k_ls, quote=FALSE)
  curve_set_u <- create_curve_set(list(r=1:ntests, obs=k_mat, is_residual=FALSE))
  # Construct the one-sided ERL central region
  if(type %in% c("qdir", "st", "unscaled")) alt2 <- "greater"
  else alt2 <- "less"
  switch(CR_or_GET,
    CR = {
      res_erl <- central_region(curve_set_u, type="erl", coverage=coverage, savedevs=TRUE, alternative=alt2)
    },
    GET = {
      res_erl <- global_envelope_test(curve_set_u, type="erl", alpha=1-coverage, savedevs=TRUE, alternative=alt2)
    }
  )

  # 3) The 100(1-alpha)% global combined ERL envelope
  distance_lexo_sorted <- sort(attr(res_erl, "k"), decreasing=TRUE)
  Nfunc <- ncol(curve_set_u$obs)
  kalpha <- distance_lexo_sorted[floor(coverage*Nfunc)]
  # Indices of the curves from which to calculate the convex hull
  curves_for_envelope_ind <- which(attr(res_erl, "k") >= kalpha)
  # Curves
  data_and_sim_curves_l <- lapply(curve_sets, function(x) { data_and_sim_curves(x) })
  # Curves from which to calculate the convex hull
  curves_for_envelope_l <- lapply(data_and_sim_curves_l, function(x) { x[curves_for_envelope_ind,] })
  # Bounding curves
  LB <- lapply(curves_for_envelope_l, FUN = function(x) { apply(x, MARGIN=2, FUN=min) })
  UB <- lapply(curves_for_envelope_l, FUN = function(x) { apply(x, MARGIN=2, FUN=max) })
  # Update the bounding curves (lo, hi) and kalpha to the first level central regions
  for(i in 1:length(curve_sets)) {
    res_ls[[i]]$lo <- LB[[i]]
    res_ls[[i]]$hi <- UB[[i]]
    attr(res_ls[[i]], "k_alpha") <- NULL
  }

  # Return
  res <- list(global_envelope_ls = res_ls,
              step2_erl = res_erl, step2_erl_curve_set = curve_set_u)
  class(res) <- "combined_global_envelope"
  res
}

#' Combined central region
#'
#' Given a list of 'curve_set' objects (see \code{\link{create_curve_set}}),
#' a combined global central region, i.e. a combined global envelope, is constructed with
#' the functions saved in the curve set objects.
#' Details of the combined central region can be found in Myllymäki and Mrkvička (2018).
#'
#' @param curve_sets A list of objects of type 'curve_set' or \code{\link[spatstat]{envelope}}.
#' @inheritParams central_region
#' @param ... Additional parameters to be passed to \code{\link{central_region}}.
#' @references
#' Myllymäki, M. and Mrkvička, T. (2018). GET: Global envelopes in R.
#' @export
#' @aliases combined_global_envelope
combined_central_region <- function(curve_sets, coverage = 0.50, ...) {
  combined_CR_or_GET(curve_sets, CR_or_GET = "CR", coverage=coverage, ...)
}

#' Combined global envelope test
#'
#' Given a list of 'curve_set' objects (see \code{\link{create_curve_set}}),
#' a combined global envelope test is done with the functions saved in the curve set objects.
#' Details of the combined central region can be found in Myllymäki and Mrkvička (2018).
#'
#' @param curve_sets A list of objects of type 'curve_set' or \code{\link[spatstat]{envelope}}.
#' @inheritParams global_envelope_test
#' @param ... Additional parameters to be passed to \code{\link{central_region}}.
#' @references
#' Myllymäki, M. and Mrkvička, T. (2018). GET: Global envelopes in R.
#' @export
#' @examples
#' # As an example test CSR of the saplings point pattern from spatstat by means of
#' # L, F, G and J functions.
#' require(spatstat)
#' data(saplings)
#' X <- saplings
#'
#' nsim <- 499 # the number of simulations for the tests
#' # Specify distances for different test functions
#' n <- 500 # the number of r-values
#' rmin <- 0; rmax <- 20; rstep <- (rmax-rmin)/n
#' rminJ <- 0; rmaxJ <- 8; rstepJ <- (rmaxJ-rminJ)/n
#' r <- seq(0, rmax, by=rstep)    # r-distances for Lest
#' rJ <- seq(0, rmaxJ, by=rstepJ) # r-distances for Fest, Gest, Jest
#'
#' # Perform simulations of CSR and calculate the L-functions
#' system.time( env_L <- envelope(X, nsim=nsim,
#'  simulate=expression(runifpoint(X$n, win=X$window)),
#'  fun="Lest", correction="translate",
#'  transform = expression(.-r), # Take the L(r)-r function instead of L(r)
#'  r=r,                         # Specify the distance vector
#'  savefuns=TRUE,               # Save the estimated functions
#'  savepatterns=TRUE) )         # Save the simulated patterns
#' # Take the simulations from the returned object
#' simulations <- attr(env_L, "simpatterns")
#' # Then calculate the other test functions F, G, J for each simulated pattern
#' system.time( env_F <- envelope(X, nsim=nsim,
#'                                simulate=simulations,
#'                                fun="Fest", correction="Kaplan", r=rJ,
#'                                savefuns=TRUE) )
#' system.time( env_G <- envelope(X, nsim=nsim,
#'                                simulate=simulations,
#'                                fun="Gest", correction="km", r=rJ,
#'                                savefuns=TRUE) )
#' system.time( env_J <- envelope(X, nsim=nsim,
#'                                simulate=simulations,
#'                                fun="Jest", correction="none", r=rJ,
#'                                savefuns=TRUE) )
#'
#' # Crop the curves to the desired r-interval I
#' curve_set_L <- crop_curves(env_L, r_min=rmin, r_max=rmax)
#' curve_set_F <- crop_curves(env_F, r_min=rminJ, r_max=rmaxJ)
#' curve_set_G <- crop_curves(env_G, r_min=rminJ, r_max=rmaxJ)
#' curve_set_J <- crop_curves(env_J, r_min=rminJ, r_max=rmaxJ)
#'
#' res <- combined_central_region(curve_sets=list(curve_set_L, curve_set_F,
#'                                                curve_set_G, curve_set_J),
#'                                type = "qdir")
#' res <- combined_global_envelope_test(curve_sets=list(curve_set_L, curve_set_F,
#'                                                curve_set_G, curve_set_J),
#'                                     type = "qdir")
#' plot(res, plot_style="ggplot2",
#'      labels=c("L(r)-r", "F(r)", "G(r)", "J(r)"),
#'      separate_yaxes=TRUE, base_size=12)
#'
combined_global_envelope_test <- function(curve_sets, alpha = 0.05, ...) {
  combined_CR_or_GET(curve_sets, CR_or_GET = "GET", coverage=1-alpha, ...)
}

#' Print method for the class 'combined_global_envelope'
#' @usage \method{print}{combined_global_envelope}(x, ...)
#'
#' @param x an 'combined_global_envelope' object
#' @param ... Ignored.
#'
#' @method print combined_global_envelope
#' @export
print.combined_global_envelope <- function(x, ...) {
  print(x$step2_erl)
}

#' Plot method for the class 'combined_global_envelope'
#' @usage \method{plot}{combined_global_envelope}(x, plot_style="basic", level = 1, max_ncols_of_plots=2, ...)
#'
#' @param x an 'combined_global_envelope' object
#' @inheritParams plot.global_envelope
#' @inheritParams env_basic_plot
#' @param level 1 for plotting the combined global envelopes or
#' 2 for plotting the second level ERL test result.
#' @param ... Additional parameters to be passed to \code{\link{plot.global_envelope}}.
#'
#' @method plot combined_global_envelope
#' @export
plot.combined_global_envelope <- function(x, plot_style="basic", level = 1,
                                          max_ncols_of_plots=2, ...) {
  plot_style <- spatstat::pickoption("ptype", plot_style, c(basic = "basic",
                                                            b = "basic",
                                                            fv = "fv",
                                                            f = "fv",
                                                            ggplot2 = "ggplot2",
                                                            ggplot = "ggplot2",
                                                            g = "ggplot2"))
  if(!(level %in% c(1,2))) stop("Unreasonable value for level.\n")
  if(level == 1) {
    if(plot_style %in% c("basic", "ggplot2")) {
      # Create a combined envelope object for plotting
      if(!is.null(x$global_envelope_ls[[1]]$obs))
        res_tmp <- structure(data.frame(r = do.call(c, lapply(x$global_envelope_ls, FUN = function(x) x$r), quote=FALSE),
                                        obs = do.call(c, lapply(x$global_envelope_ls, FUN = function(x) x$obs), quote=FALSE),
                                        central = do.call(c, lapply(x$global_envelope_ls, FUN = function(x) x$central), quote=FALSE),
                                        lo = do.call(c, lapply(x$global_envelope_ls, FUN = function(x) x$lo), quote=FALSE),
                                        hi = do.call(c, lapply(x$global_envelope_ls, FUN = function(x) x$hi), quote=FALSE)),
                             class = class(x$global_envelope_ls[[1]]))
      else
        res_tmp <- structure(data.frame(r = do.call(c, lapply(x$global_envelope_ls, FUN = function(x) x$r), quote=FALSE),
                                        central = do.call(c, lapply(x$global_envelope_ls, FUN = function(x) x$central), quote=FALSE),
                                        lo = do.call(c, lapply(x$global_envelope_ls, FUN = function(x) x$lo), quote=FALSE),
                                        hi = do.call(c, lapply(x$global_envelope_ls, FUN = function(x) x$hi), quote=FALSE)),
                             class = class(x$global_envelope_ls[[1]]))
      attr(res_tmp, "method") <- attr(x$global_envelope_ls[[1]], "method")
      attr(res_tmp, "alternative") <- attr(x$global_envelope_ls[[1]], "alternative")
      if(!is.null(attr(x$step2_erl, "p"))) attr(res_tmp, "p") <- attr(x$step2_erl, "p")
      plot(res_tmp, plot_style=plot_style, ...)
    }
    else {
      n_of_plots <- length(x$global_envelope_ls)
      ncols_of_plots <- min(n_of_plots, max_ncols_of_plots)
      nrows_of_plots <- ceiling(n_of_plots / ncols_of_plots)
      par(mfrow=c(nrows_of_plots, ncols_of_plots))
      for(i in 1:length(x$global_envelope_ls))
        plot(x$global_envelope_ls[[i]], plot_style=plot_style, ...)
    }
  }
  else {
    plot(x$step2_erl, ...)
  }
}

# Quantiles (qdir) or sds (st) for the combined scaled MAD envelope tests
#
# Quantiles (qdir) or sds (st) for the combined scaled MAD envelope tests, i.e.
# the quantities are calculated for all the curve sets provided in the list "curve_sets"
# and the results are returned in corresponding lists.
# @inheritParams combined_scaled_MAD_envelope
#' @importFrom stats quantile
#' @importFrom stats sd
combined_scaled_MAD_bounding_curves_chars <- function(curve_sets, test = c("qdir", "st"), probs = c(0.025, 0.975)) {
  curve_sets_res <- lapply(curve_sets, FUN = function(x) residual(x, use_theo = TRUE))

  switch(test,
         qdir = {
           quant_m_ls <- lapply(curve_sets_res, FUN = curve_set_quant, probs = probs)
           lower_f <- lapply(quant_m_ls, FUN = function(x) as.vector(abs(x[1,])))
           upper_f <- lapply(quant_m_ls, FUN = function(x) as.vector(abs(x[2,])))
         },
         st = {
           lower_f <- upper_f <- lapply(curve_sets_res, FUN = function(x) { as.vector(curve_set_sd(x)) })
         })

  list(lower_f = lower_f, upper_f = upper_f)
}

# Calculate the lower and upper bounding curves of a combined global scaled MAD envelope test
#
# @param central_curves_ls A list containing the central functions for different test functions.
# @param max_u The k_alpha'th largest value of the u_i, i=1,...,nsim+1 for each individual test.
# @param lower_f The first component in the object returned by \code{\link{combined_scaled_MAD_bounding_curves_chars}}.
# @param upper_f The second component in the object returned by \code{\link{combined_scaled_MAD_bounding_curves_chars}}.
combined_scaled_MAD_bounding_curves <- function(central_curves_ls, max_u, lower_f, upper_f) {
  ntests <- length(central_curves_ls)
  if(length(max_u) != ntests | length(lower_f) != ntests | length(upper_f) != ntests)
    stop("The lengths of different arguments do not match.\n")
  # The global 100(1-alpha)% envelope
  # Find the k_alpha'th largest value of the u_i, i=1,...,nsim+1 for each individual test
  # Typically max_u <- res_rank$hi, where res_rank is the combined rank envelope test done.
  # Lower and upper envelope
  lo <- function(i) { as.vector(central_curves_ls[[i]] - max_u[i]*lower_f[[i]]) }
  up <- function(i) { as.vector(central_curves_ls[[i]] + max_u[i]*upper_f[[i]]) }
  list(lower_ls = lapply(1:ntests, FUN = lo), upper_ls = lapply(1:ntests, FUN = up))
}

#' Combined global scaled maximum absolute difference (MAD) envelope tests
#'
#' Given a list of 'curve_set' objects (see \code{\link{create_curve_set}}), a combined global scaled (directional quantile
#' or studentized) MAD envelope test is performed with the test functions saved in the curve set objects.
#' Details of the combined test can be found in Mrkvicka et al. (2017)
#'
#' @param curve_sets A list of objects of type 'curve_set' or \code{\link[spatstat]{envelope}}.
#' @param test Either "qdir" for the \code{\link{qdir_envelope}} test or
#' "st" for the \code{\link{st_envelope}} test.
#' @param alpha The significance level. The 100(1-alpha)\% global envelope will be calculated.
#' @param probs A two-element vector containing the lower and upper
#' quantiles for the envelope, in that order and on the interval [0, 1].
#' The default values are 0.025 and 0.975.
#' @param ... Additional parameters to be passed to \code{\link{qdir_envelope}} (if test = "qdir")
#' or \code{\link{st_envelope}} (if test = "st").
#' @references
#' Mrkvicka, T., Myllymäki, M. and Hahn, U. (2017) Multiple Monte Carlo testing, with applications in spatial point processes.
#' Statistics & Computing 27(5): 1239–1255. DOI: 10.1007/s11222-016-9683-9
#' @export
#' @examples
#' # As an example test CSR of the saplings point pattern from spatstat by means of
#' # L, F, G and J functions.
#' require(spatstat)
#' data(saplings)
#' X <- saplings
#'
#' nsim <- 499 # the number of simulations for the tests
#' # Specify distances for different test functions
#' n <- 500 # the number of r-values
#' rmin <- 0; rmax <- 20; rstep <- (rmax-rmin)/n
#' rminJ <- 0; rmaxJ <- 8; rstepJ <- (rmaxJ-rminJ)/n
#' r <- seq(0, rmax, by=rstep)    # r-distances for Lest
#' rJ <- seq(0, rmaxJ, by=rstepJ) # r-distances for Fest, Gest, Jest
#'
#' # Perform simulations of CSR and calculate the L-functions
#' system.time( env_L <- envelope(X, nsim=nsim,
#'  simulate=expression(runifpoint(X$n, win=X$window)),
#'  fun="Lest", correction="translate",
#'  transform = expression(.-r), # Take the L(r)-r function instead of L(r)
#'  r=r,                         # Specify the distance vector
#'  savefuns=TRUE,               # Save the estimated functions
#'  savepatterns=TRUE) )         # Save the simulated patterns
#' # Take the simulations from the returned object
#' simulations <- attr(env_L, "simpatterns")
#' # Then calculate the other test functions F, G, J for each simulated pattern
#' system.time( env_F <- envelope(X, nsim=nsim,
#'                                simulate=simulations,
#'                                fun="Fest", correction="Kaplan", r=rJ,
#'                                savefuns=TRUE) )
#' system.time( env_G <- envelope(X, nsim=nsim,
#'                                simulate=simulations,
#'                                fun="Gest", correction="km", r=rJ,
#'                                savefuns=TRUE) )
#' system.time( env_J <- envelope(X, nsim=nsim,
#'                                simulate=simulations,
#'                                fun="Jest", correction="none", r=rJ,
#'                                savefuns=TRUE) )
#'
#' # Crop the curves to the desired r-interval I
#' curve_set_L <- crop_curves(env_L, r_min=rmin, r_max=rmax)
#' curve_set_F <- crop_curves(env_F, r_min=rminJ, r_max=rmaxJ)
#' curve_set_G <- crop_curves(env_G, r_min=rminJ, r_max=rmaxJ)
#' curve_set_J <- crop_curves(env_J, r_min=rminJ, r_max=rmaxJ)
#'
#' # The combined directional quantile envelope test
#' res <- combined_scaled_MAD_envelope(curve_sets=list(curve_set_L, curve_set_F,
#'                                                     curve_set_G, curve_set_J),
#'                                     test = "qdir")
#' plot(res, plot_style="ggplot2",
#'      labels=c("L(r)-r", "F(r)", "G(r)", "J(r)"),
#'      separate_yaxes=TRUE, base_size=12)
#'
combined_scaled_MAD_envelope <- function(curve_sets, test = c("qdir", "st"), alpha = 0.05, probs = c(0.025, 0.975), ...) {

    ntests <- length(curve_sets)
    test <- match.arg(test)
    curve_sets <- check_curve_set_dimensions(curve_sets)
    # Make the individual tests saving the deviations
    switch(test, 
            qdir = {
                res_ls <- lapply(curve_sets, FUN = function(x) { qdir_envelope(x, alpha=alpha, savedevs=TRUE, probs = probs, ...) })
                method <- "Combined studentised envelope test"
            },
            st = {
                res_ls <- lapply(curve_sets, FUN = function(x) { st_envelope(x, alpha=alpha, savedevs=TRUE, ...) })
                method <- "Combined directional quantile envelope test"
            })
    # Calculate quantiles (qdir) or sds (st)
    envchars <- combined_scaled_MAD_bounding_curves_chars(curve_sets, test=test, probs=probs)

    # Create a curve_set for the rank test
    u_ls <- lapply(res_ls, FUN = function(x) attr(x, "k"))
    u_mat <- do.call(rbind, u_ls, quote=FALSE)
    curve_set_u <- create_curve_set(list(r=1:ntests, obs=u_mat[,1], sim_m=u_mat[,-1], is_residual=FALSE))
    # Perform the one-sided (greater is significant) rank test
    res_rank <- rank_envelope(curve_set_u, alpha=alpha, savedevs=TRUE, alternative="greater", type="rank", ties="erl")

    central_curves_ls <- lapply(res_ls, function(x) x$central)
    bounding_curves <- combined_scaled_MAD_bounding_curves(central_curves_ls=central_curves_ls, max_u=res_rank$hi,
                                                           lower_f=envchars$lower_f, upper_f=envchars$upper_f)

    # Create a combined envelope object for plotting
    res_env <- structure(data.frame(r = do.call(c, lapply(res_ls, FUN = function(x) x$r), quote=FALSE),
                              obs = do.call(c, lapply(res_ls, FUN = function(x) x$obs), quote=FALSE),
                              central = do.call(c, lapply(res_ls, FUN = function(x) x$central), quote=FALSE),
                              lo = do.call(c, bounding_curves$lower_ls, quote=FALSE),
                              hi = do.call(c, bounding_curves$upper_ls, quote=FALSE)),
                         class = c("combined_scaled_MAD_envelope", "envelope_test"))
    attr(res_env, "method") <- method
    attr(res_env, "alternative") <- attr(res_ls[[1]], "alternative")
    attr(res_env, "p") <- attr(res_rank, "p")
    attr(res_env, "p_interval") <- attr(res_rank, "p_interval")
    # return
    res <- structure(list(rank_test = res_rank, envelope = res_env),
                     class = "combined_scaled_MAD_test")
    attr(res, "rank_test_curve_set") <- curve_set_u
    res
}

#' Print method for the class 'combined_scaled_MAD_test'
#' @usage \method{print}{combined_scaled_MAD_test}(x, ...)
#'
#' @param x An 'combined_scaled_MAD_test' object
#' @param ... Ignored.
#'
#' @method print combined_scaled_MAD_test
#' @export
print.combined_scaled_MAD_test <- function(x, ...) {
    print(x$rank_test)
}

#' Plot method for the class 'combined_scaled_MAD_test'
#'
#' Plot method for the class 'combined_scaled_MAD_test'.
#' If plot_type is "rank", then the output of the combined rank envelope test
#' is plotted. If plot_type is "envelope", then the combined scaled MAD envelope
#' is plotted.
#' @usage \method{plot}{combined_scaled_MAD_test}(x, plot_type = c("envelope", "rank"), ...)
#'
#' @param x An 'combined_scaled_MAD_test' object
#' @param plot_type "rank" or "qdir_envelopes".
#' @param ... Arguments passed to \code{\link{plot.envelope_test}}.
#'
#' @method plot combined_scaled_MAD_test
#' @export
plot.combined_scaled_MAD_test <- function(x, plot_type = c("envelope", "rank"), ...) {
    plot_type <- match.arg(plot_type)
    switch(plot_type,
           rank = {
               plot.envelope_test(x[['rank_test']], ...)
           },
           envelope = {
               plot.envelope_test(x[['envelope']], ...)
           })
}
