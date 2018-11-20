# Quantiles (qdir) or sds (st) for the combined scaled MAD envelope tests
#
# Quantiles (qdir) or sds (st) for the combined scaled MAD envelope tests, i.e.
# the quantities are calculated for all the curve sets provided in the list "curve_sets"
# and the results are returned in corresponding lists.
# @inheritParams combined_scaled_MAD_envelope
#' @importFrom stats quantile
#' @importFrom stats sd
combined_scaled_MAD_bounding_curves_chars <- function(curve_sets, type = c("qdir", "st"), probs = c(0.025, 0.975)) {
  curve_sets_res <- lapply(curve_sets, FUN = function(x) residual(x, use_theo = TRUE))

  switch(type,
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
#' @param type Either "qdir" for the \code{\link{qdir_envelope}} test or
#' "st" for the \code{\link{st_envelope}} test.
#' @param alpha The significance level. The 100(1-alpha)\% global envelope will be calculated.
#' @param probs A two-element vector containing the lower and upper
#' quantiles for the envelope, in that order and on the interval [0, 1].
#' The default values are 0.025 and 0.975.
#' @inheritParams central_region
#' @param ... Additional parameters to be passed to \code{\link{qdir_envelope}} (if type = "qdir")
#' or \code{\link{st_envelope}} (if type = "st").
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
#'                                     type = "qdir")
#' plot(res, plot_style="ggplot2",
#'      labels=c("L(r)-r", "F(r)", "G(r)", "J(r)"),
#'      separate_yaxes=TRUE, base_size=12)
#'
combined_scaled_MAD_envelope <- function(curve_sets, type = c("qdir", "st"), alpha = 0.05,
                                         probs = c(0.025, 0.975), central = "mean", ...) {
    ntests <- length(curve_sets)
    if(ntests <= 1) stop("Number of functions should be at least two.\n")
    type <- match.arg(type)
    curve_sets <- check_curve_set_dimensions(curve_sets)
    # Make the individual tests saving the deviations
    switch(type,
            qdir = {
                res_ls <- lapply(curve_sets, FUN = function(x) { central_region(x, type="qdir", coverage=1-alpha, probs = probs, central=central, ...) })
            },
            st = {
                res_ls <- lapply(curve_sets, FUN = function(x) { central_region(x, type="st", coverage=1-alpha, central=central, ...) })
            })
    # Calculate quantiles (qdir) or sds (st)
    envchars <- combined_scaled_MAD_bounding_curves_chars(curve_sets, type=type, probs=probs)

    # Create a curve_set for the rank test
    u_ls <- lapply(res_ls, FUN = function(x) attr(x, "k"))
    u_mat <- do.call(rbind, u_ls, quote=FALSE)
    curve_set_u <- create_curve_set(list(r=1:ntests, obs=u_mat[,1], sim_m=u_mat[,-1], is_residual=FALSE))
    # Perform the one-sided (greater is significant) rank test
    res_rank <- rank_envelope(curve_set_u, alpha=alpha, alternative="greater", type="rank", ties="erl")

    central_curves_ls <- lapply(res_ls, function(x) x$central)
    bounding_curves <- combined_scaled_MAD_bounding_curves(central_curves_ls=central_curves_ls, max_u=res_rank$hi,
                                                           lower_f=envchars$lower_f, upper_f=envchars$upper_f)
    # Update the bounding curves (lo, hi) and kalpha to the first level central regions
    for(i in 1:length(curve_sets)) {
      res_ls[[i]]$lo <- bounding_curves$lower_ls[[i]]
      res_ls[[i]]$hi <- bounding_curves$upper_ls[[i]]
      attr(res_ls[[i]], "k_alpha") <- NULL
    }

    # Return
    res <- list(global_envelope_ls = res_ls,
                step2_test = res_rank, step2_test_curve_set = curve_set_u)
    class(res) <- "combined_global_envelope"
    res
}
