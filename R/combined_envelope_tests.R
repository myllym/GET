# Quantiles (qdir) or sds (st) for the combined scaled MAD envelope tests
#
# Quantiles (qdir) or sds (st) for the combined scaled MAD envelope tests, i.e.
# the quantities are calculated for all the curve sets provided in the list "curve_sets"
# and the results are returned in corresponding lists.
# @inheritParams combined_scaled_MAD_envelope_test
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
# @param max_u The u_(alpha)'th largest value of the u_i, i=1,...,nsim+1 for each individual test.
# @param lower_f The first component in the object returned by \code{\link{combined_scaled_MAD_bounding_curves_chars}}.
# @param upper_f The second component in the object returned by \code{\link{combined_scaled_MAD_bounding_curves_chars}}.
combined_scaled_MAD_bounding_curves <- function(central_curves_ls, max_u, lower_f, upper_f) {
  ntests <- length(central_curves_ls)
  if(length(max_u) != ntests | length(lower_f) != ntests | length(upper_f) != ntests)
    stop("The lengths of different arguments do not match.")
  # The global 100(1-alpha)% envelope
  # Find the u_(alpha)'th largest value of the u_i, i=1,...,nsim+1 for each individual test
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
#' Details of this combined test can be found in Mrkvicka et al. (2017).
#' The implementation of this test is provided here for historical reasons:
#' we recommend now instead the use of \code{\link{global_envelope_test}} also for combined tests;
#' these combined tests are there implemented as described in Myllymäki and Mrkvička (2020).
#'
#' @inheritParams global_envelope_test
#' @param type Either "qdir" for the direction quantile envelope test or
#' "st" for the studentized envelope test.
#' @param ... Additional parameters to be passed to \code{\link{central_region}}.
#' @references
#' Mrkvička, T., Myllymäki, M. and Hahn, U. (2017) Multiple Monte Carlo testing, with applications in spatial point processes.
#' Statistics & Computing 27(5): 1239–1255. DOI: 10.1007/s11222-016-9683-9
#'
#' Myllymäki, M. and Mrkvička, T. (2020). GET: Global envelopes in R. arXiv:1911.06583 [stat.ME]
#' @export
#' @examples
#' if(require("spatstat", quietly=TRUE)) {
#'   # As an example test CSR of the saplings point pattern from spatstat by means of
#'   # L, F, G and J functions.
#'   data("saplings")
#'   X <- as.ppp(saplings, W=square(75))
#'
#'   \donttest{nsim <- 499 # Number of simulations for the tests}
#'   \dontshow{nsim <- 19 # Number of simulations for testing}
#'   # Specify distances for different test functions
#'   n <- 500 # the number of r-values
#'   rmin <- 0; rmax <- 20; rstep <- (rmax-rmin)/n
#'   rminJ <- 0; rmaxJ <- 8; rstepJ <- (rmaxJ-rminJ)/n
#'   r <- seq(0, rmax, by=rstep)    # r-distances for Lest
#'   rJ <- seq(0, rmaxJ, by=rstepJ) # r-distances for Fest, Gest, Jest
#'
#'   # Perform simulations of CSR and calculate the L-functions
#'   env_L <- envelope(X, nsim=nsim,
#'    simulate=expression(runifpoint(X$n, win=X$window)),
#'    fun="Lest", correction="translate",
#'    transform = expression(.-r), # Take the L(r)-r function instead of L(r)
#'    r=r,                         # Specify the distance vector
#'    savefuns=TRUE,               # Save the estimated functions
#'    savepatterns=TRUE)           # Save the simulated patterns
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
#'                    savefuns=TRUE)
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
#'   # The combined directional quantile envelope test
#'   res <- combined_scaled_MAD_envelope_test(
#'              curve_sets=list(L=curve_set_L, F=curve_set_F,
#'                              G=curve_set_G, J=curve_set_J),
#'              type = "qdir")
#'   plot(res)
#' }
#'
combined_scaled_MAD_envelope_test <- function(curve_sets, type = c("qdir", "st"), alpha = 0.05,
                                         probs = c(0.025, 0.975), central = "mean", ...) {
    ntests <- length(curve_sets)
    if(ntests <= 1) stop("Number of functions should be at least two.")
    type <- match.arg(type)
    if(!all(sapply(curve_sets, FUN = function(x) { inherits(x, "envelope") }))) {
      tmp <- lapply(curve_sets, FUN = convert_to_curveset)
      if(!all(sapply(tmp, curve_set_is1obs)))
        stop("The curve_set does not contain one observed function. Testing does not make sense.\n Did you want to construct a central region of your data? See the function central_region.")
    }
    curve_sets <- check_curve_set_dimensions(curve_sets)
    # Make the individual tests saving the deviations
    switch(type,
            qdir = {
                res_ls <- lapply(curve_sets,
                                 FUN = function(x) {
                                   central_region(x, type="qdir", coverage=1-alpha,
                                                  probs = probs, central=central, ...)
                                 })
            },
            st = {
                res_ls <- lapply(curve_sets,
                                 FUN = function(x) {
                                   central_region(x, type="st", coverage=1-alpha,
                                                  central=central, ...)
                                 })
            })
    # Calculate quantiles (qdir) or sds (st)
    envchars <- combined_scaled_MAD_bounding_curves_chars(curve_sets, type=type, probs=probs)

    # Create a curve_set for the rank test
    u_ls <- lapply(res_ls, FUN = function(x) attr(x, "M"))
    u_mat <- do.call(rbind, u_ls, quote=FALSE)
    curve_set_u <- create_curve_set(list(r=1:ntests, obs=u_mat[,1], sim_m=u_mat[,-1]))
    # Perform the one-sided (greater is significant) rank test
    res_rank <- global_envelope_test(curve_set_u, type="rank", alpha=alpha, alternative="greater", ties="erl")

    central_curves_ls <- lapply(res_ls, function(x) x$central)
    bounding_curves <- combined_scaled_MAD_bounding_curves(central_curves_ls=central_curves_ls, max_u=res_rank$hi,
                                                           lower_f=envchars$lower_f, upper_f=envchars$upper_f)
    # Update the bounding curves (lo, hi) and Malpha to the first level central regions
    for(i in 1:length(curve_sets)) {
      res_ls[[i]]$lo <- bounding_curves$lower_ls[[i]]
      res_ls[[i]]$hi <- bounding_curves$upper_ls[[i]]
      attr(res_ls[[i]], "M_alpha") <- NULL
    }

    # Return
    attr(res_ls, "level2_ge") <- res_rank
    attr(res_ls, "level2_curve_set") <- curve_set_u
    attr(res_ls, "method") <- "Combined global test (scaled MAD)"
    class(res_ls) <- c("combined_global_envelope", class(res_ls))
    if(!is.null(attr(res_ls[[1]], "argu")))
      res_ls <- envelope_set_labs(res_ls, xlab=attr(res_ls[[1]], "xlab"),
                                  ylab=substitute(italic(T(i)), list(i=attr(res_ls[[1]], "argu"))))
    else
      res_ls <- envelope_set_labs(res_ls, xlab=expression(italic(r)),
                                  ylab=expression(italic(T(r))))
    attr(res_ls, "alternative") <- get_alternative(res_ls[[1]])
    attr(res_ls, "type") <- type
    attr(res_ls, "alpha") <- alpha
    attr(res_ls, "p") <- attr(res_rank, "p")
    attr(res_ls, "p_interval") <- attr(res_rank, "p_interval")
    attr(res_ls, "nstep") <- 2
    res_ls
}
