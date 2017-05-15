#' Quantiles (qdir) or sds (st) for the combined scaled MAD envelope tests
#'
#' Quantiles (qdir) or sds (st) for the combined scaled MAD envelope tests, i.e.
#' the quantities are calculated for all the curve sets provided in the list "curve_sets"
#' and the results are returned in corresponding lists.
#' @inheritParams combined_scaled_MAD_envelope
#' @importFrom stats quantile
#' @importFrom stats sd
combined_scaled_MAD_bounding_curves_chars <- function(curve_sets, test = c("qdir", "st"), probs = c(0.025, 0.975)) {
    curve_sets_res <- lapply(curve_sets, FUN = function(x) residual(x, use_theo = TRUE))

    switch(test,
            qdir = {
                quant_m_ls <- lapply(curve_sets_res, FUN = function(x) apply(x[['sim_m']], MARGIN=1, FUN=stats::quantile, probs = probs))
                lower_f <- lapply(quant_m_ls, FUN = function(x) as.vector(abs(x[1,])))
                upper_f <- lapply(quant_m_ls, FUN = function(x) as.vector(abs(x[2,])))
            },
            st = {
                lower_f <- upper_f <- lapply(curve_sets_res, FUN = function(x) { as.vector(apply(x[['sim_m']], MARGIN=1, FUN=stats::sd)) })
            })

    list(lower_f = lower_f, upper_f = upper_f)
}

#' Calculate the lower and upper bounding curves of a combined global scaled MAD envelope test
#'
#' @param central_curves_ls A list containing the central functions for different test functions.
#' @param max_u The k_alpha'th largest value of the u_i, i=1,...,nsim+1 for each individual test.
#' @param lower_f The first component in the object returned by \code{\link{combined_scaled_MAD_bounding_curves_chars}}.
#' @param upper_f The second component in the object returned by \code{\link{combined_scaled_MAD_bounding_curves_chars}}.
combined_scaled_MAD_bounding_curves <- function(central_curves_ls, max_u, lower_f, upper_f) {
    ntests <- length(central_curves_ls)
    if(length(max_u) != ntests | length(lower_f) != ntests | length(upper_f) != ntests)
        stop("The lengths of different arguments do not match.\n")
    # The global 100(1-alpha)% envelope
    # Find the k_alpha'th largest value of the u_i, i=1,...,nsim+1 for each individual test
    # Typically max_u <- res_rank$upper, where res_rank is the combined rank envelope test done.
    # Lower and upper envelope
    lo <- function(i) { as.vector(central_curves_ls[[i]] - max_u[i]*lower_f[[i]]) }
    up <- function(i) { as.vector(central_curves_ls[[i]] + max_u[i]*upper_f[[i]]) }
    list(lower_ls = lapply(1:ntests, FUN = lo), upper_ls = lapply(1:ntests, FUN = up))
}


#' Combined global scaled maximum absolute difference (MAD) envelope tests
#'
#' Given a list of 'curve_set' objects (see \code{\link{create_curve_set}}), a combined global scaled (directional quantile
#' or studentized) MAD envelope test is performed with the test functions saved in the curve set objects.
#' Details of the combined test can be found in Mrkvicka et al.
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
#' Mrkvicka, T., Myllymäki, M. and Hahn, U. Multiple Monte Carlo testing, with applications in spatial point processes.
#' Revision submitted to Statistics & Computing.
#' @export
combined_scaled_MAD_envelope <- function(curve_sets, test = c("qdir", "st"), alpha = 0.05, probs = c(0.025, 0.975), ...) {

    ntests <- length(curve_sets)
    test <- match.arg(test)
    curve_sets <- check_curve_set_dimensions(curve_sets)
    # Make the individual tests saving the deviations
    switch(test, 
            qdir = {
                res_ls <- lapply(curve_sets, FUN = function(x) { qdir_envelope(x, alpha=alpha, savedevs=TRUE, probs = probs, ...) })
            },
            st = {
                res_ls <- lapply(curve_sets, FUN = function(x) { st_envelope(x, alpha=alpha, savedevs=TRUE, ...) })
            })
    # Calculate quantiles (qdir) or sds (st)
    envchars <- combined_scaled_MAD_bounding_curves_chars(curve_sets, test=test, probs=probs)

    # Create a curve_set for the rank test
    u_ls <- lapply(res_ls, FUN = function(x) x$u)
    u_mat <- do.call(cbind, u_ls, quote=FALSE)
    u_mat <- rbind(res_ls[[1]]$u, res_ls[[2]]$u, res_ls[[3]]$u, res_ls[[4]]$u)
    curve_set_u <- create_curve_set(list(r=1:ntests, obs=u_mat[,1], sim_m=u_mat[,-1], is_residual=FALSE))
    # Perform the one-sided (greater is significant) rank test
    res_rank <- rank_envelope(curve_set_u, alpha=alpha, savedevs=TRUE, alternative="greater", lexo=TRUE)

    central_curves_ls <- lapply(res_ls, function(x) x$central)
    bounding_curves <- combined_scaled_MAD_bounding_curves(central_curves_ls=central_curves_ls, max_u=res_rank$upper,
                                                           lower_f=envchars$lower_f, upper_f=envchars$upper_f)

    # Create a combined envelope object for plotting
    res_env <- structure(list(r = do.call(c, lapply(res_ls, FUN = function(x) x$r), quote=FALSE),
                              method = res_ls[[1]]$method,
                              alternative = res_ls[[1]]$alternative,
                              p = res_rank$p,
                              p_interval = res_rank$p_interval,
                              central = do.call(c, lapply(res_ls, FUN = function(x) x$central), quote=FALSE),
                              obs = do.call(c, lapply(res_ls, FUN = function(x) x$obs), quote=FALSE),
                              lo = do.call(c, bounding_curves$lower_ls, quote=FALSE),
                              hi = do.call(c, bounding_curves$upper_ls, quote=FALSE)),
                         class = c("combined_scaled_MAD_envelope", "envelope_test"))

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
