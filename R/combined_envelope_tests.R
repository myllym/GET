#' Combined global scaled maximum absolute difference (MAD) envelope tests
#'
#' @param curve_sets A list of objects of type 
#' @param test Either "qdir" for the \code{\link{qdir_envelope}} test or
#' "st" for the \code{\link{st_envelope}} test.
#' @param alpha The significance level. The 100(1-alpha)\% global envelope will be calculated.
#' @param ... Parameters to be passed to \code{\link{qdir_envelope}} or \code{\link{st_envelope}}
#' depending on the argument 'test'.
#' @export
combined_scaled_MAD_envelope <- function(curve_sets, test = c("qdir", "st"), alpha = 0.05, ...) {

    ntests <- length(curve_sets)
    test <- match.arg(test)

    switch(test, 
           qdir = {
              res_ls <- lapply(curve_sets, FUN = function(x) { qdir_envelope(x, savedevs=TRUE, save_env_details=TRUE, ...) })
           },
           st = {
              res_ls <- lapply(curve_sets, FUN = function(x) { st_envelope(x, savedevs=TRUE, save_env_details=TRUE, ...) })
           })

    # Create a curve_set for the rank test
    u_ls <- lapply(res_ls, FUN = function(x) x$u)
    u_mat <- do.call(cbind, u_ls, quote=FALSE)
    u_mat <- rbind(res_ls[[1]]$u, res_ls[[2]]$u, res_ls[[3]]$u, res_ls[[4]]$u)
    curve_set_u <- create_curve_set(list(r=1:ntests, obs=u_mat[,1], sim_m=u_mat[,-1], is_residual=FALSE))
    # Perform the one-sided (greater is significant) rank test
    res_rank <- rank_envelope(curve_set_u, alpha=alpha, savedevs=TRUE, alternative="greater", lexo=TRUE)

    # The global 100(1-alpha)% envelope
    # Find the k_alpha'th largest value of the u_i, i=1,...,nsim+1 for each individual test
    max_u <- res_rank$upper
    # Lower and upper envelope
    switch(test,
           qdir = {
               lower <- function(i) { res_ls[[i]]$central_curve - max_u[i]*attr(res_ls[[i]], "lower_quantile") }
               upper <- function(i) { res_ls[[i]]$central_curve + max_u[i]*attr(res_ls[[i]], "upper_quantile") }
           },
           st = {
               lower <- function(i) { res_ls[[i]]$central_curve - max_u[i]*attr(res_ls[[i]], "sd") }
               upper <- function(i) { res_ls[[i]]$central_curve + max_u[i]*attr(res_ls[[i]], "sd") }
           })

    LB_ls <- lapply(1:ntests, FUN = lower)
    UB_ls <- lapply(1:ntests, FUN = upper)

    # Create a combined envelope object for plotting
    res_env <- NULL
    res_env$r <- do.call(c, lapply(res_ls, FUN = function(x) x$r), quote=FALSE)
    res_env$method <- res_ls[[1]]$method
    res_env$alternative <- res_ls[[1]]$alternative
    res_env$p <- res_rank$p
    res_env$p_interval <- res_rank$p_interval
    res_env$central_curve <- do.call(c, lapply(res_ls, FUN = function(x) x$central_curve), quote=FALSE)
    res_env$data_curve <- do.call(c, lapply(res_ls, FUN = function(x) x$data_curve), quote=FALSE)
    res_env$lower <- do.call(c, LB_ls, quote=FALSE)
    res_env$upper <- do.call(c, UB_ls, quote=FALSE)
    class(res_env) <- c("combined_scaled_MAD_envelope", "envelope_test")

    # return
    res <- list(rank_test = res_rank, envelope = res_env)
    class(res) <- "combined_scaled_MAD_test"
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
print.envelope_test <- function(x, ...) {
    print(x$rank_test)
}

#' Plot method for the class 'combined_scaled_MAD_test'
#'
#' Plot method for the class 'combined_scaled_MAD_test'.
#' If plot_type is "rank", then the output of the combined rank envelope test
#' is plotted. If plot_type is "envelope", then the combined scaled MAD envelope
#' is plotted.
#' @usage \method{plot}{combined_scaled_MAD_test}(x, ...)
#'
#' @param x An 'combined_scaled_MAD_test' object
#' @param plot_type "rank" or "qdir_envelopes".
#' @param ... Arguments passed to \code{\link{plot.envelope_test}}.
#'
#' @method plot combined_scaled_MAD_test
#' @export
plot.envelope_test <- function(x, plot_type = c("rank", "envelopes"), ...) {
    plot_type <- match.arg(plot_type)
    switch(plot_type,
           rank = {
               plot(x[['rank_test']], ...)
           },
           envelopes = {
               plot(x[['envelope']], ...)
           })
}
