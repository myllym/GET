#' Residual form of the functions
#'
#' Subtract the theoretical function S_{H_0} or the mean of the functions
#' in the curve set. If the \code{curve_set} object contains already residuals
#' T_i(r) - T_0(r), use_theo ignored and the same object returned.
#'
#'
#' The mean of the functions in the \code{curve_set} is the mean of all functions.
#' If \code{use_theo = TRUE}, but the component \code{theo} does not exist in the
#' \code{curve_set}, the mean of the functions is used silently.
#' @inheritParams crop_curves
#' @param use_theo Whether to use the theoretical summary function or the
#' mean of the functions in the curve_set.
#' @return A curve set object containing residual summary functions. theo is
#'   no longer included.
#' @export
residual <- function(curve_set, use_theo = TRUE) {
    curve_set <- convert_envelope(curve_set)

    if(with(curve_set, exists('is_residual')) && curve_set[['is_residual']]) {
        # "The curve_set object contains already residuals T_i(r) - T_0(r), use_theo ignored.
        res <- curve_set
    }
    else {
        if(length(use_theo) != 1L || !inherits(use_theo, 'logical') ||
            !is.finite(use_theo)) {
            stop('use_theo must be either TRUE or FALSE.')
        }

        theo <- curve_set[['theo']]
        n_theo <- length(theo)
        if(n_theo < 1L && use_theo) {
            # warnings('use_theo == TRUE but the theoretical curve is missing. ',
            #          'Behaving as if use_theo == FALSE.')
            # Silently setting use_theo to FALSE, when 'theo' not provided:
            use_theo <- FALSE
        }

        if(use_theo) {
            mid <- theo
        } else {
            mid <- curve_set_mean(curve_set)
        }

        res <- with(curve_set, list(r = r,
                                    obs = obs - mid))
        if(!is.null(curve_set[['sim_m']])) res[['sim_m']] <- curve_set[['sim_m']] - mid
        res[['is_residual']] <- TRUE

        res <- create_curve_set(res)
    }
    res
}

# The central function T_0
#
# Calculate the central function T_0 from a curve_set object.
#
# @inheritParams crop_curves
# @param central Whether to take the mean or median of functions as the central function.
get_T_0 <- function(curve_set, central = c("mean", "median")) {
    curve_set <- convert_envelope(curve_set)
    central <- match.arg(central)
    switch(central,
           mean = { centralf <- curve_set_mean },
           median = { centralf <- curve_set_median })
    if(with(curve_set, exists('is_residual'))) {
        if(!curve_set[['is_residual']]) {
            if(with(curve_set, exists('theo'))) {
                T_0 <- curve_set[['theo']]
            }
            else {
                T_0 <- centralf(curve_set)
            }
        }
        else {
            T_0 <- rep(0, times=length(curve_set[['r']]))
        }
    }
    else { # Assume curve_set does not contain residuals
        if(with(curve_set, exists('theo'))) {
            T_0 <- curve_set[['theo']]
        }
        else {
            T_0 <- centralf(curve_set)
        }
    }
    T_0
}
