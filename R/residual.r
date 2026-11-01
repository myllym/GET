#' Residual form of the functions
#'
#' Subtract the theoretical function \eqn{S_{H_0}}{S_{H_0}} or the mean (or median) of the functions
#' in the curve set. If the \code{curve_set} object contains already residuals
#' \eqn{T_i(r) - T_0(r)}{T_i(r) - T_0(r)}, \code{use_theo} ignored and the same object returned.
#'
#'
#' The mean of the functions in the \code{curve_set} is the pointwise mean of all functions.
#' The median of the functions in the \code{curve_set} is the pointwise median of all functions.
#' If \code{use_theo = TRUE}, but the component \code{theo} does not exist in the
#' \code{curve_set}, the mean/median of the functions is used silently.
#' @inheritParams crop_curves
#' @param use_theo Whether to use the theoretical summary function or the
#' central of the functions in the curve_set.
#' @param central Either "mean" (default) or "median". If the curve sets do not contain
#' the component \code{theo} for the theoretical central function, or \code{use_theo=FALSE},
#' then the central function is calculated either as the mean or median of functions
#' provided in the curve sets. Default is the mean.
#' @param ... Additional parameters passed to \code{\link[base]{rowMeans}} if
#' \code{central="mean"} or \code{\link[stats]{median}} if \code{central="median"}.
#' @return A curve set object containing residual summary functions. theo is
#'   no longer included.
#' @export
residual <- function(curve_set, use_theo = TRUE, central=c("mean", "median"), ...) {
    curve_set <- convert_envelope(curve_set)
    central <- match.arg(central)
    switch(central,
           mean = { centralf <- curve_set_mean },
           median = { centralf <- curve_set_median })

    if(curve_set_isresidual(curve_set)) {
        # "The curve_set object contains already residuals T_i(r) - T_0(r), use_theo ignored.
        return(curve_set)
    }
    else {
        if(length(use_theo) != 1L || !is.logical(use_theo)) {
            stop('use_theo must be either TRUE or FALSE.')
        }

        theo <- curve_set[['theo']]
        if(is.null(theo) && use_theo) {
            # warnings('use_theo == TRUE but the theoretical curve is missing. ',
            #          'Behaving as if use_theo == FALSE.')
            # Silently setting use_theo to FALSE, when 'theo' not provided:
            use_theo <- FALSE
        }

        if(use_theo) mid <- theo
        else mid <- centralf(curve_set, ...)

        curve_set[['funcs']] <- curve_set[['funcs']] - mid
        curve_set[['is_residual']] <- TRUE
        curve_set[['theo']] <- NULL
        return(curve_set)
    }
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
    if(curve_set_isresidual(curve_set))
        rep(0, times=curve_set_narg(curve_set))
    else if(!is.null(curve_set[['theo']]))
        curve_set[['theo']]
    else
        centralf(curve_set)
}
