#' Subtract S_{H_0} from the summary functions.
#' @inheritParams convert_envelope
#' @param use_theo Whether to use the theoretical summary function or the
#'   mean of the simulations.
#' @return A curve set object containing residual summary functions. theo is
#'   no longer included.
residual <- function(curve_set, use_theo = TRUE) {
    curve_set <- convert_envelope(curve_set)

    if (length(use_theo) != 1L || !inherits(use_theo, 'logical') ||
        !is.finite(use_theo)) {
        stop('use_theo must be either TRUE or FALSE.')
    }

    theo <- curve_set[['theo']]
    n_theo <- length(theo)
    if (n_theo < 1L && use_theo) {
        warning('use_theo == TRUE but the theoretical curve is missing. ',
                'Behaving as if use_theo == FALSE.')
        use_theo <- FALSE
    }

    if (use_theo) {
        mid <- theo
    } else {
        mid <- apply(curve_set[['sim_m']], 1, mean)
    }

    res <- with(curve_set, list(r = r,
                                obs = obs - mid,
                                sim_m = sim_m - mid))
    res[['is_residual']] <- TRUE

    res <- create_curve_set(res)
    res
}
