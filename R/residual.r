#' Subtract S_{H_0} from the summary functions.
#' @inheritParams convert_envelope
#' @param use_theo Whether to use the theoretical summary function or the
#'   mean of the simulations.
#' @return A curve set object containing residual summary functions. theo is
#'   no longer included.
#' @export
residual <- function(curve_set, use_theo = TRUE) {
    curve_set <- convert_envelope(curve_set)

    if(with(curve_set, exists('is_residual')) && curve_set[['is_residual']]) {
        cat("The curve_set object contains already residuals T_i(r) - T_0(r), \n",
            "use_theo ignored.\n")
        res <- curve_set
    }
    else {
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
    }
    res
}

#' Define T_0 from a curve_set object
#'
#' Define T_0 from a curve_set object
#'
#'
#' T_0 is the expectation of the test function under H_0.
get_T_0 <- function(curve_set) {
    curve_set <- convert_envelope(curve_set)

    if(with(curve_set, exists('is_residual'))) {
        if(!curve_set[['is_residual']]) {
            if(with(curve_set, exists('theo'))) {
                T_0 <- curve_set[['theo']]
            }
            else {
                sim_curves <- t(curve_set[['sim_m']])
                T_0 <- colMeans(sim_curves)
            }
        }
        else {
            T_0 <- rep(0, times=length(curve_set$r))
        }
    }
    else { # Assume curve_set does not contain residuals
        if(with(curve_set, exists('theo'))) {
            T_0 <- curve_set[['theo']]
        }
        else {
            sim_curves <- t(curve_set[['sim_m']])
            T_0 <- colMeans(sim_curves)
        }
    }
    T_0
}
