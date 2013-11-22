#' Crop the curves to a certain interval.
#'
#' Crop the curves to a certain interval in preparation for a deviation
#' test.
#'
#' @inheritParams convert_envelope
#' @param r_min The minimum radius to include.
#' @param r_max The maximum radius to include.
#' @return A curve_set object containing the cropped summary functions and
#'   the cropped radius vector.
#' @export
crop_curves <- function(curve_set, r_min = NULL, r_max = NULL) {
    curve_set <- convert_envelope(curve_set)

    n_r_min <- length(r_min)
    if (n_r_min > 0L && (n_r_min != 1L || !is.finite(r_min))) {
        stop('r_min must be a finite scalar value or NULL.')
    }

    n_r_max <- length(r_max)
    if (n_r_max > 0L && (n_r_max != 1L || !is.finite(r_max))) {
        stop('r_max must be a finite scalar value or NULL.')
    }

    r <- curve_set[['r']]

    if (n_r_min == 1L) {
        if (n_r_max == 1L) {
            if (r_min >= r_max) {
                stop('r_min must be smaller than r_max.')
            }
            cut_idx <- which(r >= r_min & r <= r_max)
        } else {
            cut_idx <- which(r >= r_min)
        }
    } else {
        if (n_r_max == 1L) {
            cut_idx <- which(r <= r_max)
        } else {
            return(curve_set)
        }
    }

    if (length(cut_idx) < 1L) {
        stop('r_min and r_max cropped everything away.')
    }

    r_cut <- r[cut_idx]
    obs_cut <- curve_set[['obs']][cut_idx]
    sim_m_cut <- curve_set[['sim_m']][cut_idx, , drop = FALSE]
    theo <- curve_set[['theo']]
    n_theo <- length(theo)
    if (n_theo > 0L) {
        theo_cut <- theo[cut_idx]
    }

    res <- list(r=r_cut, obs=obs_cut, sim_m=sim_m_cut)
    if (n_theo > 0L) {
        res[['theo']] <- theo_cut
    }
    if(with(curve_set, exists('is_residual'))) res[['is_residual']] <- curve_set[['is_residual']]
    else res[['is_residual']] <- FALSE

    res <- create_curve_set(res)
    res
}
