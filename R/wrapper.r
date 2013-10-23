#' Random labeling with K_f functions
#'
#' Perform simulations under the random labelling for a marked point pattern
#' and estimate K_f- or L_f-functions for the pattern and simulations.
#'
#' @inheritParams marksummary::summ_func_random_labelling
#' @inheritParams crop_curves
#' @inheritParams residual
#' @param use_L A boolean describing whether the L_f function, L_f = sqrt(K_f/pi),
#'   should be used instead of the K_f function.
#' @seealso \code{\link{rank_envelope}}, \code{\link{st_envelope}}, \code{\link{qdir_envelope}}, \code{\link{deviation_test}}
#' @references
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2013). Deviation test construction and power comparison for marked spatial point patterns. arXiv:1306.1028
#' @export
random_labelling <- function(pattern,
        mtf_name = 'm', nsim = 999L,
        r_min = NULL, r_max = NULL, r_vec = NULL,
        use_L = TRUE,
        edge_correction = 'translate',
        use_theo = TRUE, method = 'permute',
        use_biased_lambda2 = FALSE) {
    if (use_L) {
        func <- random_labelling_L_f
    } else {
        func <- random_labelling_K_f
    }

    curve_set <- func(pattern, mtf_name = mtf_name, r_max = r_max,
            r_vec = r_vec, n_sim = nsim, calc_theo = use_theo,
            edge_correction = edge_correction, method = method,
            use_biased_lambda2 = use_biased_lambda2)
    # r_max has already been dealt with in marksummary.
    curve_set <- crop_curves(curve_set, r_min = r_min, r_max = NULL)
    curve_set <- residual(curve_set, use_theo = use_theo)
    curve_set
}

# FIXME: inheritParams does not work in the expected way here.
#' Random labeling test
#'
#' Do a random labelling deviation test for a pattern using K_f- or
#' L_f-functions.
#'
#' Given a marked point pattern, the function performs simulations under
#' the random labelling hypothesis, estimates K_f- or L_f-functions
#' for the pattern and simulations and makes a deviation test.
#'
#' @inheritParams random_labelling
#' @inheritParams scale_curves
#' @inheritParams deviation
#' @inheritParams estimate_p_value.default
#' @export
#' @examples
#' require(spatstat)
#' mpp <- spruces
#' # T(r) = \hat{L}_m(r), an estimator of the L_m(r) function
#' p <- random_labelling_test(mpp, mtf_name = 'm', nsim=999, r_min=0, r_max=9.5)
#' p
random_labelling_test <- function(pattern,
                                  mtf_name = 'm', nsim = 999L,
                                  r_min = NULL, r_max = NULL, r_vec = NULL,
                                  measure = 'max', scaling = 'qdir',
                                  use_L = TRUE,
                                  edge_correction = 'translate',
                                  use_theo = TRUE, method = 'permute',
                                  use_biased_lambda2 = FALSE,
                                  ties = 'random') {
#    if (use_L) {
#        func <- random_labelling_L_f
#    } else {
#        func <- random_labelling_K_f
#    }
#
#    curve_set <- func(pattern, mtf_name = mtf_name, r_max = r_max,
#                      r_vec = r_vec, n_sim = nsim, calc_theo = use_theo,
#                      edge_correction = edge_correction, method = method,
#                      use_biased_lambda2 = use_biased_lambda2)
#    # r_max has already been dealt with in marksummary.
#    curve_set <- crop_curves(curve_set, r_min = r_min, r_max = NULL)
#    curve_set <- residual(curve_set, use_theo = use_theo)
#    curve_set <- scale_curves(curve_set, scaling = scaling)
    curve_set <- random_labelling(pattern=pattern,
            mtf_name = mtf_name, nsim = nsim,
            r_min = r_min, r_max = r_max, r_vec = r_vec,
            use_L = use_L,
            edge_correction = edge_correction,
            use_theo = use_theo, method = method,
            use_biased_lambda2 = use_biased_lambda2)
    devs <- deviation(curve_set, measure = measure)
    p <- estimate_p_value(devs)
    p
}
