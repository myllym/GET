#' Random labeling with K_f functions
#'
#' Perform simulations under the random labelling for a marked point pattern
#' and estimate K_f- or L_f-functions for the pattern and simulations.
#'
#' @param pattern A \code{\link[spatstat]{ppp}} object as the simple marked
#'   point pattern to be analysed. The marks need to be in the form of a
#'   numeric vector. The window has to have the type "rectangle".
#' @param mtf_name A vector of mark test function names. "1" stands for the
#'   unmarked K-function. Accepted values are '1', 'm', 'mm', 'gamma',
#'   'gammaAbs' and 'morAbs'. See \code{\link[marksummary]{summ_func_random_labelling}}.
#' @param nsim The number of permutations.
#' @param r_min The minimum radius to include.
#' @param r_max A positive scalar value representing the maximum radius that
#'   should be considered. r_vec overrides r_max in calculation of functions.
#'   However, after calculation of functions, they are cropped to the interval
#'   [r_min, r_max] (if max(r_vec) > r_max loss in computing time).
#'   By default, r_max is NULL and will get a sensible default.
#' @param r_vec A monotonically increasing vector of non-negative r-values
#'   to act as the endpoints of the bins for the K_f-functions. The bins are
#'   exclusive on the left and inclusive on the right. If the first vector
#'   element has value zero, it will be regarded as the collapsed bin [0, 0],
#'   and the next bin will start from and exclude 0.
#' @param use_L A boolean describing whether the L_f function, L_f = sqrt(K_f/pi),
#'   should be used instead of the K_f function.
#' @param edge_correction The name of the edge correction to be used. Options
#'   are 'translate' and 'none'.
#' @param method The name of the method to create simulations under the null
#'   hypothesis. 'permute' results in permutations of the marks. Using
#'   'sample' will sample the marks from the empirical mark distribution
#'   with replacement. 'permute' is the default.
#' @param use_biased_lambda2 A logical scalar on whether to use the biased
#'   or the unbiased (in the Poisson case) estimate of the intensity
#'   squared.
#' @inheritParams residual
#' @seealso \code{\link{rank_envelope}}, \code{\link{st_envelope}}, \code{\link{qdir_envelope}}, \code{\link{deviation_test}}.
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
            r_vec = r_vec, nsim = nsim, calc_theo = use_theo,
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
#' @inheritParams deviation_test
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
#                      r_vec = r_vec, nsim = nsim, calc_theo = use_theo,
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
    curve_set <- scale_curves(curve_set, scaling = scaling)
    devs <- deviation(curve_set, measure = measure)
    p <- estimate_p_value(devs)
    p
}
