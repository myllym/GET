# FIXME: inheritParams does not work in the expected way here.
#' Do a random labelling deviation test for a pattern using K_f- or
#' L_f-functions.
#'
#' @inheritParams marksummary::summ_func_random_labelling
#' @inheritParams crop_curves
#' @inheritParams residual
#' @inheritParams scale_curves
#' @inheritParams deviation
#' @inheritParams estimate_p_value.default
#' @export
random_labelling_test <- function(pattern,
                                  mtf_name = 'm', nsim = 999L,
                                  r_min = NULL, r_max = NULL, r_vec = NULL,
                                  measure = 'max', scaling = 'qdir',
                                  use_L = TRUE,
                                  edge_correction = 'translate',
                                  use_theo = TRUE, method = 'permute',
                                  use_biased_lambda2 = FALSE,
                                  ties = 'random') {
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
    curve_set <- scale_curves(curve_set, scaling = scaling)
    devs <- deviation(curve_set, measure = measure)
    p <- estimate_p_value(devs)
    p
}

#' Scale the curve_set and perform a deviation test.
#'
#' @inheritParams convert_envelope
#' @inheritParams crop_curves
#' @export
#' @examples
#' ## Testing complete spatial randomness (CSR)
#' #-------------------------------------------
#' require(spatstat)
#' pp <- unmark(spruces)
#' Generate nsim simulations under CSR, calculate L-function for the data and simulations
#' env <- envelope(pp, fun="Lest", nsim=999, savefuns=TRUE, correction="translate")
#' # The deviation test using the integral deviation measure
#' res <- deviation_test(env, measure = 'int')
#' res
#' # or
#' res <- deviation_test(env, r_min=0, r_max=7, measure='int')
deviation_test <- function(curve_set, r_min = NULL, r_max = NULL,
                                  use_theo = TRUE, scaling = 'qdir',
                                  measure = 'max', savedevs=FALSE) {
    curve_set <- crop_curves(curve_set, r_min = r_min, r_max = r_max)
    curve_set <- residual(curve_set, use_theo = use_theo)
    curve_set <- scale_curves(curve_set, scaling = scaling)
    devs <- deviation(curve_set, measure = measure)
    p <- estimate_p_value(devs)
    p
}
