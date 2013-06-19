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

#' Deviation test
#'
#' Crop the curve set to the interval of distances [r_min, r_max],
#' calculate residuals, scale the residuals and perform a deviation test
#' with a chosen deviation measure.
#'
#'
#' The deviation test is based on a test function T(r) and it works as follows:
#'
#' 1) The test function estimated for the data, T_1(r), and for nsim simulations
#' from the null model, T_2(r), ...., T_{nsim+1}(r), must be saved in 'curve_set'
#' and given to the deviation_test function.
#'
#' 2) The deviation_test function then
#'\itemize{
#'   \item Crops the functions to the chosen range of distances [r_min, r_max].
#'   \item Calculates the residuals d_i(r) = T_i(r) - T_0(r), where T_0(r) is
#'     the expectation of T(r) under the null hypothesis.
#'     If use_theo = TRUE, as T_0(r) the theoretical value given in the curve_set$theo
#'     is used, if it is given. Otherwise, T_0(r) is estimated by the mean of 
#'     T_j(r), j=2,...,nsim+1.
#'   \item Scales the residuals. Options are
#'         \itemize{
#'           \item 'none' No scaling. Nothing done.
#'           \item 'q' Quantile scaling.
#'           \item 'qdir' Directional quantile scaling.
#'           \item 'st' Studentised scaling.
#'         }
#'   \item Calculates the global deviation measure u_i, i=1,...,nsim+1, see options for 'measure'.
#'   \item Calculates the p-value.
#'}
#'
#' @inheritParams convert_envelope
#' @inheritParams crop_curves
#' @inheritParams residual
#' @inheritParams deviation
#' @param scaling The name of the scaling to use. Options include 'none',
#'   'q', 'qdir' and 'st'. 'qdir' is default.
#' @return If 'savedevs=FALSE' (default), the p-value is returned.
#' If 'savedevs=TRUE', then a list containing the p-value and calculated deviation measures
#' u_i, i=1,...,nsim+1 (where u_1 corresponds to the data pattern) is returned.
#' @references
#' Myllymäki, M., Grabarnik, P., Seijo, H. and Stoyan. D. (2013). Deviation test construction and power comparison for marked spatial point patterns. arXiv:1306.1028
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
